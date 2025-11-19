%==================================================================================================
% QSM RAWデータ (3D) を読み込み、
% 3D実空間で回転 (imrotate3) をシミュレートし、
% k空間でハイブリッド化 (一部置換) を行うMATLABスクリプト
%
% [修正内容 v2]
% - user_21 のコメント指示に基づき、extend_high/back の代入を修正
% - user_21 のコメント指示に基づき、space_rotated_base の切り出し処理を追加 (推測)
% - 未定義変数 (rotation_axis, orig_matrix_x, orig_matrix_y, Mag_4D, Slice) を修正
% - 致命的なエラー (clear ;) を修正
%==================================================================================================

fprintf('スクリプトを開始します (3D実空間回転 + k空間ハイブリッド)\n');
clear variables;
close all;

%% --- 1. 撮像・シミュレーション パラメータ設定 ---
fprintf('1. パラメータを設定しています...\n');

% パス設定
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

image_file_0 = image_file_00; % slab用

% 読み込みパスと保存パスを定義
load_base_path = fullfile(image_file_0, image_file_1);
load_mask_path = fullfile(image_file_0, image_file_3);
save_path = fullfile(image_file_0, image_file_4);

% 保存先フォルダが存在しない場合は作成する
if ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf('保存フォルダを作成しました: %s\n', save_path);
end

% --- params構造体: QSMデータの撮像パラメータ ---
params = struct();
params.matrix_size = [512, 512, 23]; % !! 要変更 !! : 行列サイズ
params.TE = [0.015]; % !! 要変更 !! : エコー時間 (秒)
% (他の params は現在使用されていません)

% --- 拡張と回転のパラメータ ---
extention = 2.0/1.3; % Y方向の拡張率 (約1.54倍)
magnification = round( params.matrix_size(2) * extention); % 拡張後のYサイズ (512 -> 788)
theta = -18.6; % 回転角度 (度)
rotation_axis = [0 0 1]; % 回転軸 (Z軸) [★修正: 未定義変数を定義]
IdealSize = 512; % 画像処理の基準サイズ (現在未使用)

% --- k空間ハイブリッド化パラメータ ---
cutted_matrix_x = 224; % 実際に収集されたk空間の有効データサイズ
cutted_matrix_y = 352;
width = 112; % 置き換える行数
pix_start_row = 112; % k空間ROIの何行目から置き換えるか (1-based)

%% --- 2. RAWデータの読み込み ---
fprintf('2. RAWデータと処理済みデータを読み込んでいます...\n');

% ファイル名
mag_filename = '1st_2DGE_0deg_mag.raw'; % !! 要変更 !!
phase_filename = '1st_2DGE_0deg_phase.raw'; % !! 要変更 !!

mag_filepath = fullfile(load_base_path, mag_filename);
phase_filepath = fullfile(load_base_path, phase_filename);

% データの次元を定義 [x, y, z, echo]
dims = [params.matrix_size, length(params.TE)];
num_elements = prod(dims); 
precision = 'double=>double'; 

% --- 強度・位相データの読み込み ---
fid_mag = fopen(mag_filepath, 'rb');
if fid_mag == -1, error('強度ファイルを開けませんでした: %s', mag_filepath); end
iMag_4D = reshape(fread(fid_mag, inf, precision), dims);
fclose(fid_mag);

fid_phase = fopen(phase_filepath, 'rb');
if fid_phase == -1, error('位相ファイルを開けませんでした: %s', phase_filepath); end
iPhase_4D = reshape(fread(fid_phase, inf, precision), dims);
fclose(fid_phase);

% --- 処理済みデータ (iFreq, RDF) の読み込み ---
try
    load(fullfile(load_mask_path, 'phase.mat'), 'iFreq');
    load(fullfile(load_mask_path, 'PDF.mat'), 'RDF');
    % load(fullfile(load_mask_path, 'Mask.mat'), 'Mask');
catch ME
    fprintf('phase.mat または PDF.mat の読み込みに失敗しました。\n');
    rethrow(ME);
end

% --- 変数定義 ---
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);
echo_idx = 1; % 最初のエコーのみ使用

fprintf('データの読み込み完了。%d スライス、%d エコーのデータを処理します。\n', num_slices);
original_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * iPhase_4D(:,:,:, echo_idx));

% Highpass_img = original_img ./ fitting;  % fittingを正規化したうえで
Highpass_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * (RDF));
% Back_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * (iFreq - RDF));
Back_img =  exp(1i * (iFreq - RDF)); % ./abs(exp(1i * (iFreq - RDF))); % 正規化した

%% --- 3. 実空間データの拡張 (ゼロパディング) ---
fprintf('3. 実空間データのY方向を拡張 (ゼロパディング) しています...\n');

% 拡張後のサイズの複素数ゼロ行列を初期化 (512 x 788 x 23)
extend_org = complex(zeros(matrix_x, magnification, num_slices));
extend_high = complex(zeros(matrix_x, magnification, num_slices));
extend_back = complex(zeros(matrix_x, magnification, num_slices));

% 拡張後Yサイズの中心を計算 (788/2 + 1 = 395)
y_center_final_ext = floor(magnification / 2) + 1;
% 元のYサイズ (512) を中央に配置するための開始・終了インデックスを計算
% y_start_final = 395 - 256 = 139
y_start_final = y_center_final_ext - floor(matrix_y / 2);
% y_end_final = 139 + 512 - 1 = 650
y_end_final = y_start_final + matrix_y - 1;

% [★修正: user_21 コメント反映] ゼロ行列の中央にデータを代入
extend_org(:,y_start_final:y_end_final, :) = original_img;
extend_high(:,y_start_final:y_end_final, :) = Highpass_img; 
extend_back(:,y_start_final:y_end_final, :) = Back_img; 

fprintf('   3D実空間の準備が完了しました。\n');

%% --- 4. 3D 実空間の回転（モーションのシミュレート） ---
fprintf('4. 3D 実空間の回転 (imrotate3) を実行中...\n');

% --- 4.1. original_img (extend_org) の回転 ---
large_org_Re = real(extend_org);
large_org_Im = imag(extend_org);
fprintf('    実数部 (original) を回転中...\n');
rotated_org_Re = imrotate3(large_org_Re, theta, rotation_axis, 'linear', 'crop');
fprintf('    虚数部 (original) を回転中...\n');
rotated_org_Im = imrotate3(large_org_Im, theta, rotation_axis, 'linear', 'crop');
extend_org_rotated = complex(rotated_org_Re, rotated_org_Im);
% [★修正: 致命的エラー `clear ;` を修正]
clear large_org_Re large_org_Im rotated_org_Re rotated_org_Im;

% --- 4.2. Highpass_img (extend_high) の回転 ---
large_high_Re = real(extend_high);
large_high_Im = imag(extend_high);
fprintf('    実数部 (Highpass) を回転中...\n');
rotated_high_Re = imrotate3(large_high_Re, theta, rotation_axis, 'linear', 'crop');
fprintf('    虚数部 (Highpass) を回転中...\n');
rotated_high_Im = imrotate3(large_high_Im, theta, rotation_axis, 'linear', 'crop');
extend_high_rotated = complex(rotated_high_Re, rotated_high_Im);
clear large_high_Re large_high_Im rotated_high_Re rotated_high_Im;
% 
% % --- 4.3. Back_img (extend_back) の回転 ---
% large_back_Re = real(extend_back);
% large_back_Im = imag(extend_back);
% fprintf('    実数部 (Background) を回転中...\n');
% rotated_back_Re = imrotate3(large_back_Re, theta, rotation_axis, 'linear', 'crop');
% fprintf('    虚数部 (Background) を回転中...\n');
% rotated_back_Im = imrotate3(large_back_Im, theta, rotation_axis, 'linear', 'crop');
% extend_back_rotated = complex(rotated_back_Re, rotated_back_Im);
% clear large_back_Re large_back_Im rotated_back_Re rotated_back_Im;
% fprintf('   3D 実空間の回転が完了しました。\n'); 

% --- 4.4. 回転後実空間データの合成 ---
% 「回転後Highpass」 * 「回転前Background(実空間)」を意図していた
% ただし extend_back，extend_high_rotated は (512x788x23)
space_rotated = extend_high_rotated .* extend_back;


%% --- 5. 3D k空間データのハイブリッド化 ---
fprintf('5. (2D) k空間のハイブリッド化を実行中...\n'); % [★修正] fprintfのコメントアウトを解除
artifact_img_3D = complex(zeros(params.matrix_size));

for slice_idx = 1:num_slices  
    
    fprintf('   --- スライス %d / %d を処理中 ---\n', slice_idx, num_slices);

    % [★修正: ロジックエラー]
    %   回転「前」の extend_org をベースk空間の元データにする
    % org_slice = extend_org_rotated(:,:,slice_idx); % [★誤り]
    org_slice = extend_org(:,:,slice_idx); % [★正解] 回転前の画像 (extend_org) を使用
    direct_slice = extend_org_rotated(:,:,slice_idx);
    rotated_slice = space_rotated(:,:,slice_idx);
    
    % [★修正: 変数名を明確化]
    k_space_org = fftshift(fft2(org_slice)); % 回転前k空間 (512 x 788)
    k_space_rotated = fftshift(fft2(rotated_slice)); % 回転後k空間 (512 x 788)
    direct_k_space_org = fftshift(fft2(direct_slice)); % 回転語k空間 (512 x 788)
    
    % [★修正: k_space_org をコピーして hybrid_k_space を作成]
    hybrid_k_space = k_space_org; 
    
    x_center = floor(matrix_x / 2) + 1; % 257
    x_start_cut = x_center - floor(cutted_matrix_x / 2); % 145
    x_end_cut = x_start_cut + cutted_matrix_x - 1; % 368
    y_center_org = floor(magnification / 2) + 1; % 395
    y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2); % 219
    y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1; % 570
    hybrid_row_indices = (x_start_cut + pix_start_row - 1) : (x_start_cut + pix_start_row + width - 2); 
    
    % [★修正: コピー (hybrid_k_space) に対してハイブリッド化を実行]
    hybrid_k_space(hybrid_row_indices, :) = k_space_rotated(hybrid_row_indices, :);
    
    
    % --- k空間のROI切り出し ---
    final_k_space = complex(zeros(matrix_x, magnification));
    direct_final_k_space = complex(zeros(matrix_x, magnification));

    % [★修正: final_k_space には「ハイブリッド化済み」の k_space_hybrid を入れる]
    final_k_space(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut) = ...
        hybrid_k_space(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut);
    
    % [★修正: direct_final_k_space には「ハイブリッド化前」の k_space_org を入れる]
    direct_final_k_space(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut) = ...
        direct_k_space_org(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut);
    
    
    % [★修正: タイポと不要な変数のクリア]
    % clear hybrid_space hybrid_k_space; % [★誤り]
    clear k_space_org k_space_rotated hybrid_k_space; % [★正解]

%% --- 6. アーティファクト画像の再構成 ---
    fprintf('6. (2D) アーティファクト画像を再構成中 (ifft2)...\n');
    
    artifact_img_ext = ifft2(ifftshift(final_k_space));
    direct_img_ext = ifft2(ifftshift(direct_final_k_space));
    
    clear final_k_space direct_final_k_space; 
    
    artifact_img_slice = artifact_img_ext(:, y_start_final:y_end_final);
    direct_img_slice = direct_img_ext(:, y_start_final:y_end_final);
    
    clear artifact_img_ext direct_img_ext; 
    
    % [★修正] 3Dボリュームに結果を格納

    direct_img_slice = permute(direct_img_slice, [2 1 3]);
    artifact_img_slice = permute(artifact_img_slice, [2 1 3]);
    artifact_img_3D(:,:,slice_idx) = artifact_img_slice;



%% --- 7. スライスごとの表示 ---
    % [★修正] ループに統合。デバッグ用にスライス12のみ表示
    if slice_idx == 12
        fprintf('7. *** スライス 12 を検知、画像を表示します ***\n');
       
        original_img_slice_mag = permute(iMag_4D(:,:,slice_idx, echo_idx), [2 1 3]);
        % [★コメント] 修正後のコードでは、2列目と3列目は異なる画像になるはずです
        figure('Name', sprintf('Slice %d Comparison (2D Process)', slice_idx), 'WindowState', 'maximized');
        subplot(1, 3, 1);
        imshow(original_img_slice_mag, []);
        title(sprintf('Original (Slice %d)', slice_idx));

        subplot(1, 3, 2); 
        imshow(abs(artifact_img_slice), []);
        title(sprintf('Artifact (Slice %d)', slice_idx));

        % [★コメント] "Direct" はハイブリッド化「前」の画像 (＝Originalの複素画像)
        subplot(1, 3, 3);
        imshow(abs(direct_img_slice), []); 
        title(sprintf('Direct (non-hybrid) (Slice %d)', slice_idx));
        
        drawnow;
    else
        % fprintf('7. スライス %d の表示はスキップします。\n', slice_idx);
    end
   
end % --- スライスループ (for slice_idx) の終了 ---




%% --- 8. 保存 ---
artifact_img_3D = permute(artifact_img_3D, [2 1 3]);
     
    % ファイル名にスライス番号を含める
filename_base = sprintf('3d_2D_rotate_artifact_th%.1f', theta);
save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(artifact_img_3D));
save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(artifact_img_3D));
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(artifact_img_3D));
save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(artifact_img_3D));

fprintf('...スライス %d のシミュレーションが完了しました。\n', slice_idx);
fprintf('結果は %s に保存されました。\n', save_path);

% -------------------------------------------------------------------
% スクリプトの最後にローカル関数を定義します
% -------------------------------------------------------------------
function save_raw_data(filepath, data)
    % データを 'double' (計算時の型) で保存
    fid = fopen(filepath, 'w');
    if fid == -1
        error('ファイルが開けませんでした: %s', filepath);
    end
    fwrite(fid, data, 'double'); % 'single' ではなく 'double' を推奨
    fclose(fid);
end

