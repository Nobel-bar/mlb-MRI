%==================================================================================================
% QSM RAWデータ (3D) を読み込み、
% 3D実空間で回転 (imrotate3) をシミュレートし、
% k空間でハイブリッド化 (一部置換) を行うMATLABスクリプト
%
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
Back_img =  exp(1i * (iFreq - RDF)) ./abs(exp(1i * (iFreq - RDF))); % 正規化した


fprintf('\n--- 2Dスライスごとの処理ループを開始します ---\n');

% [★修正] 最終的な3Dアーティファクト画像を (デバッグ用に) 保存する行列を初期化
artifact_img_3D = complex(zeros(params.matrix_size));

% [★修正] 3D処理から2Dスライスごとのループ処理に変更
%          元のセクション7, 8 をこのループに統合
for slice_idx = 1:num_slices
    
    fprintf('\n--- スライス %d / %d を処理中 ---\n', slice_idx, num_slices);

    % [★修正] 現在のスライスデータを抽出 (2D)
    original_img_slice = original_img(:,:,slice_idx);
    Highpass_img_slice = Highpass_img(:,:,slice_idx);
    Back_img_slice = Back_img(:,:,slice_idx);

%% --- 3. 実空間データの拡張 (ゼロパディング) ---
    fprintf('3. (2D) 実空間データのY方向を拡張しています...\n');

    % [★修正] 拡張行列を2Dで初期化 (512 x 788)
    extend_org = complex(zeros(matrix_x, magnification));
    extend_high = complex(zeros(matrix_x, magnification));
    extend_back = complex(zeros(matrix_x, magnification));

    % 拡張後Yサイズの中心を計算 (788/2 + 1 = 395)
    y_center_final_ext = floor(magnification / 2) + 1;
    % 元のYサイズ (512) を中央に配置するための開始・終了インデックスを計算
    y_start_final = y_center_final_ext - floor(matrix_y / 2);
    y_end_final = y_start_final + matrix_y - 1;

    % [★修正] 2Dスライスをゼロ行列の中央に代入
    extend_org(:,y_start_final:y_end_final) = original_img_slice;
    extend_high(:,y_start_final:y_end_final) = Highpass_img_slice; 
    extend_back(:,y_start_final:y_end_final) = Back_img_slice; 

    % fprintf('   2D実空間の準備が完了しました。\n');

%% --- 4. 2D 実空間の回転（モーションのシミュレート） ---
    fprintf('4. (2D) 実空間の回転 (imrotate) を実行中...\n');

    % (4.1. original_img の回転はコメントアウトのまま)
    
    % --- 4.2. Highpass_img (extend_high) の回転 ---
    large_high_Re = real(extend_high);
    large_high_Im = imag(extend_high);
    
    % [★修正] imrotate3 (3D) から imrotate (2D) に変更
    %          rotation_axis は 2D 回転では不要
    % fprintf('     実数部 (Highpass) を回転中...\n');
    rotated_high_Re = imrotate(large_high_Re, theta, 'linear', 'crop');
    % fprintf('     虚数部 (Highpass) を回転中...\n');
    rotated_high_Im = imrotate(large_high_Im, theta, 'linear', 'crop');
    extend_high_rotated = complex(rotated_high_Re, rotated_high_Im);
    clear large_high_Re large_high_Im rotated_high_Re rotated_high_Im;
    
    % (4.3. Back_img の回転はコメントアウトのまま)
    % fprintf('   2D 実空間の回転が完了しました。\n'); 

    % --- 4.4. 回転後実空間データの合成 ---
    % [★ 指示通り、ロジックは変更しない]
    % 「回転後Highpass」 * 「回転前Background(実空間)」
    space_rotated = extend_high_rotated .* extend_back; % (512x788)


%% --- 5. 2D k空間データのハイブリッド化 ---
    fprintf('5. (2D) k空間のハイブリッド化を実行中...\n');

    % [★修正] fftn (3D) から fft2 (2D) に変更
    % fprintf('     回転前のk空間 (hybrid_space) を作成中 (fft2)...\n');
    hybrid_space = fftshift(fft2(extend_org)); % 回転前 (512 x 788)

    % fprintf('     回転後のk空間 (k_space_rolate) を作成中 (fft2)...\n');
    k_space_rolate = fftshift(fft2(space_rotated)); % 回転後 (512 x 788)

    % (x_... y_... のインデックス計算はそのまま流用)
    x_center = floor(matrix_x / 2) + 1; % 257
    x_start_cut = x_center - floor(cutted_matrix_x / 2); % 145
    x_end_cut = x_start_cut + cutted_matrix_x - 1; % 368
    y_center_org = floor(magnification / 2) + 1; % 395
    y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2); % 219
    y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1; % 570
    hybrid_row_indices = (x_start_cut + pix_start_row - 1) : (x_start_cut + pix_start_row + width - 2); 
    
    % fprintf('     k空間データをハイブリッド化 (置換) しています...\n');
    
    % [★修正] 3Dインデックス (:,:,:) から 2Dインデックス (:,:) に変更
    hybrid_space(hybrid_row_indices, :) = k_space_rolate(hybrid_row_indices, :);

    clear k_space_rolate; % メモリ節約

    % --- k空間のROI切り出し ---
    % fprintf('     k空間のROIを切り出しています...\n');
    % [★修正] 3D (..., num_slices) から 2D に変更
    final_k_space = complex(zeros(matrix_x, magnification));

    % [★修正] 3Dインデックス (:,:,:) から 2Dインデックス (:,:) に変更
    final_k_space(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut) = ...
        hybrid_space(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut);
    
    clear hybrid_space;

%% --- 6. アーティファクト画像の再構成 ---
    fprintf('6. (2D) アーティファクト画像を再構成中 (ifft2)...\n');
    
    % [★修正] ifftn (3D) から ifft2 (2D) に変更
    artifact_img_ext = ifft2(ifftshift(final_k_space));

    clear final_k_space; 

    % [★修正] 3Dスライス (...,:) から 2Dスライス (:,:) に変更
    % [★修正] 変数名を artifact_img_slice に変更
    artifact_img_slice = artifact_img_ext(:, y_start_final:y_end_final);

    clear artifact_img_ext; 
    
    % [★修正] 3Dボリュームに結果を格納 (デバッグ・確認用)
    artifact_img_3D(:,:,slice_idx) = artifact_img_slice;


%% --- 7. スライスごとの表示 ---
    % [★修正] ループに統合。デバッグ用にスライス12のみ表示
    if slice_idx == 12
        fprintf('7. *** スライス 12 を検知、画像を表示します ***\n');
        
        % 元の強度画像
        original_img_slice_mag = iMag_4D(:,:,slice_idx, echo_idx);

        figure('Name', sprintf('Slice %d Comparison (2D Process)', slice_idx), 'WindowState', 'maximized');
        subplot(1, 2, 1);
        imshow(original_img_slice_mag, []);
        title(sprintf('Original (Slice %d)', slice_idx));
        
        subplot(1, 2, 2);
        imshow(abs(artifact_img_slice), []);
        title(sprintf('Artifact (Slice %d, 2D Process)', slice_idx));

        drawnow;
    else
        % fprintf('7. スライス %d の表示はスキップします。\n', slice_idx);
    end
    
%% --- 8. 結果のファイル保存 (ループ内に移動) ---
    fprintf('8. スライス %d の結果をファイルに保存中...\n', slice_idx);
    
    final_artifact_img_slice = artifact_img_slice; % 2Dなので転置は不要 (元コードの仮定に従う)
        
    % [★修正] ファイル名を 2D 処理用に変更
    filename_base = sprintf('slice%03d_2D_rotate_artifact_th%.1f', slice_idx, theta);
    save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(final_artifact_img_slice));
    save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(final_artifact_img_slice));
    save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(final_artifact_img_slice));
    save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(final_artifact_img_slice));

end % --- スライスループ (for slice_idx) の終了 ---

fprintf('\n...全 %d スライスのシミュレーションが完了しました。\n', num_slices);
fprintf('結果は %s に保存されました。\n', save_path);

% -------------------------------------------------------------------
% スクリプトの最後にローカル関数を定義します
% (save_raw_data 関数は変更なし)
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