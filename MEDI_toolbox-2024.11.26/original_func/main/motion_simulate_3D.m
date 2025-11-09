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
Back_img =  exp(1i * (iFreq - RDF)); % 正規化した

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
rotated_org_Im = imrotate3(large_org_Im, theta, rotation_axis, 'linear', 'crop');
extend_org_rotated = complex(rotated_org_Re, rotated_org_Im);
clear large_org_Re large_org_Im rotated_org_Re rotated_org_Im;

% --- 4.2. Highpass_img (extend_high) の回転 ---
large_high_Re = real(extend_high);
large_high_Im = imag(extend_high);
fprintf('    実数部 (Highpass) を回転中...\n');
rotated_high_Re = imrotate3(large_high_Re, theta, rotation_axis, 'linear', 'crop');
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
fprintf('5. 3D k空間のハイブリッド化を実行中...\n');

% [★修正: `large_org` ではなく `extend_org` (回転前実空間) をk空間へ]
fprintf('    回転前のk空間 (hybrid_space) を作成中 (fftn)...\n');
hybrid_space = fftshift(fftn(extend_org)); % 回転前 (512 x 788 x 23)

direct_hybrid_space = fftshift(fftn(extend_org_rotated)); % 回転前 (512 x 788 x 23)

% [★修正: `space_rotated_base` (回転後実空間) をk空間へ]
fprintf('    回転後のk空間 (k_space_rolate) を作成中 (fftn)...\n');
k_space_rolate = fftshift(fftn(space_rotated)); % 回転後 (512 x 788 x 23)

direct_k_space_rolate = fftshift(fftn(space_rotated)); % 回転後 (512 x 788 x 23)
x_center = floor(matrix_x / 2) + 1; % 257
x_start_cut = x_center - floor(cutted_matrix_x / 2); % 257 - 112 = 145
x_end_cut = x_start_cut + cutted_matrix_x - 1; % 145 + 224 - 1 = 368

% [★修正: `orig_matrix_y` -> `magnification`]
y_center_org = floor(magnification / 2) + 1; % 788 -> 395
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2); % 395 - 176 = 219
y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1; % 219 + 352 - 1 = 570


% ハイブリッド化を行う行の絶対インデックス
% (145 + 112 - 1) = 256 から 112 行分
hybrid_row_indices = (x_start_cut + pix_start_row - 1) : (x_start_cut + pix_start_row + width - 2); 

 
fprintf('    k空間データをハイブリッド化 (置換) しています...\n');
% 指定された行 (256:367) を、回転後のk空間データで上書き
hybrid_space(hybrid_row_indices, :, :) = k_space_rolate(hybrid_row_indices, :, :);
direct_hybrid_space(hybrid_row_indices, :, :) = direct_k_space_rolate(hybrid_row_indices, :, :);

clear k_space_rolate; % メモリ節約

% --- k空間のROI切り出し (ゼロパディング解除に相当) ---
fprintf('    k空間のROIを切り出しています...\n');
% 最終k空間 (512 x 788 x 23) をゼロで初期化
final_k_space = complex(zeros(matrix_x, magnification, num_slices));

direct_final_k_space = complex(zeros(matrix_x, magnification, num_slices));
% [★修正: `y_start_cut` -> `y_start_org_cut`]
% `hybrid_space` から k空間ROI (224 x 352) を切り出し、
% `final_k_space` の同じ位置に配置　実際に読み取っている範囲
% --- k空間ROIの絶対インデックスを計算 ---
% [★修正: `orig_matrix_x` -> `matrix_x`])

% 
% final_k_space(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut, :) = ...
%     hybrid_space(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut, :);
final_k_space=hybrid_space;
 direct_final_k_space(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut, :) = ...
    direct_hybrid_space(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut, :);
clear hybrid_space;
 
%% --- 6. アーティファクト画像の再構成 ---
fprintf('6. アーティファクト画像を再構成中 (ifftn)...\n');
[max_val, max_idx] = max(abs(final_k_space(:)));
[kk, mm, nn] = ind2sub(size(final_k_space), max_idx);
fprintf('k空間の最大値は座標 (%d, %d, %d) にあります。\n', kk, mm, nn);
p0_factor = final_k_space(max_idx) / max_val;
k_space_p0 = final_k_space / p0_factor;

artifact_img_ext = ifftn(ifftshift(final_k_space));
direct_img_ext = ifftn(ifftshift(direct_final_k_space));
clear final_k_space; 

% 拡張したY次元 (788) を元のY次元 (512) に戻す (中央を切り出す)
artifact_img = artifact_img_ext(:, y_start_final:y_end_final, :);
direct_img = direct_img_ext(:, y_start_final:y_end_final, :);
clear artifact_img_ext; 


%% --- 7. スライスごとの表示と保存 ---
fprintf('7. スライスごとの表示と保存を開始します...\n');

% [★修正: `Slice` -> `num_slices`]
for slice_idx = 12:12 % (デバッグのためスライス12のみ)
    
    fprintf('   --- スライス %d / %d を処理中 ---\n', slice_idx, num_slices);

    % 現在のスライスを3Dボリュームから抽出
    % [★修正: `Mag_4D` -> `iMag_4D`]
    djrect_img_slice = direct_img(:,:,slice_idx, echo_idx); % 強度画像を表示
    artifact_img_slice = artifact_img(:,:,slice_idx);

    % (マスク適用は省略)

    % 各スライスごとに新しいFigureを作成する
    figure('Name', sprintf('Slice %d Comparison', slice_idx), 'WindowState', 'maximized');
     
    % 1行2列のグリッドの1番目（左側）
    subplot(1, 2, 1);
    imshow(abs(djrect_img_slice), []);
    title(sprintf('Original (Slice %d)', slice_idx));
     
    % 1行2列のグリッドの2番目（右側）
    subplot(1, 2, 2);
    imshow(abs(artifact_img_slice), []);
    title(sprintf('Artifact (Slice %d)', slice_idx));

    % Figureを強制的に今すぐ描画する
    drawnow;
     
    %% --- 8. (修正) 結果のファイル保存 (ループ内に移動) ---
     
    % 最終画像を permute (転置) する
    % (注: 以前の 2D スクリプトでは転置していた)
    % final_artifact_img_slice = permute(artifact_img_slice, [2 1]);
    final_artifact_img_slice = artifact_img_slice; % 3D処理では転置不要と仮定
     


end % --- スライスループ (for slice_idx) の終了 ---
% ファイル名にスライス番号を含める
filename_base = sprintf('3D_rotate_artifact_th%.1f', theta);
save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(artifact_img));

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

