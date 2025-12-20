%==================================================================================================
% QSM RAWデータ (3D) シミュレーション
% 【改変版: 3段階 階段状ROI置換 (指定座標で角度切り替え)】
%==================================================================================================

fprintf('スクリプトを開始します (3段階回転切り替え版)\n');
clear variables;
close all;

% ★★★ 計算負荷設定 ★★★
DS_FACTOR = 2; 
fprintf('ダウンサンプリング係数: %d\n', DS_FACTOR);

%% --- 1. パラメータ設定 ---

%% --- 1. パラメータ設定 (パス等は環境に合わせてください) ---
image_file_2DGE_0deg_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_2DGE_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
correct_mag_filename = '1st_2DGE_1_2_Rotate_mag.raw';
correct_phase_filename = '1st_2DGE_1_2_Rotate_phase.raw';

mag_filename = '1st_2DGE_0deg_mag.raw';
phase_filename = '1st_2DGE_0deg_phase.raw';
% image_file_2DGE_0deg_H = image_file_2DGE_Rotate_H;
% mag_filename = '2DGE_Rotate_H_mag.raw';
% phase_filename = '2DGE_Rotate_H_phase.raw';

% --- ★変更点: 回転切り替え用パラメータ ---
alpha = 3;       % 分割数
beta  = 352;     % 全列数 (cutted_matrix_y に相当)
gamma = 100;     % 開始列 (pix_start_col に相当)

pix_start_row = 116; % 開始行 (User指定: 116行目から)
pix_start_col = gamma; 

% --- 設定: 比較対象のスライス番号 ---
target_slice = 20;

% --- パラメータ ---
params = struct();
params.original_matrix_size = [512, 512, 23];
extention = 2.0/1.3;
theta = -18.6;
rotation_axis = [0 0 1];

% 保存
filename_base = sprintf('Artifact_th%.1f_StepShape_Alpha%d_Beta%d_Gamma%d', theta, alpha, beta, gamma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_base_path = fullfile(image_file_2DGE_0deg_H, image_file_1);
load_mask_path = fullfile(image_file_2DGE_0deg_H, image_file_3);
correct_load_base_path = fullfile(image_file_2DGE_1_2_Rotate_H, image_file_1);
correct_load_mask_path = fullfile(image_file_2DGE_1_2_Rotate_H, image_file_3);
save_path = fullfile(image_file_0, image_file_4);
if ~exist(save_path, 'dir'), mkdir(save_path); end

% --- ハイブリッド化パラメータ ---
cutted_matrix_x = round(224 / DS_FACTOR);
cutted_matrix_y = round(beta / DS_FACTOR); % betaを使用

% ダウンサンプリング対応
pix_start_row = round(pix_start_row / DS_FACTOR);
pix_start_col = round(pix_start_col / DS_FACTOR);
gamma = round(gamma / DS_FACTOR);

% --- 体動パラメータ ---
MOTION_PE_ERROR_MAX_RAD = 0; 

%% --- 2. データ読み込み & ダウンサンプリング ---
fprintf('\n2. データ読み込み中...\n');
% (読み込み部分は変更なし)
dims_orig = params.original_matrix_size;
precision = 'double=>double';

fid = fopen(fullfile(load_base_path, mag_filename), 'rb');
data_vector_mag = fread(fid, inf, precision); fclose(fid);
params.original_matrix_size(3) = numel(data_vector_mag) / (dims_orig(1)*dims_orig(2));
dims_orig = params.original_matrix_size;
iMag_4D_orig = reshape(data_vector_mag, dims_orig);

fid = fopen(fullfile(load_base_path, phase_filename), 'rb');
iPhase_4D_orig = reshape(fread(fid, inf, precision), dims_orig); fclose(fid);

try
    load(fullfile(load_mask_path, 'phase.mat'), 'iFreq');
    load(fullfile(load_mask_path, 'PDF.mat'), 'RDF');
catch
    iFreq = zeros(dims_orig); RDF = zeros(dims_orig);
end

iMag_4D = iMag_4D_orig(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
iPhase_4D = iPhase_4D_orig(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
RDF_small = RDF(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
iFreq_small = iFreq(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
clear iMag_4D_orig iPhase_4D_orig RDF iFreq;



correct_fid = fopen(fullfile(correct_load_base_path, correct_mag_filename), 'rb');
correct_data_vector_mag = fread(correct_fid, inf, precision); fclose(correct_fid);
params.original_matrix_size(3) = numel(correct_data_vector_mag) / (dims_orig(1)*dims_orig(2));
dims_orig = params.original_matrix_size;
correct_iMag_4D_orig = reshape(correct_data_vector_mag, dims_orig);

correct_fid = fopen(fullfile(correct_load_base_path, correct_phase_filename), 'rb');
correct_iPhase_4D_orig = reshape(fread(correct_fid, inf, precision), dims_orig); fclose(correct_fid);

try
    load(fullfile(correct_load_mask_path, 'phase.mat'), 'iFreq');
    load(fullfile(correct_load_mask_path, 'PDF.mat'), 'RDF');
catch
    iFreq = zeros(dims_orig); RDF = zeros(dims_orig);
end

correct_iMag_4D = correct_iMag_4D_orig(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
correct_iPhase_4D = correct_iPhase_4D_orig(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
correct_RDF_small = RDF(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
correct_iFreq_small = iFreq(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
clear correct_iMag_4D_orig correct_iPhase_4D_orig correct_RDF correct_iFreq;

params.matrix_size = size(iMag_4D);
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);


correct_original_img = correct_iMag_4D .* exp(1i * correct_iPhase_4D);
original_img = iMag_4D .* exp(1i * iPhase_4D);



%% --- 3. 拡張 & 4. 回転 (スライスループ内で実施に変更) ---
fprintf('\n3-4. 実空間拡張とパラメータ準備...\n');
magnification = round(matrix_y * extention);
y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(matrix_y / 2);
y_end_final = y_start_final + matrix_y - 1;

extend_org = complex(zeros(matrix_x, magnification, num_slices));
extend_org(:,y_start_final:y_end_final, :) = original_img;
correct_extend_org = complex(zeros(matrix_x, magnification, num_slices));
correct_extend_org(:,y_start_final:y_end_final, :) = correct_original_img;

% ★重要: 3つの回転角度を定義
angles = [theta / alpha, (2 * theta) / alpha, (3 * theta) / alpha];
fprintf('回転角度設定: %.2f度, %.2f度, %.2f度\n', angles(1), angles(2), angles(3));


%% --- 5. k空間シミュレーション (多段階回転対応) ---
fprintf('\n5. k空間シミュレーション (複雑な角度切り替え処理)...\n');
tic; 

x_center = floor(matrix_x / 2) + 1;
x_start_cut = x_center - floor(cutted_matrix_x / 2);
x_end_cut = x_start_cut + cutted_matrix_x - 1;
y_center_org = floor(magnification / 2) + 1;
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2);
y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1;

kx_collect_indices = x_start_cut : x_end_cut;
ky_collect_indices = y_start_org_cut : y_end_org_cut;

% 結果格納用配列 (ここに計算結果を集約します)
k_space_artifact = complex(zeros(matrix_x, magnification, num_slices));



%% --- 7. 正解データとの比較 (12スライス目) ---
fprintf('\n7. 正解データとの比較画像を生成中...\n');
% --- Log変換用無名関数 ---
% log(abs(k) + 1) を計算することで、0除算を防ぎつつダイナミックレンジを圧縮
get_log_mag = @(k) log(abs(k) + 1);



% --- B. 比較用画像の生成 (回転なし, th/3, 2th/3, th) ---
% シミュレーション空間(拡張領域)での画像を生成し、最終出力サイズ(y_start_final:y_end_final)にクロップします

% ベース画像(拡張済み複素数データ)の12スライス目
Base_Complex_Slice = extend_org(:,:,target_slice);

% 正解データのk空間 (拡張領域全体) -> Log変換
Ref_Slice_K_Raw = fftshift(fft2(correct_extend_org(:,:,target_slice)));
Ref_Slice_K_Log = get_log_mag(Ref_Slice_K_Raw);

% 回転関数 (複素数対応)
rot_func_crop = @(img, ang) complex(...
    imrotate(real(img), ang, 'bilinear', 'crop'), ...
    imrotate(imag(img), ang, 'bilinear', 'crop'));

% 4つの状態のk空間を生成 -> Log変換
% 1. 回転なし (0度)
Img_0_ext = rot_func_crop(Base_Complex_Slice, 0);
Img_0_K_Log = get_log_mag(fftshift(fft2(Img_0_ext)));

% 2. θ/3 回転
Img_1_ext = rot_func_crop(Base_Complex_Slice, theta/3);
Img_1_K_Log = get_log_mag(fftshift(fft2(Img_1_ext)));

% 3. 2θ/3 回転
Img_2_ext = rot_func_crop(Base_Complex_Slice, (2*theta)/3);
Img_2_K_Log = get_log_mag(fftshift(fft2(Img_2_ext)));

% 4. θ 回転
Img_3_ext = rot_func_crop(Base_Complex_Slice, theta);
Img_3_K_Log = get_log_mag(fftshift(fft2(Img_3_ext)));

% --- C. 差分画像の計算 (正解データとの差) ---
% ※「どれくらい一致しているか」を見るため、絶対差分 (差の絶対値) を計算します。
% 黒いほど一致、白いほど不一致です。

% --- 差分計算 (Log Scaleでの差) ---
% 黒(0)に近いほど「Logスケールでのスペクトル」が一致している
Diff_0 = abs(Ref_Slice_K_Log - Img_0_K_Log);
Diff_1 = abs(Ref_Slice_K_Log - Img_1_K_Log);
Diff_2 = abs(Ref_Slice_K_Log - Img_2_K_Log);
Diff_3 = abs(Ref_Slice_K_Log - Img_3_K_Log);

% スケーリング決定 (最大誤差基準)
max_diff = max([max(Diff_0(:)), max(Diff_1(:)), max(Diff_2(:)), max(Diff_3(:))]);
if max_diff == 0, max_diff = 1; end 

% --- Figure表示 (k-space difference in Log Scale) ---
figure('Name', 'Log K-Space Difference from Correct Data', 'WindowState', 'maximized');
colormap(gray(256));

% 1
subplot(2, 2, 1);
imshow(Diff_0, [0, max_diff]);
title('Log-Diff: 0 deg vs Correct');
colorbar;

% 2
subplot(2, 2, 2);
imshow(Diff_1, [0, max_diff]);
title('Log-Diff: \theta/3 deg vs Correct');
colorbar;

% 3
subplot(2, 2, 3);
imshow(Diff_2, [0, max_diff]);
title('Log-Diff: 2\theta/3 deg vs Correct');
colorbar;

% 4
subplot(2, 2, 4);
imshow(Diff_3, [0, max_diff]);
title('Log-Diff: \theta deg vs Correct');
colorbar;

sgtitle('K-Space Log-Magnitude Difference (Black = Match)', 'Interpreter', 'none');

%% --- ローカル関数 ---

function save_raw_data(filepath, data)
    fid = fopen(filepath, 'w');
    if fid == -1, error('ファイルが開けませんでした'); end
    fwrite(fid, data, 'double'); fclose(fid);
end

function k_space_line_signal = simulate_kx_line_fast(image_data, background_phase, kx_row_index, ky_col_indices)
    [Nx, Ny] = size(image_data);
    rho_eff = ifftshift(image_data .* background_phase); 
    x_vec = ((0:Nx-1) / Nx).'; 
    y_vec = ((0:Ny-1) / Ny);   
    k_x_scalar = (kx_row_index) - (floor(Nx/2) + 1);
    k_y_vec = (ky_col_indices(:)) - (floor(Ny/2) + 1);
    
    exp_kx_vec = exp(-1i * 2 * pi * k_x_scalar * x_vec); 
    term_x = exp_kx_vec.' * rho_eff; 
    
    exponent_y = -1i * 2 * pi * (y_vec.' * k_y_vec.'); 
    E_ky = exp(exponent_y);
    k_line_signal_smooth = term_x * E_ky;
    
    phase_shift_x = (-1) .^ k_x_scalar;
    phase_shift_y = (-1) .^ k_y_vec.';
    k_space_line_signal = k_line_signal_smooth .* phase_shift_x .* phase_shift_y;
end

function k_space_3d = fft2_3d_slice_by_slice(image_3d)
    [rows, cols, num_slices] = size(image_3d);
    k_space_3d = complex(zeros(rows, cols, num_slices));
    for slice_idx = 1:num_slices
        current_slice = image_3d(:, :, slice_idx);
        k_space_3d(:, :, slice_idx) = fftshift(fft2(current_slice));
    end
end