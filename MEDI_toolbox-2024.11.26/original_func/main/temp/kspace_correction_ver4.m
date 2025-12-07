%==================================================================================================
% QSM RAWデータ (3D) シミュレーション
% 【改変版: 3段階 階段状ROI置換 + Figure 2枚分割表示】
%==================================================================================================

fprintf('スクリプトを開始します (Figure 2枚分割版)\n');
clear variables;
close all;

% ★★★ 計算負荷設定 ★★★
DS_FACTOR = 2; 
fprintf('ダウンサンプリング係数: %d\n', DS_FACTOR);

%% --- 1. パラメータ設定 ---
image_file_2DGE_0deg_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!

image_file_1 = '1_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 

correct_mag_filename = '1st_2DGE_1_2_Rotate_mag.raw';
correct_phase_filename = '1st_2DGE_1_2_Rotate_phase.raw';

mag_filename = '1st_2DGE_0deg_mag.raw';
phase_filename = '1st_2DGE_0deg_phase.raw';

% --- パラメータ ---
alpha = 3;       % 分割数
beta  = 352;     % 全列数
gamma = 100;     % 開始列
pix_start_row = 116; 
pix_start_col = gamma; 
target_slice = 20; % 比較対象スライス

params = struct();
params.original_matrix_size = [512, 512, 23];
extention = 2.0/1.3;
theta = -18.6;
rotation_axis = [0 0 1];

filename_base = sprintf('Artifact_th%.1f_StepShape_Alpha%d_Beta%d_Gamma%d', theta, alpha, beta, gamma);

load_base_path = fullfile(image_file_2DGE_0deg_H, image_file_1);
load_mask_path = fullfile(image_file_2DGE_0deg_H, image_file_3);
correct_load_base_path = fullfile(image_file_2DGE_1_2_Rotate_H, image_file_1);
correct_load_mask_path = fullfile(image_file_2DGE_1_2_Rotate_H, image_file_3);
save_path = fullfile(image_file_0, image_file_4);
if ~exist(save_path, 'dir'), mkdir(save_path); end

% --- ハイブリッド化パラメータ ---
cutted_matrix_x = round(224 / DS_FACTOR);
cutted_matrix_y = round(beta / DS_FACTOR);
pix_start_row = round(pix_start_row / DS_FACTOR);
pix_start_col = round(pix_start_col / DS_FACTOR);
gamma = round(gamma / DS_FACTOR);

fprintf('Input Path:   %s\n', fullfile(load_base_path, mag_filename));
fprintf('Correct Path: %s\n', fullfile(correct_load_base_path, correct_mag_filename));

if strcmp(fullfile(load_base_path, mag_filename), fullfile(correct_load_base_path, correct_mag_filename))
    fprintf('【警告】入力ファイルと正解ファイルが全く同じです！\n');
else
    fprintf('ファイルパスは異なります。OKです。\n');
end

%% --- 2. データ読み込み & ダウンサンプリング ---
fprintf('\n2. データ読み込み中...\n');
dims_orig = params.original_matrix_size;
precision = 'double=>double';

% Input Data (0 deg)
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
iPhase_4D = iFreq(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
RDF_small = RDF(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
iFreq_small = iFreq(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
clear iMag_4D_orig iPhase_4D_orig RDF iFreq;

% Correct Data (Artifact)
correct_fid = fopen(fullfile(correct_load_base_path, correct_mag_filename), 'rb');
correct_data_vector_mag = fread(correct_fid, inf, precision); fclose(correct_fid);
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
correct_iPhase_4D = iFreq(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
clear correct_iMag_4D_orig correct_iPhase_4D_orig RDF iFreq; 

params.matrix_size = size(iMag_4D);
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);

correct_original_img = correct_iMag_4D .* exp(1i * correct_iPhase_4D);
original_img = iMag_4D .* exp(1i * iPhase_4D);


%% --- 3. 拡張 & 4. 回転準備 ---
fprintf('\n3-4. 実空間拡張とパラメータ準備...\n');
magnification = round(matrix_y * extention);
y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(matrix_y / 2);
y_end_final = y_start_final + matrix_y - 1;

extend_org = complex(zeros(matrix_x, magnification, num_slices));
extend_org(:,y_start_final:y_end_final, :) = original_img;
correct_extend_org = complex(zeros(matrix_x, magnification, num_slices));
correct_extend_org(:,y_start_final:y_end_final, :) = correct_original_img;

% 回転角度
angles = [theta / alpha, (2 * theta) / alpha, (3 * theta) / alpha];
fprintf('回転角度設定: %.2f度, %.2f度, %.2f度\n', angles(1), angles(2), angles(3));


%% --- 7. 正解データとの比較 (画像 & k空間) ---
fprintf('\n7. 正解データとの比較画像を生成中...\n');

% Log変換用無名関数
get_log_mag = @(k) log(abs(k) + 1);

% 回転関数
rot_func_crop = @(img, ang) complex(...
    imrotate(real(img), ang, 'bilinear', 'crop'), ...
    imrotate(imag(img), ang, 'bilinear', 'crop'));

% ----------------------------------------------------
% A. リファレンス(正解)データの準備
% ----------------------------------------------------
% 実画像 (Magnitude)
Ref_Slice_Img = abs(correct_extend_org(:,:,target_slice));

% k空間 (Log Magnitude)
Ref_Slice_K_Raw = fftshift(fft2(correct_extend_org(:,:,target_slice)));
Ref_Slice_K_Log = get_log_mag(Ref_Slice_K_Raw);


% ----------------------------------------------------
% B. 比較対象(シミュレーション)データの生成と差分計算
% ----------------------------------------------------
Base_Complex_Slice = extend_org(:,:,target_slice); % Input (0deg)

% --- 1. 回転なし (0度) ---
Img_0_ext = rot_func_crop(Base_Complex_Slice, 0);
% 実画像差分
Diff_Img_0 = imabsdiff(Ref_Slice_Img, abs(Img_0_ext));
% k空間差分
Img_0_K_Log = get_log_mag(fftshift(fft2(Img_0_ext)));
Diff_K_0 = abs(Ref_Slice_K_Log - Img_0_K_Log);

% --- 2. theta/3 回転 ---
Img_1_ext = rot_func_crop(Base_Complex_Slice, theta/3);
% 実画像差分
Diff_Img_1 = imabsdiff(Ref_Slice_Img, abs(Img_1_ext));
% k空間差分
Img_1_K_Log = get_log_mag(fftshift(fft2(Img_1_ext)));
Diff_K_1 = abs(Ref_Slice_K_Log - Img_1_K_Log);

% --- 3. 2*theta/3 回転 ---
Img_2_ext = rot_func_crop(Base_Complex_Slice, (2*theta)/3);
% 実画像差分
Diff_Img_2 = imabsdiff(Ref_Slice_Img, abs(Img_2_ext));
% k空間差分
Img_2_K_Log = get_log_mag(fftshift(fft2(Img_2_ext)));
Diff_K_2 = abs(Ref_Slice_K_Log - Img_2_K_Log);

% --- 4. theta 回転 ---
Img_3_ext = rot_func_crop(Base_Complex_Slice, theta);
% 実画像差分
Diff_Img_3 = imabsdiff(Ref_Slice_Img, abs(Img_3_ext));
% k空間差分
Img_3_K_Log = get_log_mag(fftshift(fft2(Img_3_ext)));
Diff_K_3 = abs(Ref_Slice_K_Log - Img_3_K_Log);


% ----------------------------------------------------
% C. Figure表示 (2枚に分割, 各2x2配置)
% ----------------------------------------------------
% スケーリング決定 (それぞれの最大誤差に合わせる)
max_diff_img = max([max(Diff_Img_0(:)), max(Diff_Img_1(:)), max(Diff_Img_2(:)), max(Diff_Img_3(:))]);
if max_diff_img == 0, max_diff_img = 1; end

max_diff_k = max([max(Diff_K_0(:)), max(Diff_K_1(:)), max(Diff_K_2(:)), max(Diff_K_3(:))]);
if max_diff_k == 0, max_diff_k = 1; end


% --- Figure 1: 実空間 (Image) の差分 (2行2列) ---
figure('Name', 'Comparison 1: Real Image Space', 'WindowState', 'maximized');
colormap(gray(256));

subplot(2, 2, 1); imshow(Diff_Img_0, [0, max_diff_img]);
title({'[Image Diff]', '0 deg vs Correct'}); colorbar;

subplot(2, 2, 2); imshow(Diff_Img_1, [0, max_diff_img]);
title({'[Image Diff]', '\theta/3 deg vs Correct'}); colorbar;

subplot(2, 2, 3); imshow(Diff_Img_2, [0, max_diff_img]);
title({'[Image Diff]', '2\theta/3 deg vs Correct'}); colorbar;

subplot(2, 2, 4); imshow(Diff_Img_3, [0, max_diff_img]);
title({'[Image Diff]', '\theta deg vs Correct'}); colorbar;

sgtitle(sprintf('Real Image Difference (Slice: %d)\nBlack=Match, White=Mismatch', target_slice));


% --- Figure 2: k空間 (K-Space) の差分 (2行2列) ---
figure('Name', 'Comparison 2: K-Space Log-Mag', 'WindowState', 'maximized');
colormap(gray(256));

subplot(2, 2, 1); imshow(Diff_K_0, [0, max_diff_k]);
title({'[K-Space Log Diff]', '0 deg vs Correct'}); colorbar;

subplot(2, 2, 2); imshow(Diff_K_1, [0, max_diff_k]);
title({'[K-Space Log Diff]', '\theta/3 deg vs Correct'}); colorbar;

subplot(2, 2, 3); imshow(Diff_K_2, [0, max_diff_k]);
title({'[K-Space Log Diff]', '2\theta/3 deg vs Correct'}); colorbar;

subplot(2, 2, 4); imshow(Diff_K_3, [0, max_diff_k]);
title({'[K-Space Log Diff]', '\theta deg vs Correct'}); colorbar;

sgtitle(sprintf('K-Space Log Difference (Slice: %d)\nBlack=Match, White=Mismatch', target_slice));


%% --- ローカル関数 ---
function save_raw_data(filepath, data)
    fid = fopen(filepath, 'w');
    if fid == -1, error('ファイルが開けませんでした'); end
    fwrite(fid, data, 'double'); fclose(fid);
end