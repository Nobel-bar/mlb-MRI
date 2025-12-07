%==================================================================================================
% QSM RAWデータ (3D) シミュレーション
% 【修正版: 0度画像を入力とし、18度アーチファクト画像を生成】
%==================================================================================================

fprintf('スクリプトを開始します (Input修正版)\n');
clear variables;
close all;

% ★★★ 計算負荷設定 ★★★
DS_FACTOR = 2; % 必要に応じて 2, 4 に変更
fprintf('ダウンサンプリング係数: %d\n', DS_FACTOR);

%% --- 1. パラメータ設定 ---
% フォルダパス定義
image_file_2DGE_0deg_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; 
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; 
% image_file_2DGE_Rotate_H = '...'; % 今回はInputとしては使いません

% 作業用フォルダ定義
image_file_1 = '1_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 

% ベースパス設定 (保存先など)
image_file_0 = image_file_2DGE_0deg_H;

% --- 【重要】ファイル設定の修正 ---

% 1. 正解データ (Target): アーチファクトありの18度画像
correct_load_base_path = fullfile(image_file_2DGE_1_2_Rotate_H, image_file_1);
correct_mag_filename = '1st_2DGE_1_2_Rotate_mag.raw';
correct_phase_filename = '1st_2DGE_1_2_Rotate_phase.raw';

% 2. 入力データ (Input): ★ここを「0度 (回転なし)」に戻します★
% これをベースに回転シミュレーションを行うことで、18度の画像を作ります
load_base_path = fullfile(image_file_2DGE_0deg_H, image_file_1);
mag_filename = '1st_2DGE_0deg_mag.raw';       
phase_filename = '1st_2DGE_0deg_phase.raw';

% ※ 前回間違いの原因だった「Rotateファイルを入力にする」記述は削除しました

% マスク等のパス
load_mask_path = fullfile(image_file_2DGE_0deg_H, image_file_3);
correct_load_mask_path = fullfile(image_file_2DGE_1_2_Rotate_H, image_file_3);

% 保存先
save_path = fullfile(image_file_0, image_file_4);
if ~exist(save_path, 'dir'), mkdir(save_path); end


% --- 回転切り替え用パラメータ ---
alpha = 3;       % 分割数
beta  = 352;     % 全列数
gamma = 100;     % 開始列
pix_start_row = 116; % 開始行 (116行目以降を回転させる -> 画像の大部分が回転する)
pix_start_col = gamma; 

% --- 共通パラメータ ---
params = struct();
params.original_matrix_size = [512, 512, 23];
extention = 2.0/1.3;
theta = -18.6; % 0度からこの角度へ回転させる
rotation_axis = [0 0 1];

filename_base = sprintf('Artifact_th%.1f_StepShape_Alpha%d_Beta%d_Gamma%d', theta, alpha, beta, gamma);

% --- ハイブリッド化パラメータ ---
cutted_matrix_x = round(224 / DS_FACTOR);
cutted_matrix_y = round(beta / DS_FACTOR);

% ダウンサンプリング対応
pix_start_row = round(pix_start_row / DS_FACTOR);
pix_start_col = round(pix_start_col / DS_FACTOR);
gamma = round(gamma / DS_FACTOR);


%% --- 2. データ読み込み & ダウンサンプリング ---
fprintf('\n2. データ読み込み中...\n');

dims_orig = params.original_matrix_size;
precision = 'double=>double';

% 2.1 入力データ (0度 Clean)
fid = fopen(fullfile(load_base_path, mag_filename), 'rb');
if fid == -1, error('Inputファイルが見つかりません: %s', fullfile(load_base_path, mag_filename)); end
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

% 2.2 正解データ (18度 Artifact)
fid = fopen(fullfile(correct_load_base_path, correct_mag_filename), 'rb');
if fid == -1, error('Correctファイルが見つかりません: %s', fullfile(correct_load_base_path, correct_mag_filename)); end
correct_data_vector_mag = fread(fid, inf, precision); fclose(fid);
correct_iMag_4D_orig = reshape(correct_data_vector_mag, dims_orig);

fid = fopen(fullfile(correct_load_base_path, correct_phase_filename), 'rb');
correct_iPhase_4D_orig = reshape(fread(fid, inf, precision), dims_orig); fclose(fid);

correct_iMag_4D = correct_iMag_4D_orig(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
correct_iPhase_4D = correct_iPhase_4D_orig(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
clear correct_iMag_4D_orig correct_iPhase_4D_orig;

params.matrix_size = size(iMag_4D);
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);

correct_original_img = correct_iMag_4D .* exp(1i * correct_iPhase_4D);
original_img = iMag_4D .* exp(1i * iPhase_4D); % これは0度画像


%% --- 3. 拡張 & 4. 回転準備 ---
fprintf('\n3-4. 実空間拡張とパラメータ準備...\n');
magnification = round(matrix_y * extention);
y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(matrix_y / 2);
y_end_final = y_start_final + matrix_y - 1;

extend_org = complex(zeros(matrix_x, magnification, num_slices));
extend_org(:,y_start_final:y_end_final, :) = original_img; % 0度ベース

correct_extend_org = complex(zeros(matrix_x, magnification, num_slices));
correct_extend_org(:,y_start_final:y_end_final, :) = correct_original_img; % 18度正解

% 回転角度定義 (Input=0度に対して、以下の角度を適用)
angles = [theta / alpha, (2 * theta) / alpha, (3 * theta) / alpha];
fprintf('回転角度設定: %.2f度, %.2f度, %.2f度\n', angles(1), angles(2), angles(3));


%% --- 5. k空間シミュレーション ---
fprintf('\n5. k空間シミュレーション (0deg -> 18deg Transition)...\n');
tic; 

x_center = floor(matrix_x / 2) + 1;
x_start_cut = x_center - floor(cutted_matrix_x / 2);
x_end_cut = x_start_cut + cutted_matrix_x - 1;
y_center_org = floor(magnification / 2) + 1;
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2);
y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1;

kx_collect_indices = x_start_cut : x_end_cut;
ky_collect_indices = y_start_org_cut : y_end_org_cut;

k_space_artifact = complex(zeros(matrix_x, magnification, num_slices));

% --- スライスループ ---
for slice_idx = 1 : num_slices
    fprintf('  Processing Slice %d / %d ...\n', slice_idx, num_slices);
    
    img_slice_0 = extend_org(:,:,slice_idx); % 0度画像
    
    % 回転関数
    rot_func = @(img, th) complex(...
        imrotate(real(img), th, 'bilinear', 'crop'), ...
        imrotate(imag(img), th, 'bilinear', 'crop'));
    
    % 各角度の回転画像を作成 (ここでの th は 0度からの変位)
    img_slice_rot1 = rot_func(img_slice_0, angles(1));
    img_slice_rot2 = rot_func(img_slice_0, angles(2));
    img_slice_rot3 = rot_func(img_slice_0, angles(3)); % これが theta (18度)
    
    % k空間バッファ (ベースは0度画像から生成)
    % ※ 116行目以前(maskにかからない部分)は、この0度データがそのまま残ります
    temp_k_space_slice = fftshift(fft2(img_slice_0)); 
    
    % --- 行ごとの処理 ---
    for i = 1:length(kx_collect_indices)
        kx_val = kx_collect_indices(i); 
        current_rel_row = kx_val - x_start_cut + 1; 
        rel_cols = 1:length(ky_collect_indices);
        
        mask_rot1 = false(size(rel_cols));
        mask_rot2 = false(size(rel_cols));
        mask_rot3 = false(size(rel_cols));
        
        % 条件分岐
        if current_rel_row == 116
            mask_rot1 = (rel_cols >= 100);
        elseif current_rel_row == 117
            mask_rot1 = (rel_cols <= 252);
            mask_rot2 = (rel_cols >= 253);
        elseif current_rel_row == 118
            mask_rot2 = (rel_cols <= 100);
            mask_rot3 = (rel_cols >= 101);
        elseif current_rel_row > 118
            % 119行目以降 (画像の中心を含む大部分) は Rot3 (18度) になる
            % つまり、画像のメインの見た目は18度になる
            mask_rot3 = true(size(rel_cols));
        end
        
        current_back = 1;
        
        % ベースライン (0度)
        line_base = simulate_kx_line_fast(img_slice_0, current_back, kx_val, ky_collect_indices);
        
        % 回転データの埋め込み
        if any(mask_rot1)
            line_rot1 = simulate_kx_line_fast(img_slice_rot1, current_back, kx_val, ky_collect_indices);
            line_base(mask_rot1) = line_rot1(mask_rot1);
        end
        if any(mask_rot2)
            line_rot2 = simulate_kx_line_fast(img_slice_rot2, current_back, kx_val, ky_collect_indices);
            line_base(mask_rot2) = line_rot2(mask_rot2);
        end
        if any(mask_rot3)
            line_rot3 = simulate_kx_line_fast(img_slice_rot3, current_back, kx_val, ky_collect_indices);
            line_base(mask_rot3) = line_rot3(mask_rot3);
        end
        
        temp_k_space_slice(kx_val, ky_collect_indices) = line_base;
    end
    
    k_space_artifact(:,:,slice_idx) = temp_k_space_slice;
end
elapsedTime = toc;
fprintf('計算完了: %.2f 秒\n', elapsedTime);


%% --- 6. ハイブリッド化と再構成 ---
fprintf('\n6. ハイブリッド化と再構成...\n');

k_space_org_3D = fft2_3d_slice_by_slice(extend_org); % 0度ベース
final_k_space = k_space_org_3D;
final_k_space(kx_collect_indices, ky_collect_indices, :) = k_space_artifact(kx_collect_indices, ky_collect_indices, :);

% 再構成
artifact_img_3D = complex(zeros(matrix_x, matrix_y, num_slices));
for slice_idx = 1:num_slices
    artifact_img_ext = ifft2(ifftshift(final_k_space(:,:,slice_idx)));
    artifact_img_3D(:,:,slice_idx) = artifact_img_ext(:, y_start_final:y_end_final);
end
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(artifact_img_3D));


%% --- 7. 正解データとの比較 (Log-Diff) ---
fprintf('\n7. 正解データとの比較画像を生成中...\n');

get_log_mag = @(k) log(abs(k) + 1);
target_slice = 12;

% 正解データ (18度 Artifact) のk空間
Ref_Slice_K_Raw = fftshift(fft2(correct_extend_org(:,:,target_slice)));
Ref_Slice_K_Log = get_log_mag(Ref_Slice_K_Raw);

% シミュレーション結果 (上記のStep5で作成したもの)
Sim_Result_K_Raw = final_k_space(:,:,target_slice);
Sim_Result_K_Log = get_log_mag(Sim_Result_K_Raw);


% --- 比較用の各段階画像 (今回は Input=0度 なので、回転させて18度に近づくかを見る) ---
Base_Complex_Slice = extend_org(:,:,target_slice); % 0度
rot_func_crop = @(img, ang) complex(...
    imrotate(real(img), ang, 'bilinear', 'crop'), ...
    imrotate(imag(img), ang, 'bilinear', 'crop'));

% 1. Input (0度) との差
Img_0_K_Log = get_log_mag(fftshift(fft2(rot_func_crop(Base_Complex_Slice, 0))));
Diff_Input = abs(Ref_Slice_K_Log - Img_0_K_Log);

% 2. 完全に18度回転させた場合 (理想の18度) との差
Img_Theta_K_Log = get_log_mag(fftshift(fft2(rot_func_crop(Base_Complex_Slice, theta))));
Diff_FullRot = abs(Ref_Slice_K_Log - Img_Theta_K_Log);

% 3. 今回のシミュレーション結果 (Step) との差
Diff_Sim = abs(Ref_Slice_K_Log - Sim_Result_K_Log);

max_diff = max([max(Diff_Input(:)), max(Diff_FullRot(:)), max(Diff_Sim(:))]);
if max_diff == 0, max_diff = 1; end 

figure('Name', 'Simulation Accuracy Check', 'WindowState', 'maximized');
colormap(gray(256));

subplot(1, 3, 1);
imshow(Diff_Input, [0, max_diff]);
title('Diff: Input(0deg) vs Correct(18deg)');
xlabel('Should be WHITE (Mismatch)');

subplot(1, 3, 2);
imshow(Diff_FullRot, [0, max_diff]);
title('Diff: Pure 18deg vs Correct');
xlabel('Should be mostly BLACK');

subplot(1, 3, 3);
imshow(Diff_Sim, [0, max_diff]);
title('Diff: Simulation(Step) vs Correct');
xlabel('Target Result (Should be Blackest)');

sgtitle('Comparison: Can we reproduce the Artifact?', 'Interpreter', 'none');


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
function image_3d = ifft2_3d_slice_by_slice(k_space_3d)
    [rows, cols, num_slices] = size(k_space_3d);
    image_3d = complex(zeros(rows, cols, num_slices));
    for slice_idx = 1:num_slices
        image_3d(:, :, slice_idx) = ifft2(ifftshift(k_space_3d(:, :, slice_idx)));
    end
end