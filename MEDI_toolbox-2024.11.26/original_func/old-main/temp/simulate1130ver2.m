%==================================================================================================
% QSM RAWデータ (3D) シミュレーション
% 【改変版: 3段階 階段状ROI置換 (指定座標で角度切り替え)】
%==================================================================================================

fprintf('スクリプトを開始します (3段階回転切り替え版)\n');
clear variables;
close all;

% ★★★ 計算負荷設定 ★★★
DS_FACTOR = 1; 
fprintf('ダウンサンプリング係数: %d\n', DS_FACTOR);

%% --- 1. パラメータ設定 (パス等は環境に合わせてください) ---
image_file_2DGE_0deg_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_2DGE_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_Rotate_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H_local = 'C:\Users\hamaguchi\Downloads\matlab\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
% ダミー定義（エラー回避用）
image_file_1 = '1_data'; image_file_3 = '3_output_data'; image_file_4 = '4_rolate_output_data';
image_file_5 = '5_fitting_output_data'; % スペース修正

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_file_0 = image_file_2DGE_0deg_H;
mag_filename = '1st_2DGE_0deg_mag.raw';
phase_filename = '1st_2DGE_0deg_phase.raw';


% --- ★変更点: 回転切り替え用パラメータ ---
alpha = 3;       % 分割数
beta  = 352;     % 全列数 (cutted_matrix_y に相当)
gamma = 100;     % 開始列 (pix_start_col に相当)

pix_start_row = 116; % 開始行 (User指定: 116行目から)
pix_start_col = gamma; 
width = 224 - pix_start_row; % 置き換える行数


% --- パラメータ ---
params = struct();
params.original_matrix_size = [512, 512, 23];
extention = 2.0/1.3;
theta = -18.6;
rotation_axis = [0 0 1];

% 保存
filename_base = sprintf('Artifact_th%.1f_StepShape_Alpha%d_Beta%d_Gamma%d', theta, alpha, beta, gamma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_base_path = fullfile(image_file_0, image_file_1);
load_mask_path = fullfile(image_file_0, image_file_3);
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

params.matrix_size = size(iMag_4D);
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);

original_img = iMag_4D .* exp(1i * iPhase_4D);
Highpass_img = iMag_4D .* exp(1i * (RDF_small));
Back_img = exp(1i * (iFreq_small - RDF_small));


%% --- 3. 拡張 & 4. 回転 (スライスループ内で実施に変更) ---
fprintf('\n3-4. 実空間拡張とパラメータ準備...\n');
magnification = round(matrix_y * extention);
y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(matrix_y / 2);
y_end_final = y_start_final + matrix_y - 1;

extend_org = complex(zeros(matrix_x, magnification, num_slices));
extend_org(:,y_start_final:y_end_final, :) = original_img;

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

% --- スライスループ ---
for slice_idx = 1 : num_slices
    fprintf('  Processing Slice %d / %d ...\n', slice_idx, num_slices);
    
    % ベース画像 (0度: オリジナル)
    img_slice_0 = extend_org(:,:,slice_idx);
    
    % --- 3つの回転画像を作成 (スライスごとに生成してメモリ節約) ---
    % 回転関数
    rot_func = @(img, th) complex(...
        imrotate(real(img), th, 'bilinear', 'crop'), ...
        imrotate(imag(img), th, 'bilinear', 'crop'));
    
    img_slice_rot1 = rot_func(img_slice_0, angles(1)); % theta / alpha
    img_slice_rot2 = rot_func(img_slice_0, angles(2)); % 2*theta / alpha
    img_slice_rot3 = rot_func(img_slice_0, angles(3)); % 3*theta / alpha (=theta)
    
    % このスライスのk空間バッファ
    temp_k_space_slice = fftshift(fft2(img_slice_0)); % 初期値は0度データで埋める
    
    % --- 行ごとの処理 ---
    for i = 1:length(kx_collect_indices)
        kx_val = kx_collect_indices(i); % 絶対座標 (行)
        
        % 現在の行番号 (ROI内での相対行番号: 1開始)
        current_rel_row = kx_val - x_start_cut + 1;
        
        % 列インデックスの相対値 (1 ～ beta)
        % ky_collect_indices は絶対座標だが、ロジック判定用に1からの連番を作成
        rel_cols = 1:length(ky_collect_indices);
        
        % マスクの初期化
        mask_rot1 = false(size(rel_cols));
        mask_rot2 = false(size(rel_cols));
        mask_rot3 = false(size(rel_cols));
        
        % --- ★★★ 条件分岐ロジック (指定された行・列で切り替え) ★★★ ---
        % User指定: 
        % 1. Row 116, Col 100 -> Row 117, Col 252 : Rot 1
        % 2. Row 117, Col 253 -> Row 118, Col 100 : Rot 2
        % 3. Row 118, Col 101 -> End              : Rot 3
        
        if current_rel_row == 116
            % 116行目: 100列目以降が Rot1 (それ以前は0度=デフォルトのまま)
            mask_rot1 = (rel_cols >= 100);
            
        elseif current_rel_row == 117
            % 117行目: 252列目までは Rot1, 253列目からは Rot2
            mask_rot1 = (rel_cols <= 252);
            mask_rot2 = (rel_cols >= 253);
            
        elseif current_rel_row == 118
            % 118行目: 100列目までは Rot2, 101列目からは Rot3
            mask_rot2 = (rel_cols <= 100);
            mask_rot3 = (rel_cols >= 101);
            
        elseif current_rel_row > 118
            % 118行目以降: すべて Rot3
            mask_rot3 = true(size(rel_cols));
        end
        
        % --- 信号生成と合成 ---
        % 体動シミュレーション用位相項 (簡易版: 今回はスカラー1と仮定)
        current_back = 1; 

        % 合成用ラインバッファ (初期値は0度データから計算、または既存値)
        % ここでは厳密にするため、0度画像から当該ラインを計算してベースにする
        line_base = simulate_kx_line_fast(img_slice_0, current_back, kx_val, ky_collect_indices);
        
        % Rot1 の部分を上書き
        if any(mask_rot1)
            line_rot1 = simulate_kx_line_fast(img_slice_rot1, current_back, kx_val, ky_collect_indices);
            line_base(mask_rot1) = line_rot1(mask_rot1);
        end
        
        % Rot2 の部分を上書き
        if any(mask_rot2)
            line_rot2 = simulate_kx_line_fast(img_slice_rot2, current_back, kx_val, ky_collect_indices);
            line_base(mask_rot2) = line_rot2(mask_rot2);
        end
        
        % Rot3 の部分を上書き
        if any(mask_rot3)
            line_rot3 = simulate_kx_line_fast(img_slice_rot3, current_back, kx_val, ky_collect_indices);
            line_base(mask_rot3) = line_rot3(mask_rot3);
        end
        
        % 計算結果をk空間に格納
        temp_k_space_slice(kx_val, ky_collect_indices) = line_base;
    end
    
    k_space_artifact(:,:,slice_idx) = temp_k_space_slice;
end
elapsedTime = toc;
fprintf('計算完了: %.2f 秒\n', elapsedTime);


%% --- 6. ハイブリッド化と再構成 ---
fprintf('\n6. ハイブリッド化と再構成...\n');

% 元画像(0度)のk空間
k_space_org_3D = fft2_3d_slice_by_slice(extend_org); 

% Step 5で「ROI内の全ライン」について、条件に従って0度/Rot1/Rot2/Rot3を
% 混ぜ合わせた結果(k_space_artifact)を作成済みです。
% したがって、ここでは作成したアーチファクトデータを単純にベース画像に埋め込みます。

final_k_space = k_space_org_3D; % ベースはオリジナル

% シミュレーションしたROI範囲のみを置換
final_k_space(kx_collect_indices, ky_collect_indices, :) = k_space_artifact(kx_collect_indices, ky_collect_indices, :);

% 再構成
artifact_img_3D = complex(zeros(matrix_x, matrix_y, num_slices));

for slice_idx = 1:num_slices
    artifact_img_ext = ifft2(ifftshift(final_k_space(:,:,slice_idx)));
    artifact_img_3D(:,:,slice_idx) = artifact_img_ext(:, y_start_final:y_end_final);
end

% 表示
slice_idx = round(num_slices/2);
figure('WindowState', 'maximized');
subplot(1,2,1); imshow(abs(iMag_4D(:,:,slice_idx)),[]); title('Original');
subplot(1,2,2); imshow(abs(artifact_img_3D(:,:,slice_idx)),[]); title('Multi-Step Motion Artifact');
drawnow;

save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(artifact_img_3D));
fprintf('保存完了: %s\n', save_path);


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