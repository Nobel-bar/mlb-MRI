%==================================================================================================
% QSM RAWデータ (3D) シミュレーション
% 【高速化版: 行列演算 】
%==================================================================================================

fprintf('スクリプトを開始します (高速化・並列化版)\n');
clear variables;
close all;

% --- 並列プールの起動 ---
% % 既にプールが開いている場合は何もしない、なければ開く
% if isempty(gcp('nocreate'))
%     try
%         parpool; % 自動的に利用可能なコア数で起動します
%     catch
%         fprintf('並列化ツールの起動に失敗しました。通常の計算を行います。\n');
%     end
% end

% ★★★ 計算負荷設定 ★★★
DS_FACTOR = 1; % 2=解像度1/2, 4=解像度1/4
fprintf('ダウンサンプリング係数: %d\n', DS_FACTOR);

%% --- 1. パラメータ設定 (パス等は環境に合わせてください) ---
image_file_2DGE_0deg_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_2DGE_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_Rotate_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H_local = 'C:\Users\hamaguchi\Downloads\matlab\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; % スペース修正
image_file_4 = '4_rolate_output_data'; % スペース修正
image_file_5 = '5_fitting_output_data'; % スペース修正

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_file_0 = image_file_2DGE_0deg_H;
mag_filename = '1st_2DGE_0deg_mag.raw';
phase_filename = '1st_2DGE_0deg_phase.raw';

width = 112; % 置き換える行数
pix_start_row = 112; % k空間ROIの何行目から置き換えるか (1-based)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_base_path = fullfile(image_file_0, image_file_1);
load_mask_path = fullfile(image_file_0, image_file_3);
save_path = fullfile(image_file_0, image_file_4);

if ~exist(save_path, 'dir'), mkdir(save_path); end

% --- パラメータ ---
params = struct();
params.original_matrix_size = [512, 512, 23];
extention = 2.0/1.3;
theta = -18.6;
rotation_axis = [0 0 1];

% --- ハイブリッド化パラメータ (スケーリング) ---
cutted_matrix_x = round(224 / DS_FACTOR);
cutted_matrix_y = round(352 / DS_FACTOR);

% --- 体動パラメータ ---
MOTION_START_LINE_ORIG = 257;
MOTION_START_LINE = round(MOTION_START_LINE_ORIG / DS_FACTOR);
MOTION_PE_ERROR_MAX_RAD = pi/4;


%% --- 2. データ読み込み & ダウンサンプリング ---
fprintf('\n2. データ読み込み中...\n');
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


% ダウンサンプリング
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


%% --- 3. 拡張 & 4. 回転 ---
fprintf('\n3-4. 実空間拡張と回転...\n');
magnification = round(matrix_y * extention);
y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(matrix_y / 2);
y_end_final = y_start_final + matrix_y - 1;

extend_org = complex(zeros(matrix_x, magnification, num_slices));
extend_high = complex(zeros(matrix_x, magnification, num_slices));
extend_back = complex(zeros(matrix_x, magnification, num_slices));

extend_org(:,y_start_final:y_end_final, :) = original_img;
extend_high(:,y_start_final:y_end_final, :) = Highpass_img;
extend_back(:,y_start_final:y_end_final, :) = Back_img;

rotate_complex = @(img) complex(...
    imrotate3(real(img), theta, rotation_axis, 'linear', 'crop'), ...
    imrotate3(imag(img), theta, rotation_axis, 'linear', 'crop'));

extend_org_rotated  = rotate_complex(extend_org);
extend_high_rotated = rotate_complex(extend_high);
extend_back_rotated = rotate_complex(extend_back);


%% ---%% --- 5. k空間シミュレーション (parfor未使用・高速化版) ---
fprintf('\n5. k空間シミュレーション (通常forループ・行列演算高速化)...\n');
tic; % 計測開始

x_center = floor(matrix_x / 2) + 1;
x_start_cut = x_center - floor(cutted_matrix_x / 2);
x_end_cut = x_start_cut + cutted_matrix_x - 1;
y_center_org = floor(magnification / 2) + 1;
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2);
y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1;

kx_collect_indices = x_start_cut : x_end_cut;
ky_collect_indices = y_start_org_cut : y_end_org_cut;
motion_start_line_shifted = y_center_org - floor(cutted_matrix_y / 2) + MOTION_START_LINE;

% 結果格納用配列
k_space_artifact = complex(zeros(matrix_x, magnification, num_slices));
k_space_direct_rotated = complex(zeros(matrix_x, magnification, num_slices));

% --- スライスループ ---
for slice_idx = 1 : num_slices
    fprintf('  Processing Slice %d / %d ...\n', slice_idx, num_slices);
    
    % 事前計算
    k_space_direct_rotated(:,:,slice_idx) = fftshift(fft2(extend_org_rotated(:,:,slice_idx)));
    
    % スライスデータの準備
    % 普通のforループなので一時変数へのコピーなどは不要ですが、可読性のためそのまま使います
    temp_k_space_slice = k_space_direct_rotated(:,:,slice_idx);
    
    img_slice = extend_high_rotated(:,:,slice_idx);
    back_slice = extend_back_rotated(:,:,slice_idx);
    
    % --- 通常の for ループ (parfor ではない) ---
    % 行列演算関数を使っているため、これでも十分に高速です
    for i = 1:length(ky_collect_indices)
        ky_val = ky_collect_indices(i); % 実際のky座標
        
        % 体動シミュレーション判定
        current_back_slice = back_slice;
        if ky_val >= motion_start_line_shifted
             % ランダム位相誤差
             phase_error_val = (rand() - 0.5) * MOTION_PE_ERROR_MAX_RAD;
             current_back_slice = back_slice .* exp(1i * phase_error_val);
        end

        % ★高速化シミュレーション関数呼び出し
        k_line_signal = simulate_ky_line_fast(img_slice, current_back_slice, ky_val, kx_collect_indices);

        % ★通常のforなので、直接代入が可能（エラーになりません）
        temp_k_space_slice(kx_collect_indices, ky_val) = k_line_signal;
    end
    
    % 計算結果をメイン配列に戻す
    k_space_artifact(:,:,slice_idx) = temp_k_space_slice;
end
elapsedTime = toc;
fprintf('計算完了: %.2f 秒\n', elapsedTime);


%% --- 5.3 以降: ハイブリッド化・再構成・保存 (同じ) ---
fprintf('\n6. ハイブリッド化と再構成...\n');

k_space_org_3D = k_space_direct_rotated; 

hybrid_row_indices = (x_start_cut + pix_start_row - 1) : (x_start_cut + pix_start_row + width - 2); 
hybrid_col_indices =  y_start_org_cut: y_end_org_cut;


base_k_space_artifact = k_space_org_3D;
base_k_space_artifact(hybrid_row_indices, hybrid_col_indices, :) = k_space_artifact(hybrid_row_indices, hybrid_col_indices, :);

base_k_space_direct = k_space_org_3D;
base_k_space_direct(hybrid_row_indices, hybrid_col_indices, :) = k_space_direct_rotated(hybrid_row_indices, hybrid_col_indices, :);

% 再構成
artifact_img_3D = complex(zeros(matrix_x, matrix_y, num_slices));
direct_img_3D = complex(zeros(matrix_x, matrix_y, num_slices));

for slice_idx = 1:num_slices
    artifact_img_ext = ifft2(ifftshift(base_k_space_artifact(:,:,slice_idx)));
    direct_img_ext = ifft2(ifftshift(base_k_space_direct(:,:,slice_idx)));
    
    artifact_img_3D(:,:,slice_idx) = artifact_img_ext(:, y_start_final:y_end_final);
    direct_img_3D(:,:,slice_idx) = direct_img_ext(:, y_start_final:y_end_final);
end

% 表示
slice_idx = round(num_slices/2);
figure('WindowState', 'maximized');
subplot(1,3,1); imshow(abs(iMag_4D(:,:,slice_idx)),[]); title('Original');
subplot(1,3,2); imshow(abs(artifact_img_3D(:,:,slice_idx)),[]); title('Motion Artifact (Fast)');
subplot(1,3,3); imshow(abs(direct_img_3D(:,:,slice_idx)),[]); title('Direct Rotation');
drawnow;

% 保存
filename_base = sprintf('Fast_Artifact_th%.1f', theta);
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(artifact_img_3D));
fprintf('保存完了: %s\n', save_path);


%% --- ローカル関数 (高速化の肝) ---

function save_raw_data(filepath, data)
    fid = fopen(filepath, 'w');
    if fid == -1, error('ファイルが開けませんでした'); end
    fwrite(fid, data, 'double'); fclose(fid);
end


% =========================================================
% ★★★ 行列演算版シミュレーション関数  ★★★
% =========================================================
function k_space_line_signal = simulate_ky_line_fast(image_data, background_phase, ky_line_index, kx_collect_indices)
    % image_data 実空間画像
    [Nx, Ny] = size(image_data);
        

    
    rho_eff = image_data .* background_phase; 
    
    % --- 2. 空間座標定義 (0 ～ 1) ---
    % 画像の左上を原点(0)とした座標系を作ります。
    x_vec = ((0:Nx-1) / Nx).'; % [Nx x 1] (列ベクトル)
    y_vec = ((0:Ny-1) / Ny);   % [1 x Ny] (行ベクトル)
    
    % --- 3. k空間インデックスの変換 (物理周波数へ) ---
    % 入力されているインデックスは、fftshiftされた状態(1～N)のものです。
    % これを物理的な周波数 k (DC成分が0になる値) に変換します。
    % 例: 512サイズの場合、中心の257番目 → 周波数 0
    k_y_val = (ky_line_index) - (floor(Ny/2) + 1);
    k_x_vals = (kx_collect_indices(:)) - (floor(Nx/2) + 1);
    k_x_vals = kx_collect_indices;
    
    % --- 4. 信号計算 (行列演算) ---
    % S(kx, ky) = sum_x sum_y rho(x,y) * exp(-i 2pi (kx*x + ky*y))
    
    % Y方向の積分項
    exp_ky_vec = exp(-1i * 2 * pi * k_y_val * y_vec).'; 
    term_y = rho_eff * exp_ky_vec; 
    
    % X方向の積分項
    exponent = -1i * 2 * pi * (k_x_vals * x_vec.'); 
    E_kx = exp(exponent);
    
    % 最終計算
    k_space_line_signal = E_kx * term_y;
end