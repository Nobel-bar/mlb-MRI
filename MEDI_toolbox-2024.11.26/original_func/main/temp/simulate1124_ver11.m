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
MOTION_PE_ERROR_MAX_RAD = 0;

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
motion_start_line_shifted = x_start_cut + pix_start_row;

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
    for i = 1:length(kx_collect_indices)
        kx_val = kx_collect_indices(i); % 実際のkx座標
        
        % 体動シミュレーション判定
        current_back_slice = back_slice;
        if kx_val >= motion_start_line_shifted
             % ランダム位相誤差
             phase_error_val = (rand() - 0.5) * MOTION_PE_ERROR_MAX_RAD;
             current_back_slice = back_slice .* exp(1i * phase_error_val);
             current_back_slice = back_slice;
        end

        % ★高速化シミュレーション関数呼び出し
        k_line_signal = simulate_kx_line_fast(img_slice, current_back_slice, kx_val, ky_collect_indices);

        % ★通常のforなので、直接代入が可能（エラーになりません）
        temp_k_space_slice(kx_val, ky_collect_indices) = k_line_signal;
    end
    
    % 計算結果をメイン配列に戻す
    k_space_artifact(:,:,slice_idx) = temp_k_space_slice;
end
elapsedTime = toc;
fprintf('計算完了: %.2f 秒\n', elapsedTime);


%% --- 5.3 以降: ハイブリッド化・再構成・保存 (同じ) ---
fprintf('\n6. ハイブリッド化と再構成...\n');

k_space_org_3D = fft2_3d_slice_by_slice(extend_org); 

hybrid_row_indices = (x_start_cut + pix_start_row - 1) : (x_start_cut + pix_start_row + width - 2); 
hybrid_col_indices =  ky_collect_indices;


base_k_space_artifact = k_space_org_3D;
base_k_space_artifact(hybrid_row_indices, hybrid_col_indices, :) = k_space_artifact(hybrid_row_indices, hybrid_col_indices, :);

base_k_space_direct = k_space_org_3D;
base_k_space_direct(hybrid_row_indices, hybrid_col_indices, :) = k_space_direct_rotated(hybrid_row_indices, hybrid_col_indices, :);

% 再構成

% artifact_img_3D = ifft2_3d_slice_by_slice(base_k_space_artifact);
% direct_img_3D = ifft2_3d_slice_by_slice(base_k_space_direct);

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
% ★★★ 行列演算版シミュレーション関数 (方向修正版) ★★★
% =========================================================
function k_space_line_signal = simulate_kx_line_fast(image_data, background_phase, kx_row_index, ky_col_indices)
    % image_data: 実空間画像 [Nx, Ny]
    % kx_row_index: 計算したい「行」のインデックス (スカラー) -> 物理的には kx (縦周波数)
    % ky_col_indices: 計算したい「列」の範囲 (ベクトル) -> 物理的には ky (横周波数)
    
    [Nx, Ny] = size(image_data);
    
    % --- 1. 画像データ準備 (シフトなし) ---
    rho_eff = image_data .* background_phase; 
    
    % --- 2. 空間座標定義 (0 ～ 1) ---
    x_vec = ((0:Nx-1) / Nx).'; % [Nx x 1] (列ベクトル)
    y_vec = ((0:Ny-1) / Ny);   % [1 x Ny] (行ベクトル)
    
    % --- 3. k空間インデックスの変換 (物理周波数へ) ---
    % 第3引数は「行(kx)」、第4引数は「列(ky)」として扱います
    
    % kx (縦周波数): スカラー
    k_x_scalar = (kx_row_index) - (floor(Nx/2) + 1);
    
    % ky (横周波数): ベクトル
    k_y_vec = (ky_col_indices(:)) - (floor(Ny/2) + 1);
    
    % --- 4. 信号計算 (行列演算) ---
    % 目的: 横方向のライン(Row)を計算したい
    % S(kx, ky) = sum_y [ e^{-i 2pi ky y} * ( sum_x rho(x,y) e^{-i 2pi kx x} ) ]
    
    % Step A: X方向(縦)の積分
    % まず、固定された kx に対して、縦方向の位相回転を与えて足し合わせます。
    
    % X方向の位相項 [Nx x 1]
    exp_kx_vec = exp(-1i * 2 * pi * k_x_scalar * x_vec); 
    
    % 縦方向(X)の内積をとります。
    % [1 x Nx] * [Nx x Ny] -> [1 x Ny]
    % これで「X方向につぶれた(投影された)1次元データ」ができます
    term_x = exp_kx_vec.' * rho_eff; 
    
    % Step B: Y方向(横)の変換
    % 次に、必要な ky (横周波数) の分だけ計算します。
    
    % Y方向の位相項行列 [Ny x M]
    % (y_vec: 1xNy) と (k_y_vec: Mx1) の直積で作る位相マップの転置
    % exponent: [Ny x M]
    exponent_y = -1i * 2 * pi * (y_vec.' * k_y_vec.'); 
    E_ky = exp(exponent_y);
    
    % [1 x Ny] * [Ny x M] -> [1 x M]
    % 最終的に、指定された列(ky)の分だけの信号が横ベクトルとして出ます
    k_space_line_signal = term_x * E_ky;
end


function k_space_3d = fft2_3d_slice_by_slice(image_3d)
% FFT2_3D_SLICE_BY_SLICE 3D画像に対してスライスごとの2D-FFTを実行する関数
%
%   [Input]
%       image_3d : (Nx, Ny, Nz) の3次元画像データ (実数または複素数)
%
%   [Output]
%       k_space_3d : (Nx, Ny, Nz) のk空間データ (複素数)
%                    各スライスに対して fftshift(fft2(...)) が適用されています。

    % 1. 入力サイズの取得
    [rows, cols, num_slices] = size(image_3d);

    % 2. 出力用配列の事前確保 (メモリ効率化のため)
    % 結果は必ず複素数になるため、complexで初期化します
    k_space_3d = complex(zeros(rows, cols, num_slices));

    % 3. スライスごとのループ処理
    for slice_idx = 1:num_slices
        % 各スライスを取り出し、2次元フーリエ変換 -> 中心シフトr
        current_slice = image_3d(:, :, slice_idx);
        k_space_3d(:, :, slice_idx) = fftshift(fft2(current_slice));
    end
end


function image_3d = ifft2_3d_slice_by_slice(k_space_3d)
% IFFT2_3D_SLICE_BY_SLICE 3D k空間データに対してスライスごとの2D逆FFTを実行する関数
%
%   [Input]
%       k_space_3d : (Nx, Ny, Nz) のk空間データ (通常は複素数)
%
%   [Output]
%       image_3d   : (Nx, Ny, Nz) の実空間画像 (複素数)
%                    各スライスに対して ifft2(ifftshift(...)) が適用されています。

    % 1. 入力サイズの取得
    [rows, cols, num_slices] = size(k_space_3d);

    % 2. 出力用配列の事前確保
    % 逆変換結果も複素数になるため complex で初期化
    image_3d = complex(zeros(rows, cols, num_slices));

    % 3. スライスごとのループ処理
    for slice_idx = 1:num_slices
        % 各スライスを取り出す
        current_k_slice = k_space_3d(:, :, slice_idx);
        
        % ifftshift (DC成分を中心から四隅へ戻す) -> ifft2 (逆フーリエ変換)
        image_3d(:, :, slice_idx) = ifft2(ifftshift(current_k_slice));
    end
end