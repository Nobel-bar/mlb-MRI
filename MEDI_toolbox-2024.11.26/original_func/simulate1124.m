%==================================================================================================
% QSM RAWデータ (3D) を読み込み、
% 3D実空間で回転 (imrotate3) をシミュレートし、
% k空間でハイブリッド化 (一部置換) を行うMATLABスクリプト
%
% 【修正版: 体動シミュレーション統合】
% - k_x = 257行目以降でランダムな位相誤差を適用。
% - k空間の行（k_x）ごとのループを明確化。
%==================================================================================================

fprintf('スクリプトを開始します (3D実空間回転 + k空間ハイブリッド + 体動シミュレーション)\n');
clear variables;
close all;

%% --- 1. 撮像・シミュレーション パラメータ設定 ---
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load_base_path = fullfile(image_file_0, image_file_1);
load_mask_path = fullfile(image_file_0, image_file_3);
save_path = fullfile(image_file_0, image_file_4);

if ~exist(save_path, 'dir')
	mkdir(save_path);
	fprintf('保存フォルダを作成しました: %s\n', save_path);
end

% ★★★ 計算負荷軽減のための設定 ★★★
% 各次元を何分の1にするか指定 (2を指定すると面積比で1/4、計算量は約1/16になります)
DS_FACTOR = 2; 
fprintf('計算負荷を下げるため、各次元を 1/%d にダウンサンプリングします。\n', DS_FACTOR);

% --- params構造体 ---
params = struct();
params.matrix_size = [512, 512, 23];
params.TE = [0.015]; 

% --- 拡張と回転のパラメータ ---
extention = 2.0/1.3;
magnification = round( params.matrix_size(2) * extention);
theta = -18.6;
rotation_axis = [0 0 1];

% --- k空間ハイブリッド化パラメータ ---
cutted_matrix_x = 224;
cutted_matrix_y = 352;
cutted_matrix_x = round(224 / DS_FACTOR);
cutted_matrix_y = round(352 / DS_FACTOR);
width = 112; 
pix_start_row = 112; 

% --- ★体動シミュレーション用パラメータ★ ---
% MOTION_START_LINE = 257; % 512x512行列における、体動が発生し始めるk_x行 (1-based)
% % --- ★体動シミュレーション用パラメータ (スケーリング対応)★ ---
MOTION_START_LINE_ORIG = 257;
MOTION_START_LINE = round(MOTION_START_LINE_ORIG / DS_FACTOR);
% MOTION_PE_ERROR_MAX_RAD = pi/4;


%% --- 2. RAWデータの読み込みとダウンサンプリング ---
fprintf('\n2. データを読み込み、ダウンサンプリングしています...\n');

mag_filepath = fullfile(load_base_path, mag_filename);
phase_filepath = fullfile(load_base_path, phase_filename);
dims_orig = params.original_matrix_size;
precision = 'double=>double';

% --- 強度・位相データの読み込み ---
fid_mag = fopen(mag_filepath, 'rb');
if fid_mag == -1, error('強度ファイルを開けませんでした: %s', mag_filepath); end
data_vector_mag = fread(fid_mag, inf, precision);
fclose(fid_mag);

params.original_matrix_size(3) = numel(data_vector_mag) / (dims_orig(1) * dims_orig(2));
dims_orig = params.original_matrix_size;

iMag_4D_orig = reshape(data_vector_mag, dims_orig);
clear data_vector_mag; 

fid_phase = fopen(phase_filepath, 'rb');
if fid_phase == -1, error('位相ファイルを開けませんでした: %s', phase_filepath); end
iPhase_4D_orig = reshape(fread(fid_phase, inf, precision), dims_orig);
fclose(fid_phase);

% 処理済みデータ (iFreq, RDF) の読み込み
try
    load(fullfile(load_mask_path, 'phase.mat'), 'iFreq');
    load(fullfile(load_mask_path, 'PDF.mat'), 'RDF');
    % 変数名が違う場合の対応
    if exist('iFreq_orig', 'var'), iFreq = iFreq_orig; end
    if exist('RDF_orig', 'var'), RDF = RDF_orig; end
catch
    fprintf('phase.mat / PDF.mat が見つかりません。ダミーを使用します。\n');
    iFreq = zeros(dims_orig); RDF = zeros(dims_orig);
end

% === ★★★ ダウンサンプリング実行 ★★★ ===
% 実空間で単純間引き (Nearest Neighbor 相当)
iMag_4D = iMag_4D_orig(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
iPhase_4D = iPhase_4D_orig(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
RDF_small = RDF(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);
iFreq_small = iFreq(1:DS_FACTOR:end, 1:DS_FACTOR:end, :);

% メモリ解放
clear iMag_4D_orig iPhase_4D_orig RDF iFreq;

% サイズ情報を更新
params.matrix_size = size(iMag_4D);
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);

fprintf('  元のサイズ: %d x %d x %d\n', dims_orig(1), dims_orig(2), dims_orig(3));
fprintf('  現在のサイズ: %d x %d x %d (DS_FACTOR=%d)\n', matrix_x, matrix_y, num_slices, DS_FACTOR);


original_img = iMag_4D(:,:,:) .* exp(1i * iPhase_4D(:,:,:));
Highpass_img = iMag_4D(:,:,:) .* exp(1i * (RDF));
Back_img = exp(1i * (iFreq - RDF)); % 背景磁場成分 (位相項)

fprintf('データの読み込み完了。%d スライスを処理します。\n', num_slices);


%% --- 3. 実空間データの拡張 (ゼロパディング) ---
fprintf('\n3. 実空間データのY方向を拡張 (ゼロパディング) しています...\n');
extend_org = complex(zeros(matrix_x, magnification, num_slices));
extend_high = complex(zeros(matrix_x, magnification, num_slices));
extend_back = complex(zeros(matrix_x, magnification, num_slices));

y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(matrix_y / 2);
y_end_final = y_start_final + matrix_y - 1;

extend_org(:,y_start_final:y_end_final, :) = original_img;
extend_high(:,y_start_final:y_end_final, :) = Highpass_img;
extend_back(:,y_start_final:y_end_final, :) = Back_img;

fprintf('    3D実空間の準備が完了しました。\n');


%% --- 4. 3D 実空間の回転（モーションのシミュレート） ---
fprintf('\n4. 3D 実空間の回転 (imrotate3) を実行中...\n');

% Note: imrotate3 が内部関数として定義されていないため、ここでは直接記述
rotate_complex = @(img) complex(...
    imrotate3(real(img), theta, rotation_axis, 'linear', 'crop'), ...
    imrotate3(imag(img), theta, rotation_axis, 'linear', 'crop'));

extend_org_rotated  = rotate_complex(extend_org);
extend_high_rotated = rotate_complex(extend_high);
extend_back_rotated = rotate_complex(extend_back); % 回転後の背景磁場位相マップ (Back_img)
fprintf('    3D 実空間の回転が完了しました。\n'); 

% --- 5.0. 準備 ---
Nx_total = matrix_y; 
Nx_ext = matrix_x;
k_space_artifact = complex(zeros(matrix_x, magnification, num_slices));
k_space_direct_rotated = complex(zeros(matrix_x, magnification, num_slices));

% --- ★体動シミュレーション変数★ ---
is_motion_fixed = false; % モーション発生フラグ
fixed_phase_error_factor = 1; % 発生した位相誤差の exp(i*phi) を保持

 x_center = floor(matrix_x / 2) + 1; % 257
x_start_cut = x_center - floor(cutted_matrix_x / 2); % 145
x_end_cut = x_start_cut + cutted_matrix_x - 1; % 368
y_center_org = floor(magnification / 2) + 1; % 395
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2); % 219
y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1; % 570
kx_collect_indices = x_start_cut : x_end_cut;
ky_collect_indices = y_start_org_cut : y_end_org_cut;

% 体動開始ラインも新しい座標系に合わせる
motion_start_line_shifted = y_center_org - floor(matrix_y/2) + MOTION_START_LINE;


% --- k空間の行（k_x）のインデックスは 1 to 512 ---
fprintf('    k空間シミュレーション計算中... (計算量 1/%d)\n', DS_FACTOR^2*DS_FACTOR);
for slice  = 1 : num_slices
    k_space_direct_rotated(:,:,slice_idx) = fftshift(fft2(extend_org_rotated(:,:,slice_idx)));
    k_space_artifact(:,:,slice_idx) = k_space_direct_rotated(:,:,slice_idx);
    
    % 画像データのスライス抽出
    img_slice_rotated = extend_high_rotated(:,:,slice_idx);
    back_slice_rotated = extend_back_rotated(:,:,slice_idx);

    % --- k空間の行（k_y: 位相エンコード）ごとのループ ---
    for k_y_idx = ky_collect_indices
        k_space_direct_rotated(:,:,num_slices)  = fftshift(fft2(extend_org_rotated(:,:,num_slices) ));
        k_space_artifact(:,:,num_slices) = k_space_direct_rotated(:,:,num_slices);
        [k_line_signal_partial, kx_out_indices] = simulate_ky_line_collection_by_sum(...
        extend_high_rotated, extend_back_rotated, k_x_idx, ky_collect_indices);

        % k_space_artifact の k_y_idx 列の、k_x_out_indices の行に代入
        k_space_artifact(kx_out_indices, k_y_idx, slice_idx) = k_line_signal_partial;
    end
end 
% --- k_x_idx ループの終了 ---


% --- 5.3. ハイブリッド化 (データの置換) ---
fprintf('    k空間データをハイブリッド化しています...\n');

hybrid_row_indices = (x_center - floor(cutted_matrix_x / 2)) : (x_center + floor(cutted_matrix_x / 2) - 1);
hybrid_col_indices = (y_center_final_ext - floor(cutted_matrix_y / 2)) : (y_center_final_ext + floor(cutted_matrix_y / 2) - 1);

% 1. アーチファクト結果をベースk空間にハイブリッド
base_k_space_artifact = k_space_org_3D; % 回転前のk空間をベースにする
base_k_space_artifact(hybrid_row_indices, hybrid_col_indices, :) = k_space_artifact(hybrid_row_indices, hybrid_col_indices, :);

% 2. Direct Rotation k空間もハイブリッド (比較用)
base_k_space_direct = k_space_org_3D;
base_k_space_direct(hybrid_row_indices, hybrid_col_indices, :) = k_space_direct_rotated(hybrid_row_indices, hybrid_col_indices, :);

% ... (以降のステップ6, 7, 8は再構成、表示、保存) ...

fprintf('    物理シミュレーション完了。\n');

%% --- 6. アーティファクト画像の再構成 ---
fprintf('\n6. アーティファクト画像を再構成中 (ifft2)...\n');
artifact_img_3D = complex(zeros(matrix_x, matrix_y, num_slices));
direct_img_3D = complex(zeros(matrix_x, matrix_y, num_slices));

for slice_idx = 1:num_slices
    % 逆FFTによる再構成
    artifact_img_ext = ifft2(ifftshift(base_k_space_artifact(:,:,slice_idx)));
    direct_img_ext = ifft2(ifftshift(base_k_space_direct(:,:,slice_idx)));
    
    % 拡張したY次元を元に戻す (Crop)
    artifact_img_3D(:,:,slice_idx) = artifact_img_ext(:, y_start_final:y_end_final);
    direct_img_3D(:,:,slice_idx) = direct_img_ext(:, y_start_final:y_end_final);
end

clear base_k_space_artifact base_k_space_direct artifact_img_ext direct_img_ext;


%% --- 7 & 8. 表示と保存 ---
fprintf('\n7. スライス 12 の表示を開始します。\n');
fprintf('8. 結果をファイル保存します...\n');

% スライス12のみ表示
slice_idx = 12;
original_img_slice_mag = permute(iMag_4D(:,:,slice_idx, echo_idx), [2 1 3]);
artifact_img_slice = permute(artifact_img_3D(:,:,slice_idx), [2 1 3]);
direct_img_slice = permute(direct_img_3D(:,:,slice_idx), [2 1 3]);

figure('Name', sprintf('Slice %d Motion Artifact Comparison', slice_idx), 'WindowState', 'maximized');
subplot(1, 3, 1);
imshow(original_img_slice_mag, []);
title(sprintf('Original Mag (Slice %d)', slice_idx));

subplot(1, 3, 2);
imshow(abs(artifact_img_slice), []);
title('Hybrid + Motion Artifact');

subplot(1, 3, 3);
imshow(abs(direct_img_slice), []);
title('Hybrid + Direct Rotation');
drawnow;

% ファイル保存 (3D全体)
filename_base = sprintf('1124_artifact_th%.1f', theta);
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(artifact_img_3D));
save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(artifact_img_3D));

fprintf('...シミュレーションが完了しました。\n');
fprintf('結果は %s に保存されました。\n', save_path);

% -------------------------------------------------------------------
% スクリプトの最後にローカル関数を定義します
% -------------------------------------------------------------------
function save_raw_data(~, data)
    fid = fopen(filepath, 'w');
	if fid == -1
        error('ファイルが開けませんでした: %s', filepath);
	end
	fwrite(fid, data, 'double');
	fclose(fid);
end

function k_space_3d = fft2_3d_slice_by_slice(image_3d)
% 3Dデータに対してスライスごとの2D-FFT (fftshift込み) を実行
    [rows, cols, num_slices] = size(image_3d);
    k_space_3d = complex(zeros(rows, cols, num_slices));
    for slice_idx = 1:num_slices
        k_space_3d(:, :, slice_idx) = fftshift(fft2(image_3d(:, :, slice_idx)));
    end
end

function image_3d = ifft2_3d_slice_by_slice(k_space_3d)
% 3D k空間データに対してスライスごとの2D逆FFT (ifftshift込み) を実行
    [rows, cols, num_slices] = size(k_space_3d);
    image_3d = complex(zeros(rows, cols, num_slices));
    for slice_idx = 1:num_slices
        image_3d(:, :, slice_idx) = ifft2(ifftshift(k_space_3d(:, :, slice_idx)));
    end
end

%==========================================================================
% k_y ライン収集シミュレーション関数 (Matlab) - 計算負荷削減版
% 目的: 特定の k_y 行の k空間信号を二重和で計算。k_xの計算範囲を限定。
%==========================================================================

function [k_space_line_signal_partial, kx_indices_out] = simulate_ky_line_collection_by_sum(image_data, background_phase_factor, ky_line_index, kx_collect_indices)
% INPUTS:
%   image_data (Nx x Ny): スピン密度 (複素数画像) rho(x, y)
%   background_phase_factor (Nx x Ny): 背景磁場による空間位相誤差 exp(i*Delta_phi(x,y))
%   ky_line_index (scalar): 収集したい k空間の行インデックス (1 to Ny)
%   kx_collect_indices (1 x M vector): 実際に計算/収集する k_x インデックス (1 to Nx)
%
% OUTPUTS:
%   k_space_line_signal_partial (M x 1): 収集された k空間の信号ベクトル (計算された点のみ)
%   kx_indices_out (1 x M vector): 収集された k_x インデックス

    [Nx, Ny] = size(image_data);
    M = length(kx_collect_indices);
    
    % --- 1. k空間座標の定義 ---
    k_x_indices_full = (0:Nx-1) - floor(Nx/2); 
    k_y_indices = (0:Ny-1) - floor(Ny/2); 
    
    k_y_set = k_y_indices(ky_line_index); 
    
    % --- 2. 空間座標の定義 ---
    x_coords = ((0:Nx-1) - floor(Nx/2)) / Nx; 
    y_coords = ((0:Ny-1) - floor(Ny/2)) / Ny; 
    
    % --- 3. k空間信号の計算 (二重和による積分) ---
    k_space_line_signal_partial = complex(zeros(M, 1)); % ★修正: サイズを M x 1 に変更★
    
    % 画像データに背景磁場位相誤差を適用
    signal_with_B0_error = image_data .* background_phase_factor;
    
    % FE軸 (k_x) の各サンプリング点 (k_x_set) ごとに計算を実行
    % ★修正: ループを kx_collect_indices の範囲のみに限定★
    
    for idx = 1:M
        kx_line_idx = kx_collect_indices(idx); % 1-based の kx インデックス
        k_x_set = k_x_indices_full(kx_line_idx); % 対応する k_x 値
        
        S_point = 0;
        
        % --- 空間座標ループ (二重和: Sigma_x Sigma_y) ---
        for x_idx = 1:Nx
            x = x_coords(x_idx);
            for y_idx = 1:Ny
                y = y_coords(y_idx);
                
                rho_xy_error = signal_with_B0_error(x_idx, y_idx);
                
                % 位相項: exp(-i * 2*pi * (k_x * x + k_y * y))
                phase_term = exp(-1i * 2*pi * (k_x_set * x + k_y_set * y));
                
                S_point = S_point + rho_xy_error * phase_term;
            end
        end
        
        % 計算された k空間の1点 S(k_x_set, k_y_set) を結果配列に格納
        k_space_line_signal_partial(idx) = S_point;
    end
    
    % 収集された k_x インデックスを返す
    kx_indices_out = kx_collect_indices;
end
