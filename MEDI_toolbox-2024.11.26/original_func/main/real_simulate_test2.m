%==================================================================================================
% QSM RAWデータ処理/シミュレーション統合スクリプト
% 目的: 3D実空間回転とk空間ハイブリッド化によるReadout Distortionシミュレーション
%==================================================================================================

fprintf('スクリプトを開始します (MRIアーチファクトシミュレーション)\n');
clear variables;
close all;

%% 1. パラメータとファイルパス設定
% --- 処理の根幹となるパラメータとファイルパスを定義 ---
fprintf('1. パラメータを設定しています...\n');

% --- ファイルパス設定 (環境に合わせて変更してください) ---
% (パス設定は元のまま維持)
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_000 = "C:\Users\hamaguchi\Downloads\matlab\2DGE_0deg_H";
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

image_file_0 = image_file_00;
% image_file_0 = image_file_000;

load_base_path = fullfile(image_file_0, '1_data');
load_mask_path = fullfile(image_file_0, '3_output_data');
save_path = fullfile(image_file_0, '4_rolate_output_data');

mag_filename = '1st_2DGE_0deg_mag.raw'; 
phase_filename = '1st_2DGE_0deg_phase.raw'; 

if ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf('保存フォルダを作成しました: %s\n', save_path);
end

% --- 撮像・シミュレーション パラメータ ---
params = struct();
params.matrix_size = [512, 512, 23];  % 行列サイズ [Nx, Ny, Nz]
params.TE = [0.015];                  % エコー時間 [s]
params.gamma = 2 * pi * 42.57e6;      % 磁気回転比 [rad/(s*T)]
params.FOV = 0.256;                   % 視野 [m]
params.dwell_time = 2e-5;             % k空間サンプリング間隔 [s]

% --- シミュレーション制御パラメータ ---
distortion_scale_factor = 1;          % Readout歪み位相の増幅係数 (1で元の磁場不均一性)
extention = 2.0/1.3;                  % Y方向のゼロパディング拡張率
magnification = round(512 * extention); % 拡張後のYサイズ

% --- モーション（回転）パラメータ ---
theta = -18.6;                        % 回転角度 [度]
rotation_axis = [0 0 1];              % 回転軸 (Z軸)

% --- k空間ハイブリッド化パラメータ ---
cutted_matrix_x = 224; 
cutted_matrix_y = 352;
width = 112;                          % 置き換える行数
pix_start_row = 112;                  % 置き換え開始行 (k空間ROI)

%% 2. RAWデータと処理済みデータの読み込み
% --- QSMに必要な強度、位相、処理済みデータ(iFreq, RDF)をロード ---
fprintf('\n2. RAWデータと処理済みデータを読み込んでいます...\n');

mag_filepath = fullfile(load_base_path, mag_filename);
phase_filepath = fullfile(load_base_path, phase_filename);
precision = 'double=>double'; 

% 強度データ読み込み
fid_mag = fopen(mag_filepath, 'rb');
if fid_mag == -1, error('強度ファイルを開けませんでした: %s', mag_filepath); end
data_vector_mag = fread(fid_mag, inf, precision);
fclose(fid_mag);

% スライス枚数と次元の計算
pixels_per_slice = params.matrix_size(1) * params.matrix_size(2);
num_echoes = length(params.TE);
params.matrix_size(3) = numel(data_vector_mag) / (pixels_per_slice * num_echoes);
if mod(params.matrix_size(3), 1) ~= 0
    error('ファイルサイズが不正です。行列サイズ設定を確認してください。'); 
end
dims = [params.matrix_size, num_echoes];
iMag_4D = reshape(data_vector_mag, dims);
clear data_vector_mag; 

% 位相データ読み込み
fid_phase = fopen(phase_filepath, 'rb');
if fid_phase == -1, error('位相ファイルを開けませんでした: %s', phase_filepath); end
iPhase_4D = reshape(fread(fid_phase, inf, precision), dims); 
fclose(fid_phase);

% 処理済みデータ (iFreq, RDF) の読み込み
try
    load(fullfile(load_mask_path, 'phase.mat'), 'iFreq');
    load(fullfile(load_mask_path, 'PDF.mat'), 'RDF');
catch ME
    fprintf('警告: phase.mat または PDF.mat の読み込みに失敗しました。\n');
    rethrow(ME);
end

% データ準備
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);
echo_idx = 1;

original_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * iPhase_4D(:,:,:, echo_idx)); % フル画像
Highpass_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * (RDF));                         % High-pass フィルター画像
Back_img = exp(1i * (iFreq - RDF));                                                 % 背景磁場成分 (位相項)

fprintf('    データの読み込み完了。%d スライスを処理します。\n', num_slices);

%% 3. 🖼️ 実空間データの拡張 (ゼロパディング)
% --- FFT後のk空間のハイブリッド化に備え、Y軸（位相エンコード軸）を拡張 ---
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

%% 4. 3D 実空間の回転（モーションシミュレート）
% --- 回転関数 perform_padded_rotation を使用してデータを回転 ---
fprintf('\n4. 3D 実空間の回転 (imrotate3カスタム) を実行中...\n');

% perform_padded_rotationはスクリプト末尾に定義されている関数
extend_org_rotated = perform_padded_rotation(extend_org, -theta);
extend_high_rotated = perform_padded_rotation(extend_high, -theta);
extend_back_rotated = extend_back; 

fprintf('    3D 実空間の回転が完了しました。\n'); 

%% 5. 🧬 k空間ハイブリッド化と歪みシミュレーション
% --- k空間のベース作成、歪みシミュレーション、DC補正、ハイブリッド化 ---
fprintf('\n5. k空間のハイブリッド化とReadout歪みシミュレーションを実行中...\n');

% --- 5.1 回転なしデータ (Base) のk空間作成 ---
base_and_simulate_space = fft2_3d_slice_by_slice(extend_org); 
base_and_direct_space = base_and_simulate_space; % コピー

% --- 5.2 回転ありデータ (Direct: 歪みなし) のk空間作成 ---
direct_space = fft2_3d_slice_by_slice(extend_org_rotated);

% --- 5.3 回転あり + Readout Distortion (Simulation) ---
fprintf('    物理シミュレーションによるk空間を計算中...\n');
% 静的dB歪みモデル: S_distort = FT { Image_rotated * exp(i * Δφ * Scale) }


% 背景磁場によるTEでの位相誤差を抽出
Background_Phase_Map_Rotated = angle(extend_back_rotated); 

% 位相誤差項の計算 (歪み強制可視化のためスケールファクタを使用)
% Background_Phase_Map_Rotatedは ΔB * TE に相当します
Phase_Error_Term = exp(1i * Background_Phase_Map_Rotated * distortion_scale_factor);
fprintf('    位相誤差を %d 倍に拡大して適用しました。\n', distortion_scale_factor);

% 空間ドメインで位相誤差を乗算し、一度でFFT (k空間歪み) を生成
Image_with_Distortion = extend_high_rotated .* Phase_Error_Term;
simulate_space = fft2_3d_slice_by_slice(Image_with_Distortion);

fprintf('    物理シミュレーション完了。\n');

% --- 5.4 ハイブリッド化 (データの置換・DC成分補正) ---
fprintf('    k空間データをハイブリッド化 (DC成分による複素数補正) しています...\n');

x_center = floor(matrix_x / 2) + 1; 
x_start_cut = x_center - floor(cutted_matrix_x / 2); 
y_center_org = floor(magnification / 2) + 1; 

y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2); 
y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1; 

hybrid_row_indices = (x_start_cut + pix_start_row - 1) : (x_start_cut + pix_start_row + width - 2); 
hydrid_y_indices = y_start_org_cut: y_end_org_cut;
 
% DC成分による複素数補正とデータの置換
target_region_base = base_and_simulate_space(hybrid_row_indices, hydrid_y_indices, :);
target_region_sim = simulate_space(hybrid_row_indices, hydrid_y_indices, :);

center_x_local = floor(size(target_region_base, 1) / 2) + 1;
center_y_local = floor(size(target_region_base, 2) / 2) + 1;

dc_val_base = target_region_base(center_x_local, center_y_local, :);
dc_val_sim = target_region_sim(center_x_local, center_y_local, :);

complex_scale_factor = dc_val_base ./ dc_val_sim;
complex_scale_factor(isnan(complex_scale_factor) | isinf(complex_scale_factor)) = 0;

scale_map = repmat(complex_scale_factor, [size(target_region_sim, 1), size(target_region_sim, 2), 1]);
sim_data_corrected = target_region_sim .* scale_map;

% データの置換
base_and_simulate_space(hybrid_row_indices, hydrid_y_indices, :) = sim_data_corrected;
base_and_direct_space(hybrid_row_indices, hydrid_y_indices, :) = direct_space(hybrid_row_indices, hydrid_y_indices, :);

clear direct_space target_region_sim sim_data_corrected scale_map; 


%% 6. 🖼️ 再構成とトリミング
% --- k空間データをIFFTし、拡張したY軸を元に戻す ---
fprintf('\n6. アーティファクト画像を再構成中...\n');

% DC位相オフセットの最終調整
base_and_simulate_space = correct_global_phase(base_and_simulate_space);
base_and_direct_space = correct_global_phase(base_and_direct_space);

% 逆FFTによる再構成
artifact_img_ext = ifft2_3d_slice_by_slice(base_and_simulate_space);
direct_img_ext = ifft2_3d_slice_by_slice(base_and_direct_space);

% 拡張したY次元を元に戻す (Crop)
artifact_img = artifact_img_ext(:, y_start_final:y_end_final, :);
direct_img = direct_img_ext(:, y_start_final:y_end_final, :);

clear artifact_img_ext direct_img_ext; 


%% 7. 👁️ 結果の表示
% --- 特定のスライスを表示し、シミュレーション結果を比較 ---
fprintf('\n7. スライスごとの表示を開始します...\n');

for slice_idx = 12:12 
    fprintf('    --- スライス %d / %d を処理中 ---\n', slice_idx, num_slices);
    
    direct_img_slice = permute(direct_img(:,:,slice_idx), [2 1]); 
    artifact_img_slice = permute(artifact_img(:,:,slice_idx), [2 1]);
    
    figure('Name', sprintf('Slice %d Comparison', slice_idx), 'WindowState', 'maximized');
    subplot(1, 2, 1);
    imshow(abs(direct_img_slice), []);
    title(sprintf('Direct Rotation (Slice %d)', slice_idx));
    
    subplot(1, 2, 2);
    imshow(abs(artifact_img_slice), []);
    title(sprintf('Hybrid + Distortion (Slice %d)', slice_idx));
    drawnow;
end

%% 8. 💾 データ保存
% --- 再構成された複素数画像データを保存 ---
fprintf('\n8. 結果をファイル保存します...\n');
filename_base = sprintf('real_2D_rotate_artifact_th%.1f', theta);

save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(artifact_img));

fprintf('結果は %s に保存されました。\n', save_path);


% ===================================================================
% 🛠️ Local Functions (補助関数)
% ===================================================================

function save_raw_data(filepath, data)
    % RAWデータをダブル精度で保存
    fid = fopen(filepath, 'w');
    if fid == -1, error('ファイルが開けませんでした: %s', filepath); end
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

function corrected_kspace = correct_global_phase(kspace_3d)
% グローバル位相の補正: 最大強度の点の位相を0度にする
    [max_val, max_idx] = max(abs(kspace_3d(:)));
    if max_val == 0
        corrected_kspace = kspace_3d;
        return;
    end
    p0_factor = kspace_3d(max_idx) / max_val;
    corrected_kspace = kspace_3d / p0_factor;
end

function vol_out = perform_padded_rotation(vol_in, theta)
% 3Dボリュームの各スライスに、パディング -> 回転 -> 中央切り出しを行うカスタム回転関数
    [rows, cols, num_slices] = size(vol_in);
    vol_out = complex(zeros(rows, cols, num_slices));
    
    % 座標グリッドの作成とパディングサイズ計算 (元のロジックを維持)
    pad_rows_total = rows * 3;
    pad_cols_total = cols * 3;
    [X_in, Y_in] = meshgrid(1:pad_cols_total, 1:pad_rows_total);
    centerX = (pad_cols_total + 1) / 2;
    centerY = (pad_rows_total + 1) / 2;
    
    X_shifted = X_in - centerX;
    Y_shifted = Y_in - centerY;
    
    theta_rad = theta * (pi/180);
    cosT = cos(theta_rad);
    sinT = sin(theta_rad);
    
    % 逆回転座標の計算
    X_orig = X_shifted * cosT + Y_shifted * sinT + centerX;
    Y_orig = -X_shifted * sinT + Y_shifted * cosT + centerY;
    
    % 画像を埋め込む位置 (中央)
    start_r = floor((pad_rows_total - rows)/2) + 1;
    end_r = start_r + rows - 1;
    start_c = floor((pad_cols_total - cols)/2) + 1;
    end_c = start_c + cols - 1;
    
    % スライスごとの処理ループ
    for z = 1:num_slices
        slice_data = vol_in(:,:,z);
        
        % パディングと配置
        padded_img = complex(zeros(pad_rows_total, pad_cols_total));
        padded_img(start_r:end_r, start_c:end_c) = slice_data;
        
        % interp2 による回転 (実部・虚部に分離して処理)
        rotated_img_Re = interp2(X_in, Y_in, real(padded_img), X_orig, Y_orig, 'linear', 0);
        rotated_img_Im = interp2(X_in, Y_in, imag(padded_img), X_orig, Y_orig, 'linear', 0);
        
        rotated_padded_img = complex(rotated_img_Re, rotated_img_Im);
        
        % 中央切り出し
        vol_out(:,:,z) = rotated_padded_img(start_r:end_r, start_c:end_c);
    end
end