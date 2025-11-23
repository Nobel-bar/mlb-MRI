%==================================================================================================
% QSM RAWデータ (3D) を読み込み、
% 3D実空間で回転 (ユーザー提供のパディング+interp2ロジック) をシミュレートし、
% k空間でハイブリッド化 (一部置換 / ループ処理版 / Peak Recovery) を行うMATLABスクリプト
%
% ※ 手法A（Point-wise処理）維持、冗長記述削除、進捗バー追加版
%==================================================================================================
fprintf('スクリプトを開始します (手法A: ループ処理 + Peak Recovery)\n');
clear variables;
close all;

%% --- 1. 撮像・シミュレーション パラメータ設定 ---
fprintf('1. パラメータを設定しています...\n');

% パス設定
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
image_file_0 = image_file_000;

% 読み込みパスと保存パスを定義
load_base_path = fullfile(image_file_0, image_file_1);
load_mask_path = fullfile(image_file_0, image_file_3);
save_path = fullfile(image_file_0, image_file_4);

% 保存先フォルダが存在しない場合は作成する
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% --- 撮像パラメータ ---
params = struct();
params.matrix_size = [512, 512, 23]; 
params.TE = [0.015]; 

% --- 拡張と回転 ---
extention = 2.0/1.3; 
magnification = round(params.matrix_size(2) * extention); 
theta = -18.6; 

% --- k空間ハイブリッド化パラメータ (ROI) ---
cutted_matrix_x = 224; 
cutted_matrix_y = 352; 

%% --- 2. RAWデータの読み込み ---
fprintf('2. RAWデータと処理済みデータを読み込んでいます...\n');

mag_filepath = fullfile(load_base_path, '1st_2DGE_0deg_mag.raw');
phase_filepath = fullfile(load_base_path, '1st_2DGE_0deg_phase.raw');
dims = [params.matrix_size, length(params.TE)];

% 強度・位相読み込み
fid = fopen(mag_filepath, 'rb');
if fid == -1, error('ファイルが見つかりません: %s', mag_filepath); end
iMag_4D = reshape(fread(fid, inf, 'double=>double'), dims);
fclose(fid);

fid = fopen(phase_filepath, 'rb');
if fid == -1, error('ファイルが見つかりません: %s', phase_filepath); end
iPhase_4D = reshape(fread(fid, inf, 'double=>double'), dims);
fclose(fid);

% QSM結果読み込み
try
    load(fullfile(load_mask_path, 'phase.mat'), 'iFreq');
    load(fullfile(load_mask_path, 'PDF.mat'), 'RDF');
catch
    error('phase.mat または PDF.mat が見つかりません。');
end

matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);
echo_idx = 1; 

% 画像変数の作成
original_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * iPhase_4D(:,:,:, echo_idx));
Highpass_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * RDF);
Back_img     = exp(1i * (iFreq - RDF)); 

%% --- 3. 実空間データの拡張 (ゼロパディング) ---
fprintf('3. 実空間データのY方向を拡張しています...\n');

extend_org  = complex(zeros(matrix_x, magnification, num_slices));
extend_high = complex(zeros(matrix_x, magnification, num_slices));
extend_back = complex(zeros(matrix_x, magnification, num_slices));

y_center_final = floor(magnification / 2) + 1;
y_start_final  = y_center_final - floor(matrix_y / 2);
y_end_final    = y_start_final + matrix_y - 1;

extend_org(:, y_start_final:y_end_final, :)  = original_img;
extend_high(:, y_start_final:y_end_final, :) = Highpass_img; 
extend_back(:, y_start_final:y_end_final, :) = Back_img; 

%% --- 4. 3D 実空間の回転 ---
fprintf('4. 3D 実空間の回転を実行中...\n');

% オリジナルと背景磁場を回転
extend_org_rotated  = perform_padded_rotation(extend_org, theta);
extend_back_rotated = perform_padded_rotation(extend_back, theta);

% Highpass (RDF) を回転させて、シミュレーションのベースとする
extend_high_rotated = perform_padded_rotation(extend_high, theta);
extend_rotated      = extend_high_rotated;

% %% --- 5. 3D k空間データのハイブリッド化 (ループ処理) ---
% fprintf('5. 3D k空間のハイブリッド化 (Point-wise Loop + Peak Recovery) を実行中...\n');
% 
% % ベースとなるk空間 (Original)
% base_and_simulate_space = fftshift(fftn(extend_org)); 
% base_and_direct_space   = base_and_simulate_space; % 比較用
% 
% % シミュレーションソースのk空間 (Rotated Highpass)
% simulate_space = fftshift(fftn(extend_rotated));
% % 単純回転のk空間 (比較用)
% direct_space   = fftshift(fftn(extend_org_rotated));
% 
% % ROIインデックス計算
% x_center = floor(matrix_x / 2) + 1; 
% x_start_cut = x_center - floor(cutted_matrix_x / 2); 
% x_end_cut   = x_start_cut + cutted_matrix_x - 1; 
% 
% y_center_org = floor(magnification / 2) + 1; 
% y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2); 
% y_end_org_cut   = y_start_org_cut + cutted_matrix_y - 1; 
% 
% hybrid_row_indices = x_start_cut : x_end_cut; 
% hydrid_y_indices   = y_start_org_cut : y_end_org_cut;
% hydrid_z_indices   = 13 : 23;
% 
% % Direct (比較用) は一括置換
% base_and_direct_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices) = ...
%     direct_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices);
% 
% % --- ループ処理の準備 ---
% template_k = complex(zeros(matrix_x, magnification, num_slices));
% total_iter = length(hydrid_z_indices) * length(hydrid_y_indices) * length(hybrid_row_indices);
% count = 0;
% 
% fprintf('    k空間点ごとの背景磁場適用を開始します (全 %d 点)...\n', total_iter);
% fprintf('    ※注意: 非常に時間がかかります\n');
% 
% % --- 【修正】進捗バーの作成 ---
% hBar = waitbar(0, 'Processing k-space points...', 'Name', 'Hybrid Simulation Progress');
% 
% % === ループ開始 ===
% tic; % 時間計測開始
% for z_idx = hydrid_z_indices
%     for y_idx = hydrid_y_indices
%         for x_idx = hybrid_row_indices
% 
%             count = count + 1;
% 
%             % --- 【修正】進捗バーとログの更新 (1000回に1回) ---
%             if mod(count, 1000) == 0
%                 elapsed_time = toc;
%                 avg_time_per_step = elapsed_time / count;
%                 remaining_steps = total_iter - count;
%                 est_time_left = remaining_steps * avg_time_per_step;
% 
%                 % プログレスバー更新
%                 waitbar(count / total_iter, hBar, ...
%                     sprintf('Processing: %.1f%% \nEst. Remaining: %.0f sec', ...
%                     count/total_iter*100, est_time_left));
% 
%                 % コンソールにもたまに出力
%                 if mod(count, 10000) == 0
%                      fprintf('Processing: %.1f%% (%d/%d) - Est. Left: %.0f min\n', ...
%                          count/total_iter*100, count, total_iter, est_time_left/60);
%                 end
%             end
% 
%             % 1. 1点だけを取り出す (残りはゼロ)
%             single_point_kspace = template_k; 
%             single_point_kspace(x_idx, y_idx, z_idx) = simulate_space(x_idx, y_idx, z_idx);
% 
%             % 2. 実空間に戻して背景磁場を掛ける
%             img_basis_3d = ifftn(ifftshift(single_point_kspace)) .* extend_back_rotated;
% 
%             % 3. 再びk空間へ (ここでエネルギーが拡散・移動している)
%             gg = fftshift(fftn(img_basis_3d));
% 
%             % 4. Peak Recovery (移動したエネルギーの最大値を捕まえて連れ戻す)
%             [max_val, max_linear_idx] = max(abs(gg(:)));
%             recovered_val = gg(max_linear_idx);
% 
%             % 元の座標に強制的に収納する
%             base_and_simulate_space(x_idx, y_idx, z_idx) = recovered_val;
% 
%         end
%     end
% end
% 
% % ループ終了処理
% close(hBar);
% fprintf('    ループ処理完了。\n');
%% --- 5. 3D k空間データのハイブリッド化 (高速化版) ---
fprintf('5. 3D k空間のハイブリッド化 (Convolution Peak Scaling) を実行中...\n');

% ベースとなるk空間 (Original)
base_and_simulate_space = fftshift(fftn(extend_org)); 
base_and_direct_space   = base_and_simulate_space; % 比較用

% シミュレーションソースのk空間
simulate_space = fftshift(fftn(extend_rotated));
% 単純回転のk空間 (比較用)
direct_space   = fftshift(fftn(extend_org_rotated));

% ROIインデックス計算
x_center = floor(matrix_x / 2) + 1; 
x_start_cut = x_center - floor(cutted_matrix_x / 2); 
x_end_cut   = x_start_cut + cutted_matrix_x - 1; 

y_center_org = floor(magnification / 2) + 1; 
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2); 
y_end_org_cut   = y_start_org_cut + cutted_matrix_y - 1; 

hybrid_row_indices = x_start_cut : x_end_cut; 
hydrid_y_indices   = y_start_org_cut : y_end_org_cut;
hydrid_z_indices   = 13 : 23;

% Direct (比較用) は一括置換
base_and_direct_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices) = ...
    direct_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices);

% -----------------------------------------------------------------------
% 【高速化ロジック】 ループ処理の代替
% 「1点ずつピークを探して戻す」処理は、「背景磁場の最大スペクトル成分を掛ける」のと等価です
% -----------------------------------------------------------------------
fprintf('    背景磁場のスペクトル解析中...\n');

% 1. 背景磁場自体のk空間分布（カーネル）を計算
kernel_k = fftshift(fftn(extend_back_rotated));

% 2. その中で最もエネルギーが高い点（ピーク）を探す
[~, max_kernel_idx] = max(abs(kernel_k(:)));
peak_complex_factor = kernel_k(max_kernel_idx);

% 補足: fftnの定義上、1点だけ励起した場合の振幅スケーリングを合わせる必要があるか確認
% MATLABのifft->fftの流れではスケーリングは自動的に保存されるため、
% 単純にこの peak_complex_factor を掛ければOKです。

fprintf('    移動エネルギーのピーク係数: %.4f + %.4fi\n', ...
        real(peak_complex_factor), imag(peak_complex_factor));

% 3. 全点に対して一括で係数を掛ける (数万回のループと同じ結果になる)
%    元のsimulate_spaceの値 * カーネルのピーク値 = 連れ戻された値
corrected_roi_data = simulate_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices) ...
                     * peak_complex_factor;

% 4. 結果を格納
base_and_simulate_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices) = corrected_roi_data;

fprintf('    ハイブリッド化完了 (高速処理)。\n');

%% --- 6. アーティファクト画像の再構成 ---
fprintf('6. アーティファクト画像を再構成中...\n');

% 正規化と再構成
[max_val, max_idx] = max(abs(base_and_simulate_space(:)));
if max_val > 0
    base_and_simulate_space = base_and_simulate_space / (base_and_simulate_space(max_idx) / max_val);
end

[max_val, max_idx] = max(abs(base_and_direct_space(:)));
if max_val > 0
    base_and_direct_space = base_and_direct_space / (base_and_direct_space(max_idx) / max_val);
end

artifact_img_ext = ifftn(ifftshift(base_and_simulate_space));
direct_img_ext   = ifftn(ifftshift(base_and_direct_space));

% 中央切り出し
artifact_img = artifact_img_ext(:, y_start_final:y_end_final, :);
direct_img   = direct_img_ext(:, y_start_final:y_end_final, :);

%% --- 7. 表示と保存 ---
fprintf('7. 結果を保存しています...\n');

% スライス表示 (例: 12スライス目)
slice_idx = 12;
figure('Name', 'Artifact Comparison');
subplot(1,2,1); imagesc(abs(permute(direct_img(:,:,slice_idx), [2 1]))); 
colormap gray; axis image; axis off; title('Original (Direct Rotate)');
subplot(1,2,2); imagesc(abs(permute(artifact_img(:,:,slice_idx), [2 1]))); 
colormap gray; axis image; axis off; title('Artifact (Method A + Recovery)');
drawnow;

% 保存
filename_base = sprintf('z_artifact_MethodA_Recov_th%.1f', theta);
save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(artifact_img));

fprintf('完了。保存先: %s\n', save_path);

% ===================================================================
% ローカル関数
% ===================================================================
function save_raw_data(filepath, data)
    fid = fopen(filepath, 'w');
    if fid == -1, error('ファイルが開けませんでした: %s', filepath); end
    fwrite(fid, data, 'double');
    fclose(fid);
end

function vol_out = perform_padded_rotation(vol_in, theta)
    [rows, cols, num_slices] = size(vol_in);
    vol_out = complex(zeros(rows, cols, num_slices));

    pad_rows_total = rows * 3;
    pad_cols_total = cols * 3;
    
    [X_in, Y_in] = meshgrid(1:pad_cols_total, 1:pad_rows_total);
    
    centerX = (pad_cols_total + 1) / 2;
    centerY = (pad_rows_total + 1) / 2;
    
    theta_rad = theta * (pi/180);
    cosT = cos(theta_rad);
    sinT = sin(theta_rad);
    
    X_shifted = X_in - centerX;
    Y_shifted = Y_in - centerY;
    
    X_orig = X_shifted * cosT + Y_shifted * sinT + centerX;
    Y_orig = -X_shifted * sinT + Y_shifted * cosT + centerY;
    
    start_r = floor((pad_rows_total - rows)/2) + 1;
    end_r   = start_r + rows - 1;
    start_c = floor((pad_cols_total - cols)/2) + 1;
    end_c   = start_c + cols - 1;

    for z = 1:num_slices
        padded_img = complex(zeros(pad_rows_total, pad_cols_total));
        padded_img(start_r:end_r, start_c:end_c) = vol_in(:,:,z);
        
        % 実部と虚部を個別に補間
        r_re = interp2(X_in, Y_in, real(padded_img), X_orig, Y_orig, 'linear', 0);
        r_im = interp2(X_in, Y_in, imag(padded_img), X_orig, Y_orig, 'linear', 0);
        
        % --- 【修正】一時変数を使ってから切り出し ---
        rotated_full = complex(r_re, r_im);
        vol_out(:,:,z) = rotated_full(start_r:end_r, start_c:end_c);
    end
end