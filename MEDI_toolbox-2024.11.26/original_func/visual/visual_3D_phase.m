%==================================================================================================
% 実際のMRI RAWデータ (3D) を読み込み、P0補正を行い、
% [★変更★] k空間 スライス12 の (X=位相, Y=ky, Z=振幅) の関係を
% [★変更★] 3D散布図 (scatter3) で可視化するスクリプト
%
% [★注意★]
% 512x768x23 の全点は多すぎるため、k空間スライス12 のみを
% 8x8 ごとにサンプリング（間引き）してプロットします。
%==================================================================================================
fprintf('スクリプトを開始します (k空間 3D散布図: X=位相, Z=振幅)\n');
clear variables;
close all;

%% --- 1. 初期設定 ---
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H_local = 'C:\Users\hamaguchi\Downloads\matlab\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_000 = "C:\Users\hamaguchi\Downloads\matlab\2DGE_0deg_H'";
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

image_file_0 = image_file_00;
image_file_0 = image_file_000;
image_file_0 = image_file_2DGE_1_2_Rotate_H_local;


% 読み込みパスを定義
load_base_path = fullfile(image_file_0, image_file_2);
% 入力ファイル名 (拡張子なし)
input_Re_name = 'Real_0ch__1_1_1_1_1_0_0_1_23_1_1_1';
input_Im_name = 'Imgn_0ch__1_1_1_1_1_0_0_1_23_1_1_1';

% サイズに関するパラメータ
orig_matrix_x = 512;
orig_matrix_y = 768;

%% --- 2. データの読み込み ---
fprintf('2. RAWデータを読み込んでいます...\n');
filename_input_Re = fullfile(load_base_path, input_Re_name);
filename_input_Im = fullfile(load_base_path,input_Im_name);

% ベクトルとして読み込み、3D配列に変換
fileID_Re = fopen(filename_input_Re, 'r');
if fileID_Re == -1, error('ファイルが開けませんでした: %s', filename_input_Re); end
data_vector_re = fread(fileID_Re, inf, 'single');
fclose(fileID_Re);
Slice = numel(data_vector_re) / (orig_matrix_x * orig_matrix_y);
if mod(Slice, 1) ~= 0, error('実数部のファイルサイズが不正です。'); end
original_img_Re = reshape(data_vector_re, [orig_matrix_x, orig_matrix_y, Slice]);

fileID_Im = fopen(filename_input_Im, 'r');
if fileID_Im == -1, error('ファイルが開けませんでした: %s', filename_input_Im); end
data_vector_im = fread(fileID_Im, inf, 'single');
fclose(fileID_Im);
original_img_Im = reshape(data_vector_im, [orig_matrix_x, orig_matrix_y, Slice]);
orig_img = complex(original_img_Re, original_img_Im);
fprintf('%d x %d x %d の画像を正常に読み込みました。\n', orig_matrix_x, orig_matrix_y, Slice);

%% --- 3. k空間への変換とP0補正 ---
fprintf('3. 3D FFT と P0補正 を実行中...\n');
% P0補正 (fftshift "前" のデータで最大値を探す)
k_space_orig_unsh = fftn(orig_img);
[max_val, max_idx] = max(abs(k_space_orig_unsh(:)));
[kk, mm, nn] = ind2sub(size(k_space_orig_unsh), max_idx);
fprintf('k空間の最大値は座標 (%d, %d, %d) にあります。\n', kk, mm, nn);
p0_factor = k_space_orig_unsh(max_idx) / max_val;

% P0補正を適用し、その後 fftshift を実行
k_space_p0_unsh = k_space_orig_unsh / p0_factor;
k_space_p0_shifted = fftshift(k_space_p0_unsh); % P0補正済みの k空間 (中心化)

% [★変更★] 振幅(A) と 位相(Phase) の両方を計算
A_values = abs(k_space_p0_shifted);
Phase_values = angle(k_space_p0_shifted);

%% --- 4. [★変更★] スライス12の 振幅・位相・ky座標 を抽出・サンプリング ---
fprintf('4. k空間のスライス12 (振幅・位相) を抽出し、サンプリング中...\n');

% k空間のサイズと中心オフセットを取得
matrix_size_k = size(A_values); % [512, 768, 23]
center_offset_i = floor(matrix_size_k(1) / 2) + 1; % 257 (ky)
center_offset_j = floor(matrix_size_k(2) / 2) + 1; % 385 (kx)
center_offset_k = floor(matrix_size_k(3) / 2) + 1; % 12 (kz)

% ユーザー指定の「スライス12」 (中心スライス) を抽出
slice_to_plot = center_offset_k;

% [★変更★] 振幅(A) と 位相(Phase) データを抽出
A_slice = A_values(:, :, slice_to_plot);
Phase_slice = Phase_values(:, :, slice_to_plot);

% kx, ky の座標軸を作成
ky_coords_vec = (1:matrix_size_k(1)) - center_offset_i;
kx_coords_vec = (1:matrix_size_k(2)) - center_offset_j;

% [★サンプリング★] 
% 8x8ピクセルごとに間引く
sample_step = 8;
ky_coords_sampled = ky_coords_vec(1:sample_step:end);
kx_coords_sampled = kx_coords_vec(1:sample_step:end);

% [★変更★] 振幅と位相をサンプリング
A_slice_sampled = A_slice(1:sample_step:end, 1:sample_step:end);
Phase_slice_sampled = Phase_slice(1:sample_step:end, 1:sample_step:end);

% [★変更★] Y軸用の ky グリッドを作成
[~, KY_grid_sampled] = meshgrid(kx_coords_sampled, ky_coords_sampled);

% [★変更★] scatter3 のためにデータを 1次元ベクトルに変換
X_data_vec = Phase_slice_sampled(:); % X軸 = 位相
Y_data_vec = KY_grid_sampled(:);     % Y軸 = ky座標
Z_data_vec = A_slice_sampled(:);     % Z軸 = 振幅 (強度)

%% --- 5. [★変更★] 3D散布図 (scatter3) プロット ---
fprintf('5. 3D散布図 (scatter3) を描画中...\n');
fig = figure('Name', 'MRI k-space (X=Phase, Y=ky, Z=Amplitude)', 'WindowState', 'maximized');
t = tiledlayout(fig, 1, 1, 'TileSpacing', 'compact');
title(t, sprintf('k空間 (スライス %d, %dx%d サンプリング)', slice_to_plot, sample_step, sample_step));

ax1 = nexttile;

% [★変更★] scatter3(X, Y, Z, Size, Color)
% X=位相, Y=ky, Z=振幅
% 点の「色」も Z=振幅 にマッピング (CData)
scatter3(ax1, X_data_vec, Y_data_vec, Z_data_vec, ...
    10, ... % 点のサイズ
    Z_data_vec, ... % 点の色
    'filled'); 

% グラフの設定
title(ax1, 'k空間データの関係 (X=位相, Y=ky, Z=振幅)');

% [★変更★] 軸ラベルをご要望に合わせて設定
xlabel(ax1, '位相 (Phase) [rad]');
ylabel(ax1, 'ky (周波数エンコード方向)');
zlabel(ax1, '振幅 (A) [強度]');

% [★変更★] X軸の範囲を -pi から +pi に設定
xlim(ax1, [-pi, pi]);

axis(ax1, 'tight');
grid(ax1, 'on');
box(ax1, 'on'); 

% カラーバー
h_bar1 = colorbar(ax1);
ylabel(h_bar1, '振幅 (A)');

% [★変更★] Z軸(振幅)は非常にダイナミックレンジが広いため、
% Z軸とカラーマップを「対数スケール」にすると見やすくなります
set(ax1, 'ZScale', 'log');
set(ax1, 'ColorScale', 'log');
ylabel(h_bar1, '振幅 (A) - log scale');

colormap(ax1, 'turbo'); 
view(ax1, 3); 
rotate3d(ax1, 'on');

fprintf('完了。マウスでグラフを回転させて確認してください。\n');
fprintf('Z軸とカラーは対数スケールで表示しています。\n');