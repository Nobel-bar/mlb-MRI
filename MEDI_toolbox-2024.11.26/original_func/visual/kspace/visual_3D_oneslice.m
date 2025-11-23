%==================================================================================================
% 実際のMRI RAWデータ (3D) を読み込み、P0補正を行い、
% [★修正★] k空間の *スライス12* の振幅(A)を3Dサーフェスプロットで可視化するスクリプト
% [★修正★] 1つのFigureに線形スケールと対数スケールの2つを表示
%% [★変更★] k空間 スライス12 の (X=位相, Y=ky, Z=振幅) の関係を
% [★変更★] 3D散布図 (scatter3) で可視化するスクリプト
% k空間の A vs. (X=位相, Y=ky, Z=振幅) を描画
% 2Dスライスデータも [512, 768] と巨大なため、
% surf で描画するために 8x8 ごとにサンプリング（間引き）します。
%==================================================================================================

fprintf('スクリプトを開始します (実際のMRI k空間 3Dサーフェス描画)\n');
clear variables;
close all;

%% --- 1. 初期設定 ---
fprintf('1. パラメータを設定しています...\n');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_file_0 = image_file_2DGE_1_2_Rotate_H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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
% [★修正★] fftn の直後に fftshift を追加し、k空間の中心を中央に移動
k_space_orig_shifted = fftn(orig_img); 

% P0補正 (fftshift "前" のデータで最大値を探す)
k_space_orig_unsh = fftn(orig_img);
[max_val, max_idx] = max(abs(k_space_orig_unsh(:)));
[kk, mm, nn] = ind2sub(size(k_space_orig_unsh), max_idx);
fprintf('k空間の最大値は座標 (%d, %d, %d) にあります。\n', kk, mm, nn);
p0_factor = k_space_orig_unsh(max_idx) / max_val;

% P0補正を「fftshift 済み」の k空間に適用
k_space_p0_shifted = k_space_orig_shifted / p0_factor; 

% 振幅(A)を計算
A_values = abs(k_space_p0_shifted);

%% --- 4. [★修正★] スライス12のデータを抽出・サンプリング ---
fprintf('4. k空間のスライス12を抽出し、サンプリング中...\n');

% k空間のサイズと中心オフセットを取得
matrix_size_k = size(A_values); % [512, 768, 23]
center_offset_i = floor(matrix_size_k(1) / 2) + 1; % 257 (ky)
center_offset_j = floor(matrix_size_k(2) / 2) + 1; % 385 (kx)
center_offset_k = floor(matrix_size_k(3) / 2) + 1; % 12 (kz)

% ユーザー指定の「スライス12」 (中心スライス) を抽出
% (Slice が 23 の場合、center_offset_k は 12 になります)
slice_to_plot = center_offset_k;
A_slice = A_values(:, :, slice_to_plot);

% kx, ky の座標軸を作成 (例: -384 ... +383)
ky_coords_vec = (1:matrix_size_k(1)) - center_offset_i;
kx_coords_vec = (1:matrix_size_k(2)) - center_offset_j;

% [★サンプリング★] 
% 512x768 のメッシュは重すぎるため、8x8ピクセルごとに間引く
sample_step = 2;
ky_coords_sampled = ky_coords_vec(1:sample_step:end);
kx_coords_sampled = kx_coords_vec(1:sample_step:end);
A_slice_sampled = A_slice(1:sample_step:end, 1:sample_step:end);

% surf のためのグリッドを作成 (meshgrid は [X, Y] の順)
[KX, KY] = meshgrid(kx_coords_sampled, ky_coords_sampled);

% A_slice_sampled は (64, 96)
% KX, KY も (64, 96) となり、surf(KX, KY, Z) の Z とサイズが一致
A_slice_sampled_for_surf = A_slice_sampled;


%% --- 5. [★修正★] 3Dサーフェスプロット (線形 vs 対数) ---
fprintf('5. 3Dサーフェスプロット (surf) を描画中...\n');
fig = figure('Name', 'Actual MRI 3D k-space Surface (Slice 12)', 'WindowState', 'maximized');

% 1x2 のグリッドレイアウトを作成
t = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact');
title(t, sprintf('k空間 振幅 (スライス %d, %dx%d サンプリング)', slice_to_plot, sample_step, sample_step));

% --- グラフ1: 線形スケール (Linear Scale) ---
ax1 = nexttile;
surf(ax1, KX, KY, A_slice_sampled_for_surf, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
shading(ax1, 'interp'); % 滑らかにシェーディング

% グラフの設定
title(ax1, '線形スケール (Linear Scale)');
xlabel(ax1, 'kx (位相)');
ylabel(ax1, 'ky (周波数)');
zlabel(ax1, '振幅 (A)');
axis(ax1, 'tight');
grid(ax1, 'on');
box(ax1, 'on'); 

% カラーバー
h_bar1 = colorbar(ax1);
ylabel(h_bar1, '振幅 (A)');
colormap(ax1, 'turbo'); % [★修正★] parula から jet に変更
view(ax1, 3); 
rotate3d(ax1, 'on');


% --- グラフ2: 対数スケール (Log Scale) ---
ax2 = nexttile;

% [★修正★] Z軸のデータ (ZData) と 色のデータ (CData) を準備
ZData_log = A_slice_sampled_for_surf + eps; % Z軸の形状データ (線形)
CData_log = log10(A_slice_sampled_for_surf + eps); % 色のデータ (対数)

% [★修正★] surf(X, Y, Z, C) を使用し、色データ(C)に CData_log を指定
surf(ax2, KX, KY, ZData_log, CData_log, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
shading(ax2, 'interp'); 

% [★修正★] Z軸を対数スケールに設定 (形状)
set(ax2, 'ZScale', 'log');

% グラフの設定
title(ax2, '対数スケール (Log Scale) - 色も対数');
xlabel(ax2, 'kx (位相)');
ylabel(ax2, 'ky (周波数)');
zlabel(ax2, '振幅 (A) - log scale');
axis(ax2, 'tight');
grid(ax2, 'on');
box(ax2, 'on'); 

% カラーバー
h_bar2 = colorbar(ax2);
% [★修正★] カラーバーのラベルも対数を反映
ylabel(h_bar2, 'log10(振幅)'); 
colormap(ax2, 'turbo'); 
view(ax2, 3); 
rotate3d(ax2, 'on');


fprintf('完了。マウスでグラフを回転させて確認してください。\n');