%==================================================================================================
% 実際のMRI RAWデータ (3D) を読み込み、P0補正を行い、
% k空間の振幅(A)を3Dで可視化するスクリプト

% k空間の A vs. Phase グラフ（Y軸対数）を描画
% 3D表示
%==================================================================================================
%==================================================================================================

fprintf('スクリプトを開始します (実際のMRI k空間 3D描画)\n');
clear variables;
close all;

%% --- 1. 初期設定 ---
fprintf('1. パラメータを設定しています...\n');

image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_000 = "C:\Users\hamaguchi\Downloads\matlab";
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_file_0 = image_file_00; % 
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
k_space_orig_shifted = fftshift(fftn(orig_img)); 

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

%% --- 4. プロットするデータの選別 (★重要: サンプリング) ---
fprintf('4. k空間の強い信号を選別・サンプリング中...\n');

% 振幅が最大値の 5% 未満の点はノイズとして無視する (5%に引き上げ)
threshold_level = 0.05 * max(A_values(:)); 
indices_to_plot = find(A_values >= threshold_level);

% プロットする点の最大数を 500,000 点に制限
max_plot_points = 500000;
if numel(indices_to_plot) > max_plot_points
    fprintf('  %d 点が閾値を超えましたが、%d 点にランダムサンプリングします。\n', ...
            numel(indices_to_plot), max_plot_points);
    % datasample を使ってインデックスをランダムに間引く
    indices_to_plot = datasample(indices_to_plot, max_plot_points, 'Replace', false);
else
    fprintf('  %d 点をプロットします。\n', numel(indices_to_plot));
end

% 選別した点のインデックス (i,j,k) と振幅 A を取得
A_plot = A_values(indices_to_plot);

% インデックス (i,j,k) を kx, ky, kz 座標に変換
matrix_size_k = size(A_values); % [512, 768, 23]
center_offset_i = floor(matrix_size_k(1) / 2) + 1; % 257 (ky)
center_offset_j = floor(matrix_size_k(2) / 2) + 1; % 385 (kx)
center_offset_k = floor(matrix_size_k(3) / 2) + 1; % 12 (kz)

[i, j, k] = ind2sub(matrix_size_k, indices_to_plot);

% i が ky, j が kx, k が kz に対応 (fftshift 後の座標系)
ky_coords = i - center_offset_i;
kx_coords = j - center_offset_j;
kz_coords = k - center_offset_k;

%% --- 5. 3D k空間のプロット (scatter3) ---
fprintf('5. 3Dプロット (scatter3) を描画中...\n');
figure('Name', 'Actual MRI 3D k-space Visualization', 'WindowState', 'maximized');

% S (サイズ): 振幅(A)が小さい点は小さく
% C (色): 振幅(A)が小さい点は暗く
scatter_size = (A_plot / max(A_plot)) * 50 + 5; % サイズを調整
scatter_color = A_plot; % 振幅をそのまま色データとして使用

scatter3(kx_coords, ky_coords, kz_coords, scatter_size, scatter_color, ...
         'filled', 'MarkerFaceAlpha', 0.2); % Alphaを下げて透明度を上げる

% --- グラフの設定 ---
title('実際のMRI 3D k空間データの可視化 (振幅)');
xlabel('kx (位相)');
ylabel('ky (周波数)');
zlabel('kz (スライス)');

axis equal; % 各軸のスケールを等しく
grid on;
box on; 

% カラーバー (色のスケール) を表示
h_bar = colorbar;
ylabel(h_bar, '振幅 (A)');
colormap(parula); % または 'jet'

% 視点を調整
view(3); 
rotate3d on; % マウスでグリグリ回せるようにする

fprintf('完了。マウスでグラフを回転させて確認してください。\n');