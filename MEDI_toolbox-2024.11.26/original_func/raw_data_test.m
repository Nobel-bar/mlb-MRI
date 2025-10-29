clear variables;
close all;

% --- 1. 初期設定 ---
image_file_1 = '/Users/nori/Downloads/matlab'; % !! 要変更 !!
image_file_2 = '2_orignal_data';
load_base_path = fullfile(image_file_1, image_file_2);

% 入力ファイル名 (拡張子なし)
input_Re_name = 'Real_0ch__1_1_1_1_1_0_0_1_23_1_1_1';
input_Im_name = 'Imgn_0ch__1_1_1_1_1_0_0_1_23_1_1_1';

% サイズに関するパラメータ
orig_matrix_x = 512; % 元データのマトリクスサイズ
orig_matrix_y = 768;
cutted_matrix_x = 224; % 実際に収集されたk空間の有効データサイズ
cutted_matrix_y = 352;
final_matrix_x = 512; % 最終的に出力する画像のサイズ
final_matrix_y = 512; 
extention = 2.0/1.3;
magnification = round(orig_matrix_y * extention);

% --- 2. データの読み込み ---
fprintf('2. RAWデータを読み込んでいます...\n');
filename_input_Re = fullfile(load_base_path, input_Re_name);
filename_input_Im = fullfile(load_base_path, input_Im_name);
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
clear original_img_Re original_img_Im data_vector_re data_vector_im;

% --- 3. k空間への変換とP0補正 ---
k_space_orig = fftshift(fftn(orig_img));
[max_val, max_idx] = max(abs(k_space_orig(:)));
[kk, mm, nn] = ind2sub(size(k_space_orig), max_idx);
fprintf('k空間の最大値は座標 (%d, %d, %d) にあります。\n', kk, mm, nn);
p0_factor = k_space_orig(max_idx) / max_val;
k_space_p0 = k_space_orig / p0_factor;

% P0補正後の実空間画像
img_shifted = ifftn(ifftshift(k_space_p0));

% (k空間の1Dプロットは省略... )
% figure; 
% y_center_idx = floor(orig_matrix_y / 2) + 1;
% data_to_plot_1d = log(abs(k_space_p0(:, y_center_idx, 1)) + 1);
% plot(data_to_plot_1d);
% title('k空間の1Dプロファイル (Y軸中心, logスケール)');
% xlabel('X (周波数インデックス)');
% ylabel('信号強度 (log)');
% grid on;


% --- 3.7 1D k空間プロファイルの逆フーリエ変換 (★ここを修正) ---
figure;

% k空間のY軸中心のインデックス
y_center_idx = floor(orig_matrix_y / 2) + 1;

% (1) 青色の線: 1D k空間スライスをIFFT (＝プロジェクション)
k_space_1d_profile = k_space_p0(:, y_center_idx, 1);
k_space_1d_unshifted = ifftshift(k_space_1d_profile);
image_1d_profile_from_k = ifft(k_space_1d_unshifted);
plot(abs(image_1d_profile_from_k));

hold on; % グラフを重ね書きモードに

% (2) 赤色の線: 実空間画像をY軸方向に「合計」(＝プロジェクション)
%     img_shifted(:,:,1) は [512 x 768] の画像
%     sum(..., 2) でY軸方向(dimension 2)に合計し、[512 x 1] のベクトルにします
image_projection_from_image = sum(img_shifted(:,:,1), 2);

% (3) 赤色の破線でプロット
plot(abs(image_projection_from_image), 'r--');
hold off; 

% グラフのタイトルと軸ラベル
title('プロジェクション（投影）の比較');
xlabel('X (空間位置インデックス)');
ylabel('信号強度 (絶対値)');
legend('k空間スライスのIFFT', '実空間画像の合計 (Projection)');
grid on;
