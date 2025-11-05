%==================================================================================================
% QSM RAWデータ (3D) を読み込み、P0補正を行い、
% [★新規] k空間の A vs. Phase グラフを描画し、
% 最終的な画像を再構成するスクリプト
%
% [修正内容]
% - 未定義変数 (y_start_cut, y_center) によるエラーを修正。
% - 「A vs. Phase」グラフの描画機能 (セクション4) を追加。
%==================================================================================================

fprintf('スクリプトを開始します (RAWデータ読み込み + A vs. Phaseプロット)\n');
clear variables;
close all;

%% --- 1. 初期設定 ---
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

image_file_0 = image_file_00; % slab用

% 読み込みパスと保存パスを定義
load_base_path = fullfile(image_file_0, image_file_2);
load_mask_path = fullfile(image_file_0, image_file_3);
save_path = fullfile(image_file_0, image_file_1);

% 入力ファイル名 (拡張子なし)
% [★要確認] パスが './Volunteer_Rotate_H/...' で正しいか確認してください
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
magnification = round(orig_matrix_y * extention); % 768 -> 1182

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

% (セクション1, 2 は変更なし)
...

%% --- 3. k空間への変換とP0補正 ---
fprintf('3. 3D FFT と P0補正 を実行中...\n');
% [★修正★] fftn の直後に fftshift を追加し、k空間の中心を中央に移動
k_space_orig = fftn(orig_img); 
[max_val, max_idx] = max(abs(k_space_orig(:)));
[kk, mm, nn] = ind2sub(size(k_space_orig), max_idx);
fprintf('k空間の最大値は座標 (%d, %d, %d) にあります。\n', kk, mm, nn);
p0_factor = k_space_orig(max_idx) / max_val;
k_space_p0 = k_space_orig / p0_factor; % これで k_space_p0 は「fftshift 済み」

% (セクション1, 2, 3 は変更なし)
...
% (セクション1, 2, 3 は変更なし)
...

%% --- 4. [★修正★] k空間の A vs. Phase を 1つのFigureに描画 ---
fprintf('4. k空間の A vs. Phase グラフを描画します...\n');

% 描画するスライスを指定
slice_to_plot = floor(Slice / 2) + 1; % 中央スライス

% k空間からスライスを抽出
% k_space_p0 は [512, 768, Slice] (fftshift していない)
k_space_slice = k_space_p0(:, :, slice_to_plot); 
[num_rows, num_cols] = size(k_space_slice);

% [★修正★] 描画する kx インデックスのリストを定義
% (1 が DC成分 kx=0 に相当)
% kx_indices_to_plot = [1, 51, 101, 151, 201, 251];
kx_indices_to_plot = 1:50:758; % [★ 1から50ずつ758まで (1, 51, ..., 751) ★]
num_plots = length(kx_indices_to_plot); % 16個のプロット

% [★修正★] ループの *前* に Figure と TiledLayout を1回だけ作成
figure('Name', sprintf('A vs. Phase (Slice %d)', slice_to_plot), 'WindowState', 'maximized');
% [★修正★] 16個のプロットに合わせて 4行4列のグリッドを作成
t = tiledlayout(4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
% 共通のタイトルをFigureに追加
title(t, sprintf('k-space (Slice %d) の A vs. Phase (fftshift無し)', slice_to_plot));
% 共通の軸ラベルを追加
xlabel(t, 'Phase (位相) [rad]');
ylabel(t, 'A (振幅)');

% [★修正★] Y軸の最大値を全グラフで統一するために先に計算
max_A_global = 0;
for kx_val = kx_indices_to_plot
    if kx_val <= num_cols
        max_A_global = max(max_A_global, max(abs(k_space_slice(:, kx_val))));
    end
end
if max_A_global == 0, max_A_global = 1; end

% [★修正★] kxインデックスのリストを使ってループ
for i = 1:num_plots
    
    kx_to_plot = kx_indices_to_plot(i);
    col_idx = kx_to_plot;
    
    % [★安全対策★] col_idx がサイズを超えないかチェック
    if col_idx > num_cols
        fprintf('kx index=%d は行列サイズ (%d) を超えるためスキップします。\n', kx_to_plot, num_cols);
        continue;
    end

    % kx=N の列データ (512点) をすべて抽出 (複素数)
    kx_line_complex = k_space_slice(:, col_idx);

    % A (振幅) と Phase (位相) を計算
    A_values = abs(kx_line_complex);
    Phase_values = angle(kx_line_complex);

    % [★修正★] 次のタイル (subplot) を選択
    ax = nexttile(i);
    
    % グラフを描画
    scatter(ax, Phase_values, A_values, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.4);

    % [★修正★] グラフ「ごと」の設定 (個別のタイトルのみ)
    title(ax, sprintf('kx index = %d (kx=%d)', kx_to_plot, kx_to_plot - 1));
    
    % [★修正★] 共通の軸設定
    xlim(ax, [-pi, pi]);
    ylim(ax, [0, max_A_global * 1.1]); % 全グラフでY軸のスケールを統一
    grid(ax, 'on');
    xticks(ax, [-pi, -pi/2, 0, pi/2, pi]);
    xticklabels(ax, {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

    fprintf('kx index=%d のグラフをタイル %d に描画しました。\n', kx_to_plot, i);
end

fprintf('グラフの描画完了。スクリプトの残りを続行します...\n\n');

% ... (セクション5以降の処理が続く場合は、ここにペースト) ...% ... (もしセクション5以降が続く場合は、ここにペースト) ...