%==================================================================================================
% 実際のMRI RAWデータ (3D) を読み込み、P0補正を行い、
% k空間の A vs. Phase グラフ（Y軸対数）を描画
% 2D表示
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
k_space_orig = fftn(orig_img);
[max_val, max_idx] = max(abs(k_space_orig(:)));
[kk, mm, nn] = ind2sub(size(k_space_orig), max_idx);
fprintf('k空間の最大値は座標 (%d, %d, %d) にあります。\n', kk, mm, nn);
p0_factor = k_space_orig(max_idx) / max_val;
k_space_p0 = k_space_orig / p0_factor; % これで k_space_p0 は「fftshift 済み」

%% --- 4. [★修正★] k空間の A vs. Phase を 1つのFigureに描画 (Y軸対数) ---
fprintf('4. k空間の A vs. Phase グラフ（Y軸対数）を描画します...\n');

% 描画するスライスを指定
slice_to_plot = floor(Slice / 2) + 1; % 中央スライス

% k空間からスライスを抽出
k_space_slice = k_space_p0(:, :, slice_to_plot); 
[num_rows, num_cols] = size(k_space_slice);

% 描画する kx インデックスのリストを定義
kx_indices_to_plot = 1:20:orig_matrix_y; % [★ 1から50ずつ758まで (1, 51, ..., 751) ★]
num_plots = length(kx_indices_to_plot); % 16個のプロット

% [★修正★] ループの *前* に Figure と TiledLayout を1回だけ作成
figure('Name', sprintf('A vs. Phase (Slice %d) - Y-Log Scale', slice_to_plot), 'WindowState', 'maximized');
% [★修正★] 16個のプロットに合わせて 4行4列のグリッドを作成
t = tiledlayout(6, round(num_plots/6)+1, 'TileSpacing', 'compact', 'Padding', 'compact');
% 共通のタイトルをFigureに追加
title(t, sprintf('k-space (Slice %d) の A vs. Phase (fftshift無し) - Y軸対数スケール', slice_to_plot));
% 共通の軸ラベルを追加
xlabel(t, 'Phase (位相) [rad]');
ylabel(t, 'A (振幅) - log scale'); % [★修正★] Y軸ラベル変更

% [★修正★] Y軸の最大値・最小値を全グラフで統一するために先に計算
max_A_global = 0;
min_A_global = Inf;
for kx_val = kx_indices_to_plot
    if kx_val <= num_cols
        line_data = abs(k_space_slice(:, kx_val));
        max_A_global = max(max_A_global, max(line_data));
        % 0を除外した最小値を計算
        min_val = min(line_data(line_data > 0)); 
        if ~isempty(min_val)
            min_A_global = min(min_A_global, min_val);
        end
    end
end
if max_A_global == 0, max_A_global = 1; end
if isinf(min_A_global), min_A_global = eps; end % 全て0だった場合のフォールバック

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
    % [★修正★] 0は対数プロットできないため、epsを加算
    A_values = abs(kx_line_complex) + eps; 
    Phase_values = angle(kx_line_complex);

    % [★修正★] 次のタイル (subplot) を選択
    ax = nexttile(i);
    
    % グラフを描画
    scatter(ax, Phase_values, A_values, 20, 'b', 'filled', 'MarkerFaceAlpha', 0.4);
    
    % [★追加★] Y軸を対数スケールに設定
    set(ax, 'YScale', 'log');

    % [★修正★] グラフ「ごと」の設定 (個別のタイトルのみ)
    title(ax, sprintf('kx index = %d (kx=%d)', kx_to_plot, kx_to_plot - 1));
    
    % [★修正★] 共通の軸設定
    xlim(ax, [-pi, pi]);
    % [★修正★] Y軸のスケールを、計算した最小・最大値で統一
    ylim(ax, [min_A_global * 0.9, max_A_global * 1.1]); 
    grid(ax, 'on');
    xticks(ax, [-pi, -pi/2, 0, pi/2, pi]);
    xticklabels(ax, {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

    fprintf('kx index=%d のグラフをタイル %d に描画しました。\n', kx_to_plot, i);
end

fprintf('グラフの描画完了。スクリプトの残りを続行します...\n\n');
% ... (セクション5以降の処理が続く場合は、ここにペースト) ...