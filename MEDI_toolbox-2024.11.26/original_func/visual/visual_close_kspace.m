clear;
close all;
clc;
fprintf('kスペース再構成シミュレーションを開始します...\n');
%% 1. ユーザー定義のパス設定
% -------------------------------------------------
% ユーザーが提供したパス情報をそのまま使用
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_000 = "C:\Users\hamaguchi\Downloads\matlab";
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 
image_file_0 = image_file_000; % local用

% 読み込みパスを定義
load_base_path = fullfile(image_file_0, image_file_1);

% 入力ファイル名 (拡張子なし)
%% --- 2. rawファイルのパスを指定 ---
% ご自身のファイルが保存されている実際のパスに書き換えてください。
mag_filename = '1st_2DGE_0deg_mag.raw'; % !! 要変更 !!
phase_filename = '1st_2DGE_0deg_phase.raw'; % !! 要変更 !!

mag_filepath = fullfile(load_base_path, mag_filename);
phase_filepath = fullfile(load_base_path, phase_filename);

% --- params構造体: QSMデータの撮像パラメータ ---
params = struct();
params.matrix_size = [512, 512, 23]; % !! 要変更 !! : 行列サイズ
params.TE = 0.015; % !! 要変更 !! : エコー時間 (秒)
dims = [params.matrix_size, length(params.TE)];
precision = 'double=>double'; 

% --- 強度・位相データの読み込み ---
fid_mag = fopen(mag_filepath, 'rb');
if fid_mag == -1, error('強度ファイルを開けませんでした: %s', mag_filepath); end
iMag_4D = reshape(fread(fid_mag, inf, precision), dims);
fclose(fid_mag);
fid_phase = fopen(phase_filepath, 'rb');
if fid_phase == -1, error('位相ファイルを開けませんでした: %s', phase_filepath); end
iPhase_4D = reshape(fread(fid_phase, inf, precision), dims);
fclose(fid_phase);

%% 2. kスペースデータの読み込み
% -------------------------------------------------
try
    fprintf('kスペースデータを読み込んでいます...\n');
    % .matファイルを構造体として読み込む
    echo_idx = 1;
    % 3Dデータの中央スライスを選択
    slice_idx = round(params.matrix_size(3) / 2);
    original_img = iMag_4D(:,:,slice_idx, echo_idx) .* exp(1i * iPhase_4D(:,:,slice_idx, echo_idx));
    
    % 画像空間からkスペースへ変換
    k_space_complex = fft2(original_img);
    
    fprintf('データの読み込みに成功しました。\n');
    fprintf('  対象スライス: %d / %d\n', slice_idx, params.matrix_size(3));

catch ME
    fprintf('\n!! エラー: データの読み込みまたは変換に失敗しました !!\n');
    fprintf('エラーメッセージ: %s\n', ME.message);
    fprintf('パスやパラメータ（特に matrix_size）を確認してください。\n');
    return;
end

%% 3. kスペースの準備と再構成ステップ
% -------------------------------------------------
fprintf('画像再構成ステップを開始します...\n');

% kスペースのゼロ周波数成分を中央にシフト
k_space_shifted = fftshift(k_space_complex);
[rows, cols] = size(k_space_shifted);
center_r = floor(rows / 2) + 1;
center_c = floor(cols / 2) + 1;

% --- マスクステップの定義 (ユーザーの要求に応じて変更) ---
% ご指定のピクセル数（面積）
pixel_areas = [1, 4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144];
% マスクの一辺の長さに変換
width_steps = round(sqrt(pixel_areas));

% 最終ステップ（kスペース全体）を追加する処理を削除
% max_width = max(rows, cols);
% width_steps = [width_steps, max_width];
width_steps = unique(width_steps); % 重複を削除

num_steps = length(width_steps); % 実際のステップ数に更新
fprintf('  カスタムステップを生成しました。合計 %d ステップ。\n', num_steps);
fprintf('  使用するマスク幅 (ピクセル): %s\n', mat2str(width_steps));
% -------------------------------------------------

% ステップ数に応じてタイルレイアウトを動的に決定
n_rows = ceil(sqrt(num_steps));
n_cols = ceil(num_steps / n_rows);

% タイルレイアウトを持つFigureを準備
figure('Name', 'kスペース中心からのデータ追加による画像再構成', 'Position', [100, 100, 1200, 900]);
t_images = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_images, 'kスペース中心からのデータ追加による画像再構成', 'FontSize', 16, 'FontWeight', 'bold');

% (オプション) マスク自体の形状を表示するFigure
figure('Name', '使用したkスペース領域 (マスク)', 'Position', [1350, 100, 1200, 900]);
t_masks = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_masks, '使用したkスペース領域 (マスク)', 'FontSize', 16, 'FontWeight', 'bold');

% ステップごとに処理
for i = 1:num_steps
    current_width = width_steps(i);
    
    % 1. kスペースマスクの作成
    mask = zeros(rows, cols);
    
    % マスクの幅に基づいて開始/終了インデックスを計算
    % (奇数幅でも偶数幅でも中心に配置されるように調整)
    r_start = max(1, center_r - floor((current_width - 1) / 2));
    r_end   = min(rows, center_r + ceil((current_width - 1) / 2));
    c_start = max(1, center_c - floor((current_width - 1) / 2));
    c_end   = min(cols, center_c + ceil((current_width - 1) / 2));
    
    mask(r_start:r_end, c_start:c_end) = 1;
    
    % 2. kスペースにマスクを適用
    k_space_masked = k_space_shifted .* mask;
    
    % 3. 画像空間へ逆フーリエ変換
    % (ifftshiftでゼロ周波数をコーナーに戻してからifft2)
    image_recon = ifft2(ifftshift(k_space_masked));
    
    % 4. 再構成画像（振幅画像）を表示
    figure(t_images.Parent); % 画像表示用のFigureをアクティブに
    nexttile;
    imshow(abs(image_recon), []);
    title_str = sprintf('Step %d/%d (Mask: %dx%dpx)', i, num_steps, (c_end-c_start)+1, (r_end-r_start)+1);
    title(title_str);
    axis off;

    % 5. (オプション) マスクを表示
    figure(t_masks.Parent); % マスク表示用のFigureをアクティブに
    nexttile;
    imshow(mask, []);
    title(sprintf('Step %d/%d', i, num_steps));
    axis off;
    
    % 進行状況をコマンドウィンドウに表示
    fprintf('  %s 完了\n', title_str);
end

fprintf('シミュレーションが完了しました。\n');