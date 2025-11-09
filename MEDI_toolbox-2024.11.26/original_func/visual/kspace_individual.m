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

image_file_0 = image_file_00; % slab用
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
% 5x5 グリッドを定義
grid_size = 5;
radius = floor(grid_size / 2); % 中心から 2 ピクセル ( -2, -1, 0, 1, 2 )
num_steps = grid_size * grid_size;
fprintf('  カスタムステップを生成しました。\n');
fprintf('  kスペース中心 %dx%d (%d ピクセル) を個別再構成します。\n', grid_size, grid_size, num_steps);
% -------------------------------------------------
% ステップ数に応じてタイルレイアウトを動的に決定
n_rows = grid_size;
n_cols = grid_size;
% タイルレイアウトを持つFigureを準備
figure('Name', 'kスペース中心 5x5ピクセルの個別再構成 (振幅/Logスケール)', 'Position', [100, 100, 1000, 1000]);
t_images = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_images, sprintf('kスペース中心 %dx%dピクセルの個別再構成画像 (振幅 / Logスケール)', grid_size, grid_size), 'FontSize', 16, 'FontWeight', 'bold');
% (オプション) マスク自体の形状を表示するFigureは削除
% ステップごとに処理
tile_idx = 1;
for r_offset = -radius:radius  % -2, -1, 0, 1, 2
    for c_offset = -radius:radius  % -2, -1, 0, 1, 2
        
        % 現在のkスペースピクセル座標
        current_r = center_r + r_offset;
        current_c = center_c + c_offset;
        
        % 1. kスペースマスクの作成 (1ピクセルのみ)
        mask = zeros(rows, cols);
        mask(current_r, current_c) = 1;
        
        % 2. kスペースにマスクを適用
        %    (特定の1ピクセルの値だけを抽出)
        k_space_masked = k_space_shifted .* mask;
        
        % 3. 画像空間へ逆フーリエ変換
        % (ifftshiftでゼロ周波数をコーナーに戻してからifft2)
        image_recon = ifft2(ifftshift(k_space_masked));
        
        % 4. 再構成画像（振幅画像）を表示
        figure(t_images.Parent); % 画像表示用のFigureをアクティブに
        nexttile(tile_idx);
        
        % 
        % [変更]
        % 振幅(abs)を対数スケール(log)で表示する
        % log(1 + ...) で 0 の対数を防ぐ
        % [] で全タイルの log(1+abs) の最小/最大値でスケーリングする
        imshow(log(1 + imag(image_recon)), []);
        
        title_str = sprintf('k-Pixel (%d, %d)', current_r, current_c);
        title(title_str, 'FontSize', 8);
        axis off;
        % 5. (オプション) マスクを表示 ... 削除
        
        % 進行状況をコマンドウィンドウに表示
        fprintf('  %s 完了\n', title_str);
        
        tile_idx = tile_idx + 1;
    end
end
fprintf('シミュレーションが完了しました。\n');
