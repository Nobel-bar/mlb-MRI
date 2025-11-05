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
num_elements = prod(dims); 
precision = 'double=>double'; 
% --- 強度・位相データの読み込み ---
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
% .matファイルを読み込むと仮定します。
% .matファイル以外 (DICOM, P-fileなど) の場合は、このセクションの変更が必要です。

try
    fprintf('kスペースデータを読み込んでいます...\n');
    % .matファイルを構造体として読み込む
    echo_idx= 1;
    original_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * iPhase_4D(:,:,:, echo_idx));
    k_space_complex = fft2(original_img(:,:, round(params.matrix_size(3) /2)));
    
    fprintf('データの読み込みに成功しました。\n');

end

%% 3. kスペースの準備と再構成ステップ
% -------------------------------------------------
fprintf('画像再構成ステップを開始します...\n');

% kスペースのゼロ周波数成分を中央にシフト
k_space_shifted = fftshift(k_space_complex);

[rows, cols] = size(k_space_shifted);
center_r = floor(rows / 2) + 1;
center_c = floor(cols / 2) + 1;

% 25ステップでマスクを広げる
num_steps = 25;

% マスクの「半径」(正方形の一辺の半分の長さ)を定義
% 最小(ほぼゼロ)から最大(kスペース全体をカバー)まで25ステップに分割
max_radius = floor(max(rows, cols) / 2);
radius_steps = round(linspace(max_radius / num_steps, max_radius, num_steps));
radius_steps = unique(radius_steps); % 重複ステップを削除
num_steps = length(radius_steps); % 実際のステップ数に更新

% 5x5のタイルレイアウトを持つFigureを準備
figure('Name', 'kスペース中心からのデータ追加による画像再構成', 'Position', [100, 100, 1200, 900]);
t_images = tiledlayout(5, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_images, 'kスペース中心からのデータ追加による画像再構成', 'FontSize', 16, 'FontWeight', 'bold');

% (オプション) マスク自体の形状を表示するFigure
figure('Name', '使用したkスペース領域 (マスク)', 'Position', [1350, 100, 1200, 900]);
t_masks = tiledlayout(5, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t_masks, '使用したkスペース領域 (マスク)', 'FontSize', 16, 'FontWeight', 'bold');

% ステップごとに処理
for i = 1:num_steps
    current_radius = radius_steps(i);
    
    % 1. kスペースマスクの作成
    mask = zeros(rows, cols);
    r_start = max(1, center_r - current_radius);
    r_end   = min(rows, center_r + current_radius);
    c_start = max(1, center_c - current_radius);
    c_end   = min(cols, center_c + current_radius);
    
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
    title(sprintf('Step %d/%d (Mask: %dpx)', i, num_steps, current_radius*2 + 1));
    axis off;

    % 5. (オプション) マスクを表示
    figure(t_masks.Parent); % マスク表示用のFigureをアクティブに
    nexttile;
    imshow(mask, []);
    title(sprintf('Step %d/%d', i, num_steps));
    axis off;
    
    % 進行状況をコマンドウィンドウに表示
    if mod(i, 5) == 0 || i == num_steps
        fprintf('  Step %d/%d 完了 (k-space mask size: %d x %d)\n', i, num_steps, (r_end-r_start)+1, (c_end-c_start)+1);
    end
end

fprintf('シミュレーションが完了しました。\n');