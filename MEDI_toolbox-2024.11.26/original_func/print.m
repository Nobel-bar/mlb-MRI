%==================================================================================================
% [★新規★] 保存された「処理済みMRI画像」(.raw) を読み込み、
% 3Dスライス断面で表示するスクリプト
%
% [★注意★]
% 直前のスクリプト [cite: 該当コード] と同じファイルパス、ファイル名、画像サイズ設定を使用します。
%==================================================================================================

fprintf('スクリプトを開始します (処理済みMRI画像 3D表示)\n');
clear variables;
close all;

%% --- 1. 初期設定 (直前のスクリプト [cite: 該当コード] と設定を合わせる) ---
fprintf('1. パラメータを設定しています...\n');

% --- 直前のスクリプト [cite: 該当コード] から設定をコピー ---
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!% 2Drotate 
image_file_000 = "C:\Users\hamaguchi\Downloads\matlab";

% ★ 読み込むデータに合わせて、この2行をアクティブ化してください
% image_file_0 = image_file_00; % 0deg データの場合
image_file_0 = image_file_2DGE_1_2_Rotate_H; % 1-2_Rotate データの場合
image_file_0 = image_file_000; % local用

image_file_1 = '1_data';
image_file_2 = '2_original_data';

% ★ 読み込むデータに合わせて、filename_base を選択してください
% filename_base = sprintf('1st_2DGE_0deg_15'); 
filename_base = sprintf('1st_2DGE_1_2_Rotate'); 
filename_base = sprintf('spin');

% filename_base = sprintf('3d_3D_rotate_direct_th-18.6'); 


% 読み込みパスを定義 (直前のスクリプト [cite: 該当コード] の `save_path` [cite: 該当コード])
load_path = fullfile(image_file_0, image_file_1); % '1_data' フォルダから読み込む
% load_path = "F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H\4_rolate_output_data\direct";

% 読み込むファイル名
load_mag_name = [filename_base, '_mag.raw']; % 強度画像を読み込む

% 読み込む画像のサイズ (直前のスクリプト [cite: 該当コード] の `final_img` [cite: 該当コード] のサイズ)
% params.matrix_size = [512, 512, 23];
final_matrix_x = 512;
final_matrix_y = 512;
% --- ここまで設定 ---


%% --- 2. 処理済み画像データ (.raw) の読み込み ---
fprintf('2. 処理済み画像 (%s) を読み込んでいます...\n', load_mag_name);
filename_input_Mag = fullfile(load_path, load_mag_name);

fileID_Mag = fopen(filename_input_Mag, 'r');
if fileID_Mag == -1, error('ファイルが開けませんでした: %s', filename_input_Mag); end

% final_img [cite: 該当コード] は 'double' で保存された [cite: 該当コード]
data_vector_mag = fread(fileID_Mag, inf, 'double'); 
fclose(fileID_Mag);
Slice = numel(data_vector_mag) / (final_matrix_x * final_matrix_y);
if mod(Slice, 1) ~= 0, error('ファイルサイズが不正です。'); end


% 3D配列に変換
% V のサイズは [512 (X), 512 (Y), 23 (Z)]
V = reshape(data_vector_mag, [final_matrix_x, final_matrix_y, Slice]);
fprintf('%d x %d x %d の処理済み画像を正常に読み込みました。\n', final_matrix_x, final_matrix_y, Slice);


%% --- 3. [★修正★] 2Dスライスを 4x4 グリッドでプロット (複数Figure) ---
fprintf('3. 全%dスライスを 4x4 グリッドで一覧表示します...\n', Slice);

num_rows = 4;
num_cols = 4;
tiles_per_figure = num_rows * num_cols; % 1Figureあたり16スライス

% 必要なFigureの数を計算
num_figures = ceil(Slice / tiles_per_figure);

for fig_idx = 1:num_figures
    % 新しいFigureを作成
    figure('Name', sprintf('Processed MRI Slices (Figure %d/%d)', fig_idx, num_figures), 'WindowState', 'maximized');
    t = tiledlayout(num_rows, num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % このFigureで表示するスライスの開始番号と終了番号を計算
    slice_start = (fig_idx - 1) * tiles_per_figure + 1;
    slice_end = min(fig_idx * tiles_per_figure, Slice); % 最後のスライスを超えないように
    
    title(t, sprintf('処理済みMRI画像 スライス %d-%d (強度)', slice_start, slice_end));

    % このFigureにスライスを表示
    tile_idx = 1;
    for s = slice_start:slice_end
        % 次のタイル (subplot) を選択
        ax = nexttile(tile_idx);
        
        % スライス (:, :, s) を抽出して表示
        % [★構文エラー修正★] imshow(..., 'Parent', ax) の形式を使用
        imshow(V(:, :, s), [], 'Parent', ax);
        
        % 軸を消し、タイトルにスライス番号を付ける
        axis(ax, 'off');
        title(ax, sprintf('Slice %d / %d', s, Slice), 'FontSize', 8);
        
        tile_idx = tile_idx + 1;
    end
    
    % --- 空のタイルを非表示にする (最後のFigure用) ---
    for t_idx = tile_idx : tiles_per_figure
        ax = nexttile(t_idx);
        axis(ax, 'off');
    end

    % カラーマップをグレースケールに設定
    colormap(gray); 
end

fprintf('完了。%d 枚のFigureに必要なスライスをすべて表示しました。\n', num_figures);
fprintf('完了。\n');% カラーマップをグレースケールに設定
colormap(gray); 

fprintf('完了。\n');