% マスクの結果を，元画像と一覧で比較するためのプログラム．2D
% QSM_processingを使用した結果に適用．

%% --- 1. データの読み込み ---
% --- 1. 初期設定 ---
% --- 1. 初期設定 ---
image_file_1 = 'F:/hamaguchi/copy/20241205_RawData_H/Volunteer_Rotate_H/2DGE_0deg_H/total_slice';
image_file_2 = 'output_data';

save_path = fullfile(image_file_1, image_file_2);

load(fullfile(save_path, 'other.mat'));

fprintf('読み込んでいます...');

fprintf('データの読み込みが完了しました。\n');

fprintf('振幅画像 (iMag) と脳マスク (Mask_CSF) の比較表示を開始します...\n');

%% --- 全スライス比較プロット (4セット/Figure) ---

% 1. 基本設定
total_slices = size(iMag, 3);     % Z軸の総スライス数を取得 (iMagから)
slices_per_figure = 4;            % 1つのFigureに表示する比較セット数
rows = slices_per_figure;         % Figureの行数 (4行)
cols = 2;                         % Figureの列数 (2列: iMag/Mask_CSF)

% 画面サイズを取得してFigureの位置を決定
screen_size = get(groot, 'ScreenSize');
figure_width = 800;  % Figureの幅 (ピクセル単位)
figure_height = 900; % Figureの高さ (ピクセル単位)
% 画面中央に配置
pos_x = (screen_size(3) - figure_width) / 2;
pos_y = (screen_size(4) - figure_height) / 2;

% 2. 全スライスをループ処理
for slice_idx = 1:total_slices
    
    % 3. Figureを新規作成するタイミングか判定
    %    (slice_idxが 1, 5, 9, ... の時に新しいFigureを開く)
    if mod(slice_idx - 1, slices_per_figure) == 0
        figure; % 新しいFigureウィンドウを作成
        
        % Figureのサイズと位置を設定
        set(gcf, 'Position', [pos_x, pos_y, figure_width, figure_height]);
        
        figure_num = floor((slice_idx - 1) / slices_per_figure) + 1;
        
        % Figure全体にタイトルを追加 (MATLAB R2018b以降)
        sgtitle(['Magnitude (iMag) vs. Brain Mask_CSF (Mask_CSF) - Figure ' num2str(figure_num) ...
                 ' (Slices ' num2str(slice_idx) ' - ' ...
                 num2str(min(slice_idx + slices_per_figure - 1, total_slices)) ')']);
    end
    
    % 4. 現在のFigure内での行インデックス (1から4) を計算
    row_idx = mod(slice_idx - 1, slices_per_figure) + 1;
    
    % 5. Subplotの位置を計算
    plot_idx_left = (row_idx - 1) * cols + 1; % 左側の列
    plot_idx_right  = (row_idx - 1) * cols + 2; % 右側の列
    
    % --- 左側: 振幅画像 (iMag) ---
    subplot(rows, cols, plot_idx_left);
    imagesc(iMag(:, :, slice_idx));
    colormap(gca, 'gray'); % グレースケールに設定
    axis equal tight; % 軸をデータに合わせる
    axis off;         % 軸の目盛りを非表示
    
    % 各行の左側にスライス番号を表示
    ylabel(['Slice ' num2str(slice_idx)], 'Visible', 'on', 'FontWeight', 'bold', 'Color', 'k');
    
    % 各Figureの1行目にだけ列タイトルを表示
    if row_idx == 1
        title('Magnitude (iMag)');
    end

    % --- 右側: 脳マスク (Mask_CSF) ---
    subplot(rows, cols, plot_idx_right);
    imagesc(Mask_CSF(:, :, slice_idx));
    colormap(gca, 'gray'); % マスクもグレースケール (0=黒, 1=白)
    axis equal tight;
    axis off;
    
    if row_idx == 1
        title('Brain Mask_CSF (Mask_CSF)');
    end
    
end

fprintf('全 %d スライスの比較プロットが完了しました。\n', total_slices);