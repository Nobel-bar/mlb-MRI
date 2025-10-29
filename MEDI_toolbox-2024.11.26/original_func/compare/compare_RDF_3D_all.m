%================================================================
% PDF (背景磁場除去) 結果 全スライス3D(Mesh)比較プログラム
% 
% 概要:
%   iFreq (PDF適用前) と RDF (PDF適用後) を、全スライスにわたって
%   3Dメッシュプロットで並べて表示します。
%
%   - 1つのFigureに4スライス分の比較 (計8プロット) を表示します。
%   - 脳マスク(Mask)を NaN に変換して適用します。
%   - Z軸(高さ)を daspect で強調表示します。
%
% 依存ファイル:
%   - 'RDF.mat' (QSM_processing.m で生成。iFreq, RDF, Mask を含む)
%================================================================

%% --- 1. データの読み込み ---
% --- 1. 初期設定 ---
% --- 1. 初期設定 ---
image_file_1 = 'F:/hamaguchi/copy/20241205_RawData_H/Volunteer_Rotate_H/2DGE_0deg_H/total_slice';
image_file_2 = 'output_data';

save_path = fullfile(image_file_1, image_file_2);

load(fullfile(save_path, 'phase.mat'));
load(fullfile(save_path, 'PDF.mat'));
load(fullfile(save_path, 'Mask.mat'));

fprintf('読み込んでいます...');
% PDFの入力(iFreq), 出力(RDF), およびマスク(Mask)を読み込む
%  load(data_file, 'iFreq', 'RDF', 'Mask');

fprintf('データの読み込みが完了しました。\n');
% 必要な変数がロードされたかチェック
if ~exist('iFreq') || ~exist('RDF')
    error('RDF.mat ファイルに iFreq または RDF 変数が含まれていません。');
end
if ~exist('Mask')
    warning('変数 "Mask" が見つかりません。マスクなしで続行します。');
    % マスクがない場合は、すべて1のダミーマスクを作成
    Mask = ones(size(iFreq, 1), size(iFreq, 2), size(iFreq, 3)); 
end

fprintf('全スライスの比較画像の生成を開始します...\n');


%% --- 2. 全スライス比較プロット (2セット/Figure) ---

% 1. 基本設定
total_slices = size(iFreq, 3); % Z軸の総スライス数を取得
slices_per_figure = 2;           % 1つのFigureに表示する比較セット数
rows = slices_per_figure;        % Figureの行数 (4行)
cols = 2;                        % Figureの列数 (2列: Before/After)

% 画面サイズを取得してFigureの位置を決定
screen_size = get(groot, 'ScreenSize');
figure_width = 800;  % Figureの幅 (ピクセル単位)
figure_height = 900; % Figureの高さ (ピクセル単位)
pos_x = (screen_size(3) - figure_width) / 2; % 画面中央X
pos_y = (screen_size(4) - figure_height) / 2; % 画面中央Y

% 2. 全スライスをループ処理
for slice_idx = 1:total_slices
    
    % 3. Figureを新規作成するタイミングか判定
    if mod(slice_idx - 1, slices_per_figure) == 0
        fig = figure; 
        set(fig, 'Position', [pos_x, pos_y, figure_width, figure_height]);
        set(fig, 'Color', 'w'); % Figureの背景色を白に
        
        figure_num = floor((slice_idx - 1) / slices_per_figure) + 1;
        
        % Figure全体にタイトルを追加
        sgtitle(['3D Mesh PDF Comparison - Figure ' num2str(figure_num) ...
                 ' (Slices ' num2str(slice_idx) ' - ' ...
                 num2str(min(slice_idx + slices_per_figure - 1, total_slices)) ')']);
    end
    
    % 4. 現在のFigure内での行インデックス (1から4) を計算
    row_idx = mod(slice_idx - 1, slices_per_figure) + 1;
    
    % 5. Subplotの位置を計算
    plot_idx_before = (row_idx - 1) * cols + 1; % 左側の列
    plot_idx_after  = (row_idx - 1) * cols + 2; % 右側の列
    
    % 現在のスライスのマスクを取得
    current_mask = Mask(:, :, slice_idx);
    
    % --- 左側: PDF適用前 (iFreq) ---
    ax_before = subplot(rows, cols, plot_idx_before);
    
    % ★変更点: データを準備 (マスク外をNaNに)
    data_before = iFreq(:, :, slice_idx);
    data_before(current_mask == 0) = NaN;
    
    % ★変更点: meshでプロット
    mesh(ax_before, data_before);
    
    axis tight;
    axis on; % 3Dでは軸をONに
    daspect([50 50 1]); % ★変更点: Z軸を強調 (50 50 1)
    colormap(ax_before, 'default');
    colorbar(ax_before);
    
    % 各行の左側にスライス番号を表示 (2Dの時と同じ)
    ylabel(['Slice ' num2str(slice_idx)], 'Visible', 'on', 'FontWeight', 'bold', 'Color', 'k');
    
    % 各Figureの1行目にだけ列タイトルを表示
    if row_idx == 1
        title('Before PDF (iFreq)');
    end

    % --- 右側: PDF適用後 (RDF) ---
    ax_after = subplot(rows, cols, plot_idx_after);
    
    % ★変更点: データを準備 (マスク外をNaNに)
    data_after = RDF(:, :, slice_idx);
    data_after(current_mask == 0) = NaN;
    
    % ★変更点: meshでプロット
    mesh(ax_after, data_after);

    axis tight;
    axis on;
    daspect([50 50 1]); % ★変更点: Z軸を強調 (50 50 1)
    colormap(ax_after, 'default');
    colorbar(ax_after);
    
    if row_idx == 1
        title('After PDF (RDF)');
    end
    
end

fprintf('全スライスの 3D (mesh) 比較表示が完了しました。\n');