% アンラップした結果を一度にすべて，前後で比較するためのプログラム．3D
% QSM_processingを使用した結果に適用．


%% 1. iFreq_raw の表示 (中央のスライスを表示する例)

% --- 1. 初期設定 ---
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H_local = 'C:\Users\hamaguchi\Downloads\matlab\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
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
% 読み込みパスと保存パスを定義
load_base_path = fullfile(image_file_0, image_file_3);

load(fullfile(load_base_path, 'phase.mat'));
load(fullfile(load_base_path, 'Mask.mat'));

fprintf('データの読み込みが完了しました。\n');

% PDFの入力(iFreq), 出力(RDF), およびマスク(Mask)を読み込む
%  load(data_file, 'iFreq', 'RDF', 'Mask');



%% 1. マスクの確認
if ~exist('Mask', 'var')
    warning('変数 "Mask" が見つかりません。');
    fprintf('3D表示には脳マスク "Mask" が必要です。\n');
    fprintf('マスクなしで表示を試みますが、背景ノイズも表示されるため推奨されません。\n');
else
    % マスクを利用して、脳領域外を 0 に設定
    % volshow は 0 の値を自動的に透明として扱います
    iFreq= iFreq .* Mask;
    iFreq_raw= iFreq_raw .* Mask;
end

fprintf('読み込んでいます...');

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
        sgtitle(['3D Mesh Wrapping Comparison - Figure ' num2str(figure_num) ...
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
    
    % --- 左側: Wrap適用前 (iFreq_raw) ---
    ax_before = subplot(rows, cols, plot_idx_before);
    
    % ★変更点: データを準備 (マスク外をNaNに)
    data_before = iFreq_raw(:, :, slice_idx);
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
        title('Before Wrap (iFreq_raw)');
    end

    % --- 右側: Wrapping適用後 (iFreq) ---
    ax_after = subplot(rows, cols, plot_idx_after);
    
    % ★変更点: データを準備 (マスク外をNaNに)
    data_after = iFreq(:, :, slice_idx);
    data_after(current_mask == 0) = NaN;
    
    % ★変更点: meshでプロット
    mesh(ax_after, data_after);

    axis tight;
    axis on;
    daspect([50 50 1]); % ★変更点: Z軸を強調 (50 50 1)
    colormap(ax_after, 'default');
    colorbar(ax_after);
    
    if row_idx == 1
        title('After Wrap (iFreq)');
    end
    
end

fprintf('全スライスの 3D (mesh) 比較表示が完了しました。\n');