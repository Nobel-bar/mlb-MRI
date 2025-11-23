% アンラップした結果を一度にすべて，前後で比較するためのプログラム．2D
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

fprintf('データの読み込みが完了しました。\n');


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

%% --- 全スライス比較プロット (4セット/Figure) ---

% 1. 基本設定
total_slices = size(iFreq_raw, 3); % Z軸の総スライス数を取得
slices_per_figure = 4;            % 1つのFigureに表示する比較セット数
rows = slices_per_figure;         % Figureの行数 (4行)
cols = 2;                         % Figureの列数 (2列: Before/After)

% 画面サイズを取得してFigureの位置を決定 (★追加)
screen_size = get(groot, 'ScreenSize');
figure_width = 800;  % Figureの幅 (ピクセル単位, お好みで調整してください)
figure_height = 900; % Figureの高さ (ピクセル単位, お好みで調整してください)
% 画面中央に配置
pos_x = (screen_size(3) - figure_width) / 2;
pos_y = (screen_size(4) - figure_height) / 2;

% 2. 全スライスをループ処理
for slice_idx = 1:total_slices
    
    % 3. Figureを新規作成するタイミングか判定
    %    (slice_idxが 1, 5, 9, ... の時に新しいFigureを開く)
    if mod(slice_idx - 1, slices_per_figure) == 0
        figure; % 新しいFigureウィンドウを作成
        
        % ★Figureのサイズと位置を設定 (★追加)
        set(gcf, 'Position', [pos_x, pos_y, figure_width, figure_height]);
        
        figure_num = floor((slice_idx - 1) / slices_per_figure) + 1;
        
        % Figure全体にタイトルを追加 (MATLAB R2018b以降)
        % もし 'sgtitle' がエラーになる場合は、この行をコメントアウトしてください
        sgtitle(['Unwrapping Comparison - Figure ' num2str(figure_num) ...
                 ' (Slices ' num2str(slice_idx) ' - ' ...
                 num2str(min(slice_idx + slices_per_figure - 1, total_slices)) ')']);
    end
    
    % 4. 現在のFigure内での行インデックス (1から4) を計算
    row_idx = mod(slice_idx - 1, slices_per_figure) + 1;
    
    % 5. Subplotの位置を計算
    plot_idx_before = (row_idx - 1) * cols + 1; % 左側の列
    plot_idx_after  = (row_idx - 1) * cols + 2; % 右側の列
    
    % --- 左側: アンラップ前 (iFreq_raw) ---
    subplot(rows, cols, plot_idx_before);
    imagesc(iFreq_raw(:, :, slice_idx));
    axis equal tight; % 軸をデータに合わせる
    axis off;         % 軸の目盛りを非表示
    
    % 各行の左側にスライス番号を表示
    ylabel(['Slice ' num2str(slice_idx)], 'Visible', 'on', 'FontWeight', 'bold', 'Color', 'k');
    
    % 各Figureの1行目にだけ列タイトルを表示
    if row_idx == 1
        title('Before (iFreq\_raw)');
    end

    % --- 右側: アンラップ後 (iFreq) ---
    subplot(rows, cols, plot_idx_after);
    imagesc(iFreq(:, :, slice_idx));
    axis equal tight;
    axis off;
    
    if row_idx == 1
        title('After (iFreq)');
    end
    
end

fprintf('アンラッピング前後の画像を比較表示します...\n');

%{
%% アンラップ前後の比較表示

% 3次元データの真ん中のスライス番号を取得
% size(iFreq_raw, 3) は、z軸の枚数を取得します
slice = round(size(iFreq_raw, 3));

slice_to_display = 5;
if slice_to_display < 1 || slice_to_display > slice
    fprintf('(%d)番目のスライス を表示しました。\n', slice_to_display)
end

% 新しい図（Figure）ウィンドウを開く
figure; 

% 3. 1行2列の表示レイアウトを作成 (MATLAB R2019b以降)
tiledlayout(1, 2, 'TileSpacing', 'compact'); 

% --- 左側のプロット (アンラップ前) ---
nexttile; % 1番目のタイルを選択
imagesc(iFreq_raw(:, :, slice_to_display));
colorbar;
axis equal tight; % 軸のスケールをデータに合わせる
title(['Before Unwrapping (iFreq\_raw)' newline 'Slice ' num2str(slice_to_display)]);
xlabel('X'); ylabel('Y');

% --- 右側のプロット (アンラップ後) ---
nexttile; % 2番目のタイルを選択
imagesc(iFreq(:, :, slice_to_display));
colorbar;
axis equal tight;
title(['After Unwrapping (iFreq)' newline 'Slice ' num2str(slice_to_display)]);
xlabel('X');

fprintf('比較ウィンドウを表示しました。\n');
%}

%{
%% 2. iFreq_raw の保存 (MATLABの .mat 形式)

% save(ファイル名, 保存したい変数名)
save('iFreq_raw_output.mat', 'iFreq_raw');

fprintf('iFreq_raw 変数を ''iFreq_raw_output.mat'' として保存しました。\n');

% 読み込むときは load('iFreq_raw_output.mat') を使います。

% --- (オプション) NIfTI形式での保存 ---
% FSLやSPMなど、他の画像解析ソフトで使いたい場合は .nii 形式が便利です。
% この場合、Image Processing Toolbox が必要です。
% また、正しい表示のために voxel_size が必要です。
%{
if exist('niftiwrite', 'file') && exist('voxel_size', 'var')
    niftiwrite(iFreq_raw, 'iFreq_raw.nii', 'VoxelSize', voxel_size);
    fprintf('iFreq_raw 変数を ''iFreq_raw.nii'' としてNIfTI形式でも保存しました。\n');
else
    fprintf('NIfTI形式での保存はスキップされました。(niftiwrite関数またはvoxel_size変数がありません)\n');
end
%}
%}