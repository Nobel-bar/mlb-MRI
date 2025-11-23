%% 3Dアンラップ結果の比較プログラム（体動なし vs. 体動あり 2種）
% QSM_processingを使用した結果（RDF）を比較します。


% !! 以下 3 行のファイルパスは、ユーザーの環境に合わせて適切に設定してください !!
image_file_2DGE_0deg_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_2DGE_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_Rotate_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H_local = 'C:\Users\hamaguchi\Downloads\matlab\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; % スペース修正
image_file_4 = '4_rolate_output_data'; % スペース修正
image_file_5 = '5_fitting_output_data'; % スペース修正

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 体動なし (No Motion) データセットのベースパス
base_path_no_motion = image_file_2DGE_0deg_H;
% 体動あり (Motion) データセットのベースパス
base_path_with_artifact = image_file_2DGE_1_2_Rotate_H; % スペース修正
base_path_with_motion = image_file_2DGE_Rotate_H; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更点%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 読み込みパスの定義
load_path_no_motion = fullfile(base_path_no_motion, output_folder);
load_path_with_artifact = fullfile(base_path_with_artifact, output_folder);
load_path_with_motion = fullfile(base_path_with_motion, output_folder);


% --- 2. データの読み込み ---

% 1. 体動なし (No Motion) データの読み込み
fprintf('1/3: 体動なし (No Motion) データを読み込み中: %s\n', load_path_no_motion);
try
    load(fullfile(load_path_no_motion, 'PDF.mat'), 'RDF'); 
    load(fullfile(load_path_no_motion, 'Mask.mat'), 'Mask');
    RDF_no_motion = RDF;
    Mask_no_motion = Mask;
    clear RDF Mask;
    fprintf('読み込み完了。\n');
catch ME
    fprintf(2, 'エラー: 体動なしデータ (No Motion) の読み込みに失敗しました。\n');
    disp(ME.message);
    return;
end

% 2. 体動あり 1 (With Artifact) データの読み込み
fprintf('2/3: 体動あり 1 (With Artifact) データを読み込み中: %s\n', load_path_with_artifact);
try
    load(fullfile(load_path_with_artifact, 'PDF.mat'), 'RDF');
    load(fullfile(load_path_with_artifact, 'Mask.mat'), 'Mask');
    RDF_with_artifact = RDF;
    Mask_with_artifact = Mask;
    clear RDF Mask;
    fprintf('読み込み完了。\n');
catch ME
    fprintf(2, 'エラー: 体動あり 1 データ (With Artifact) の読み込みに失敗しました。\n');
    disp(ME.message);
    return;
end

% 3. 体動あり 2 (With Motion) データの読み込み
fprintf('3/3: 体動あり 2 (With Motion) データを読み込み中: %s\n', load_path_with_motion);
try
    load(fullfile(load_path_with_motion, 'PDF.mat'), 'RDF');
    load(fullfile(load_path_with_motion, 'Mask.mat'), 'Mask');
    RDF_with_motion = RDF;
    Mask_with_motion = Mask;
    clear RDF Mask;
    fprintf('読み込み完了。\n');
catch ME
    fprintf(2, 'エラー: 体動あり 2 データ (With Motion) の読み込みに失敗しました。\n');
    disp(ME.message);
    return;
end

% --- 3. データのチェックと初期設定 ---

% データサイズの確認
if ~(isequal(size(RDF_no_motion), size(RDF_with_artifact)) && isequal(size(RDF_no_motion), size(RDF_with_motion)))
    warning('読み込んだ3つのデータセットのサイズが異なります。続行できません。');
    return;
end

total_slices = size(RDF_no_motion, 3); % Z軸の総スライス数を取得
slices_per_figure = 1;          % 1つのFigureに表示するスライス数 (1スライス分)
rows = slices_per_figure;       % Figureの行数 (1行)
cols = 3;                       % Figureの列数 (3列: No Motion / Artifact / Motion)

% 画面サイズを取得してFigureの位置を決定
screen_size = get(groot, 'ScreenSize');
figure_width = 1400; % Figureの幅を拡大 (3列表示のため)
figure_height = 500; % Figureの高さ (1行表示のため)
pos_x = (screen_size(3) - figure_width) / 2;
pos_y = (screen_size(4) - figure_height) / 2;

% Z軸の強調比率
daspect_ratio = [50 50 1];

% --- 4. マスク処理（表示用） ---
% マスクを利用して、脳領域外を NaN に設定

RDF_no_motion_masked = RDF_no_motion;
RDF_no_motion_masked(Mask_no_motion == 0) = NaN;

RDF_with_artifact_masked = RDF_with_artifact;
RDF_with_artifact_masked(Mask_with_artifact == 0) = NaN;

RDF_with_motion_masked = RDF_with_motion;
RDF_with_motion_masked(Mask_with_motion == 0) = NaN;


%% --- 5. 全スライス比較プロット (1スライス/Figure) ---

fprintf('全スライスの 3D (mesh) 比較表示を開始します...\n');

% Z軸の全スライスをループ処理
for slice_idx = 1:total_slices
    
    % Figureを新規作成（今回はスライスごとに1つのFigureを作成）
    fig = figure; 
    set(fig, 'Position', [pos_x, pos_y, figure_width, figure_height]);
    set(fig, 'Color', 'w'); % Figureの背景色を白に
    
    % Figure全体にタイトルを追加
    sgtitle(['3D RDF Comparison - Slice ' num2str(slice_idx)]);

    
    % --- 1. 左側: 体動なし (No Motion) ---
    ax1 = subplot(rows, cols, 1);
    
    data1 = RDF_no_motion_masked(:, :, slice_idx);
    mesh(ax1, data1);
    
    axis tight;
    axis on; 
    daspect(daspect_ratio); 
    colormap(ax1, 'default');
    view(3); 
    title('1. No Motion');
    ylabel(['Slice ' num2str(slice_idx)], 'Visible', 'on', 'FontWeight', 'bold', 'Color', 'k');

    % --- 2. 中央: 体動あり 1 (With Artifact) ---
    ax2 = subplot(rows, cols, 2);
    
    data2 = RDF_with_artifact_masked(:, :, slice_idx);
    mesh(ax2, data2);
    
    axis tight;
    axis on; 
    daspect(daspect_ratio); 
    colormap(ax2, 'default');
    view(3); 
    title('2. With Artifact');
    
    % --- 3. 右側: 体動あり 2 (With Motion) ---
    ax3 = subplot(rows, cols, 3);
    
    data3 = RDF_with_motion_masked(:, :, slice_idx);
    mesh(ax3, data3);

    axis tight;
    axis on;
    daspect(daspect_ratio);
    colormap(ax3, 'default');
    view(3); 
    title('3. With Motion');

    % 全てのプロットに共通のカラーバーを右端に配置
    % *注意: meshプロットの値の範囲が大きく異なると、この共通カラーバーの表示は難しい場合があります。
    %          ここでは、一旦個別のカラーバーを外し、軸は保持しています。
    
end

fprintf('全スライスの 3D (mesh) 3データ比較表示が完了しました。\n');