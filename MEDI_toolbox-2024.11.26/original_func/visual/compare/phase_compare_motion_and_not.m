%% 3Dアンラップ結果の比較プログラム（体動なし vs. 体動あり）
% QSM_processingを使用した結果（iFreq: アンラップ後の位相）を比較します。

% --- 1. 初期設定とパスの定義 ---

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
load_path_no_motion = fullfile(base_path_no_motion, output_folder); % スペース修正
load_path_with_motion = fullfile(base_path_with_motion, output_folder);


% --- 2. データの読み込み ---

% 1. 体動なし (No Motion) データの読み込み
fprintf('体動なし (No Motion) データを読み込み中: %s\n', load_path_no_motion);
try
    % iFreq, RDF, Mask を読み込む想定 (phase.mat, Mask.matに分かれている場合、loadを調整)
    load(fullfile(load_path_no_motion, 'phase.mat'), 'iFreq', 'iFreq_raw'); % 'iFreq' (アンラップ後) を取得
    load(fullfile(load_path_no_motion, 'Mask.mat'), 'Mask'); % 'Mask' を取得 % スペース修正
    iFreq_no_motion = iFreq;
    Mask_no_motion = Mask; % スペース修正
    clear iFreq iFreq_raw Mask; % 変数をクリアして混同を防ぐ
    fprintf('体動なしデータの読み込みが完了しました。\n');
catch ME
    fprintf(2, 'エラー: 体動なしデータの読み込みに失敗しました。\n');
    disp(ME.message);
    return;
end

% 2. 体動あり (Motion) データの読み込み
fprintf('体動あり (Motion) データを読み込み中: %s\n', load_path_with_motion);
try
    load(fullfile(load_path_with_motion, 'phase.mat'), 'iFreq', 'iFreq_raw');
    load(fullfile(load_path_with_motion, 'Mask.mat'), 'Mask');
    iFreq_with_motion = iFreq;
    Mask_with_motion = Mask; % スペース修正
    clear iFreq iFreq_raw Mask;
    fprintf('体動ありデータの読み込みが完了しました。\n');
catch ME
    fprintf(2, 'エラー: 体動ありデータの読み込みに失敗しました。\n');
    disp(ME.message);
    return;
end

% --- 3. データのチェックと初期設定 ---

% データサイズの確認（両データセットでZ軸スライス数が同じであることを前提とする）
if ~isequal(size(iFreq_no_motion), size(iFreq_with_motion))
    warning('体動なしデータと体動ありデータのサイズが異なります。続行できません。');
    return;
end

total_slices = size(iFreq_no_motion, 3); % Z軸の総スライス数を取得
slices_per_figure = 2;         % 1つのFigureに表示する比較セット数 (2スライス分) % スペース修正
rows = slices_per_figure;      % Figureの行数 (2行) % スペース修正
cols = 2;                      % Figureの列数 (2列: No Motion/With Motion) % スペース修正

% 画面サイズを取得してFigureの位置を決定
screen_size = get(groot, 'ScreenSize');
figure_width = 1000; % Figureの幅 (ピクセル単位) % スペース修正
figure_height = 900; % Figureの高さ (ピクセル単位)
pos_x = (screen_size(3) - figure_width) / 2;
pos_y = (screen_size(4) - figure_height) / 2;

% --- 4. マスク処理（表示用） ---
% マスクを利用して、脳領域外を NaN に設定し、meshプロットで非表示にする

% iFreq_no_motion を マスク処理
iFreq_no_motion_masked = iFreq_no_motion;
iFreq_no_motion_masked(Mask_no_motion == 0) = NaN;

% iFreq_with_motion を マスク処理
iFreq_with_motion_masked = iFreq_with_motion;
iFreq_with_motion_masked(Mask_with_motion == 0) = NaN;


%% --- 5. 全スライス比較プロット (2セット/Figure) ---

fprintf('全スライスの 3D (mesh) 比較表示を開始します...\n');

% Z軸の全スライスをループ処理
for slice_idx = 1:total_slices
    
    % 3. Figureを新規作成するタイミングか判定
    if mod(slice_idx - 1, slices_per_figure) == 0
        fig = figure; % スペース修正
        set(fig, 'Position', [pos_x, pos_y, figure_width, figure_height]);
        set(fig, 'Color', 'w'); % Figureの背景色を白に
        
        figure_num = floor((slice_idx - 1) / slices_per_figure) + 1;
        
        % Figure全体にタイトルを追加
        sgtitle(['3D Unwrapped Phase (iFreq) Comparison: No Motion vs. With Motion - Figure ' num2str(figure_num) ...
                 ' (Slices ' num2str(slice_idx) ' - ' ...
                 num2str(min(slice_idx + slices_per_figure - 1, total_slices)) ')']);
    end
    
    % 4. 現在のFigure内での行インデックス (1から2) を計算
    row_idx = mod(slice_idx - 1, slices_per_figure) + 1;
    
    % 5. Subplotの位置を計算
    plot_idx_no_motion = (row_idx - 1) * cols + 1; % 左側の列
    plot_idx_with_motion = (row_idx - 1) * cols + 2; % 右側の列 % スペース修正
    
    % --- 左側: 体動なし (No Motion) の iFreq ---
    ax_no_motion = subplot(rows, cols, plot_idx_no_motion);
    
    % データを取得
    data_no_motion = iFreq_no_motion_masked(:, :, slice_idx);
    
    % meshでプロット
    mesh(ax_no_motion, data_no_motion);
    
    axis tight;
    axis on; % 3Dでは軸をONに
    daspect([50 50 1]); % Z軸を強調 (50 50 1)
    colormap(ax_no_motion, 'default');
    view(3); % 3D表示を有効にする
    % colorbar(ax_no_motion); % カラーバーはスペース節約のため省略可能
    
    % 各行の左側にスライス番号を表示
    ylabel(['Slice ' num2str(slice_idx)], 'Visible', 'on', 'FontWeight', 'bold', 'Color', 'k');
    
    % 各Figureの1行目にだけ列タイトルを表示
    if row_idx == 1
        title('No Motion (iFreq)');
    end

    % --- 右側: 体動あり (With Motion) の iFreq ---
    ax_with_motion = subplot(rows, cols, plot_idx_with_motion);
    
    % データを取得
    data_with_motion = iFreq_with_motion_masked(:, :, slice_idx);
    
    % meshでプロット
    mesh(ax_with_motion, data_with_motion);

    axis tight;
    axis on;
    daspect([50 50 1]); % Z軸を強調 (50 50 1)
    colormap(ax_with_motion, 'default');
    view(3); % 3D表示を有効にする
    % colorbar(ax_with_motion); % スペース修正
    
    if row_idx == 1
        title('With Motion (iFreq)');
    end
    
end

fprintf('全スライスの 3D (mesh) 比較表示が完了しました。\n'); 