%================================================================
% fitting 結果 確認用プログラム　シングルスライス
% 
% 左側: PDF適用前 (iFreq)% 3. 中央: PDF適用後 (RDF)，右側　homo_RDF = iFreq - fitting;
%
% 概要:
%   QSM_processing.m で計算された、PDF適用前の磁場マップ (iFreq) と
%   適用後の局所磁場マップ (RDF)、および多項式補正後のマップ (homo_RDF) を
%   スライスごとに、個別のFigureウィンドウで並べて表示します。
%
% 依存ファイル:
%   - 'phase.mat' (iFreq が必要)
%   - 'PDF.mat' (RDF が必要)
%   - 'Mask.mat' (Mask が必要)
%   - 'fitting.mat' (fitting が必要)
%
% 使い方:
%   1. 必要な .mat ファイルが所定のパスにあることを確認します。
%   2. このスクリプトを実行します。
%================================================================
clear variables;
close all;

%% --- 1. 初期設定とデータ読み込み ---
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
image_file_0 = image_file_00;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_path = fullfile(image_file_0, image_file_3);
load_fitting_path = fullfile(image_file_0, image_file_5);

try
    fprintf('phase.mat を読み込んでいます...\n');
    load(fullfile(load_path, 'phase.mat'), 'iFreq'); % 'iFreq' のみ読み込む
    fprintf('PDF.mat を読み込んでいます...\n');
    load(fullfile(load_path, 'PDF.mat'), 'RDF');     % 'RDF' のみ読み込む
    fprintf('Mask.mat を読み込んでいます...\n');
    load(fullfile(load_path, 'Mask.mat'), 'Mask');   % 'Mask' のみ読み込む
    fprintf('fitting.mat を読み込んでいます...\n');
    load(fullfile(load_fitting_path, 'fitting.mat'), 'fitting'); % 'fitting' のみ読み込む
catch ME
    fprintf('ファイルの読み込み中にエラーが発生しました。\n');
    fprintf('エラー: %s\n', ME.message);
    fprintf('必要な .mat ファイル (phase.mat, PDF.mat, Mask.mat, fitting.mat) がパスに存在するか確認してください。\n');
    return;
end

fprintf('データの読み込みが完了しました。\n');

% 多項式フィッティングによる背景除去結果を計算
homo_RDF = iFreq - fitting;

%% --- 2. マスクの適用 ---
if ~exist('Mask', 'var')
    warning('変数 "Mask" が見つかりません。');
    fprintf('マスクなしで表示を試みますが、背景ノイズも表示されます。\n');
    iFreq_to_show = iFreq;
    RDF_to_show = RDF;
    homo_RDF_show = homo_RDF;
else
    fprintf('脳マスクを適用しています...\n');
    % logical型に変換して適用
    Mask = logical(Mask); 
    iFreq_to_show = iFreq .* Mask;
    RDF_to_show = RDF .* Mask;
    homo_RDF_show = homo_RDF .* Mask;
end


%% --- 3. imshow (2D) での比較 (全スライスループ) ---

total_slices = size(iFreq_to_show, 3);

if total_slices == 0
    fprintf('エラー: 読み込んだデータにスライスがありません。\n');
    return;
end

fprintf('全 %d スライスの 2D (imshow) 比較をスライスごとに表示します...\n', total_slices);

% ★★★ 変更点: 全スライスをループ処理 ★★★
for slice_idx = 1:total_slices
    
    fprintf('スライス %d / %d を表示中...\n', slice_idx, total_slices);

    % 2. スライスごとの新しいFigureを作成
    figure_name = sprintf('2D PDF Comparison - Slice %d / %d', slice_idx, total_slices);
    % 'WindowState', 'maximized' で最大化表示
    figure('Name', figure_name, 'WindowState', 'maximized');
    sgtitle(sprintf('2D PDF Comparison - Slice %d', slice_idx), 'FontWeight', 'bold', 'FontSize', 14);

    % 3. 左側: PDF適用前 (iFreq)
    ax1 = subplot(1, 3, 1);
    img_slice_before = iFreq_to_show(:, :, slice_idx);
    imshow(img_slice_before, []);
    colormap(ax1, 'gray');
    axis on;
    daspect([1,1,1]);
    title('Before PDF (iFreq)');
    xlabel('X Index'); ylabel('Y Index');
    colorbar;

    % 4. 中央: PDF適用後 (RDF)
    ax2 = subplot(1, 3, 2);
    img_slice_after = RDF_to_show(:, :, slice_idx);
    imshow(img_slice_after, []);
    colormap(ax2, 'gray');
    axis on;
    daspect([1,1,1]);
    title('After PDF (RDF)');
    xlabel('X Index');
    colorbar;

    % 5. 右側: 多項式補正後 (homo_RDF)
    ax3 = subplot(1, 3, 3);
    % (変数名を img_slice_homo に変更)
    img_slice_homo = homo_RDF_show(:, :, slice_idx);
    imshow(img_slice_homo, []);
    colormap(ax3, 'gray');
    axis on;
    daspect([1,1,1]);
    title('PolyFit Corrected (homo_RDF)');
    xlabel('X Index');
    colorbar;
    
    % (オプション) ウィンドウの更新を強制
    drawnow;
    
end % --- スライスループの終了 ---

fprintf('全スライスの表示が完了しました。\n');
