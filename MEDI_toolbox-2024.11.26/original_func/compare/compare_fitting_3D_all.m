%================================================================
% PDF (背景磁場除去) 結果 確認用プログラム (全スライス 3D Mesh 表示版)
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
%================================================================
clear variables;
close all;

%% --- 1. 初期設定とデータ読み込み ---
fprintf('1. パラメータを設定し、データを読み込んでいます...\n');

% パス設定
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_000 = "C:\Users\hamaguchi\Downloads\matlab";
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

image_file_0 = image_file_000; % local用

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


%% --- 3. 3D 比較 (全スライスループ) ---

total_slices = size(iFreq_to_show, 3);

if total_slices == 0
    fprintf('エラー: 読み込んだデータにスライスがありません。\n');
    return;
end

fprintf('全 %d スライスの 3D (mesh) 比較をスライスごとに表示します...\n', total_slices);

% ★★★ 全スライスをループ処理 ★★★
for slice_idx = 7:total_slices-7
    
    fprintf('スライス %d / %d を表示中...\n', slice_idx, total_slices);

    % --- 3.1. mesh (3D) での比較 ---
    
    % 3D表示用のスライスデータを取得
    img_slice_before = iFreq_to_show(:, :, slice_idx);
    img_slice_after = RDF_to_show(:, :, slice_idx);
    img_slice_homo = homo_RDF_show(:, :, slice_idx);

    % 3D用の新しいFigureを作成
    figure_name_3D = sprintf('3D Mesh Comparison - Slice %d / %d', slice_idx, total_slices);
    figure('Name', figure_name_3D, 'WindowState', 'maximized');
    sgtitle(sprintf('3D Mesh Comparison - Slice %d', slice_idx), 'FontWeight', 'bold', 'FontSize', 14);

    % 左側: PDF適用前 (iFreq)
    ax4 = subplot(1, 3, 1);
    img_slice_before_mesh = img_slice_before;
    img_slice_before_mesh(img_slice_before_mesh == 0) = NaN; % 0の値を非表示に

    mesh(ax4, img_slice_before_mesh);
    axis tight;
    daspect([1,1,1/50]); % Z軸のスケールを調整
    axis on;
    colormap(ax4, 'default');
    xlabel('X Index');
    ylabel('Y Index');
    zlabel('Field Map (a.u.)');
    title('Before PDF (iFreq)');
    colorbar;

    % 中央: PDF適用後 (RDF)
    ax5 = subplot(1, 3, 2);
    img_slice_after_mesh = img_slice_after;
    img_slice_after_mesh(img_slice_after_mesh == 0) = NaN; % 0の値を非表示に

    mesh(ax5, img_slice_after_mesh);
    axis tight;
    daspect([1,1,1/50]);
    axis on;
    colormap(ax5, 'turbo');
    xlabel('X Index');
    ylabel('Y Index');
    zlabel('Local Field (a.u.)');
    title('After PDF (RDF)');
    colorbar;

    % 右側: 多項式補正後 (homo_RDF)
    ax6 = subplot(1, 3, 3);
    % (★変数名を img_slice_homo に修正)
    img_slice_homo_mesh = img_slice_homo;
    img_slice_homo_mesh(img_slice_homo_mesh == 0) = NaN; % 0の値を非表示に

    mesh(ax6, img_slice_homo_mesh);
    axis tight;
    daspect([1,1,1/50]);
    axis on;
    colormap(ax6, 'turbo');
    xlabel('X Index');
    ylabel('Y Index');
    zlabel('Local Field (a.u.)');
    % (★タイトルを 2D 側と統一)
    title('PolyFit Corrected (homo_RDF)');
    colorbar;
    
    % (オプション) 3Dウィンドウの更新を強制
    drawnow;
    
end % --- スライスループの終了 ---

fprintf('全スライスの表示が完了しました。\n');

