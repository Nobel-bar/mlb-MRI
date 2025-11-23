%================================================================
% fitting 結果 確認用プログラム　シングルスライス
% 
% 左側: PDF適用前 (iFreq)% 3. 中央: PDF適用後 (RDF)，右側　homo_RDF = iFreq - fitting;
%
% 使い方:
%   1. QSM_processing.m を実行し、'RDF.mat' を生成します。
%   2. このスクリプトを 'RDF.mat' と同じディレクトリに保存して実行します。
%================================================================
clear variables;
close all;


%% --- 1. データ読み込み ---
%% --- 1. 初期設定 ---
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

load(fullfile(load_path, 'phase.mat'));
load(fullfile(load_path, 'PDF.mat'));
load(fullfile(load_path, 'Mask.mat'));
load(fullfile(load_fitting_path, 'fitting.mat'));

fprintf('読み込んでいます...');
% PDFの入力(iFreq), 出力(RDF), およびマスク(Mask)を読み込む
%  load(data_file, 'iFreq', 'RDF', 'Mask');

fprintf('データの読み込みが完了しました。\n');

homo_RDF = iFreq - fitting;


%% 2. マスクの適用
if ~exist('Mask', 'var')
    warning('変数 "Mask" が RDF.mat 内に見つかりません。');
    fprintf('マスクなしで表示を試みますが、背景ノイズも表示されます。\n');
    iFreq_to_show = iFreq;
    RDF_to_show = RDF;
    homo_RDF_show = homo_RDF;
else
    fprintf('脳マスクを適用しています...\n');
    iFreq_to_show = iFreq .* Mask;
    RDF_to_show = RDF .* Mask;
    homo_RDF_show = homo_RDF .* Mask;
end


%% --- 3. imshow (2D) での比較 ---

% 1. スライス番号の指定
total_slices = size(iFreq_to_show, 3);
slice_to_display = round(total_slices / 2); % 中央スライス

if slice_to_display < 1 || slice_to_display > total_slices
    fprintf('エラー: スライス番号 %d は無効です。1から %d の間で指定してください。\n', ...
            slice_to_display, total_slices);
    return;
end

fprintf('スライス %d の 2D (imshow) 比較を表示します。\n', slice_to_display);

% 2. 比較用の新しいFigureを作成
% --- 変更点: 'WindowState', 'maximized' を追加 ---
figure('Name', '2D PDF Comparison (Maximized)', 'WindowState', 'maximized');
sgtitle(sprintf('2D PDF Comparison - Slice %d', slice_to_display), 'FontWeight', 'bold');

% 3. 左側: PDF適用前 (iFreq)
% --- 変更点: subplotのハンドル(ax1)を取得 ---
ax1 = subplot(1, 3, 1);
img_slice_before = iFreq_to_show(:, :, slice_to_display);
imshow(img_slice_before, []);
colormap(ax1, 'gray');
axis on;
daspect([1,1,1]);
title('Before PDF (iFreq)');
xlabel('X Index'); ylabel('Y Index');
colorbar;

% 4. 右側: PDF適用後 (RDF)
% --- 変更点: subplotのハンドル(ax2)を取得 ---
ax2 = subplot(1, 3, 2);
img_slice_after = RDF_to_show(:, :, slice_to_display);
imshow(img_slice_after, []);
colormap(ax2, 'gray');
axis on;
daspect([1,1,1]);
title('After PDF (RDF)');
xlabel('X Index');
colorbar;

ax3 = subplot(1, 3, 3);
img_slice_tomorrow = homo_RDF_show(:, :, slice_to_display);
imshow(img_slice_tomorrow, []);
colormap(ax3, 'gray');
axis on;
daspect([1,1,1]);
title('homo (homo_RDF)');
xlabel('X Index'); ylabel('Y Index');
colorbar;

%% --- 4. mesh (3D) での比較 ---

fprintf('スライス %d の 3D (mesh) 比較を表示します。\n', slice_to_display);

% 1. メッシュ比較用の新しいFigureを作成
% --- 変更点: 'WindowState', 'maximized' を追加 ---
figure('Name', '3D Mesh Comparison (Maximized)', 'WindowState', 'maximized');
sgtitle(sprintf('3D Mesh Comparison - Slice %d', slice_to_display), 'FontWeight', 'bold');

% 2. 左側: PDF適用前 (iFreq)
% --- 変更点: subplotのハンドル(ax3)を取得 ---
ax4 = subplot(1, 3, 1);
img_slice_before_mesh = img_slice_before;
img_slice_before_mesh(img_slice_before_mesh == 0) = NaN; % 0の値を非表示に

mesh(ax4, img_slice_before_mesh);
axis tight;
daspect([1,1,1/50]);
axis on;
colormap(ax4, 'default');
xlabel('X Index');
ylabel('Y Index');
zlabel('Field Map (a.u.)');
title('Before PDF (iFreq)');
colorbar;

% 3. 右側: PDF適用後 (RDF)
% --- 変更点: subplotのハンドル(ax4)を取得 ---
ax5 = subplot(1, 3, 2);
img_slice_after_mesh = img_slice_after;
img_slice_after_mesh(img_slice_after_mesh == 0) = NaN; % 0の値を非表示に

mesh(ax5, img_slice_after_mesh);
axis tight;
daspect([1,1,1/50]);
axis on;
colormap(ax5, 'default');
xlabel('X Index');
ylabel('Y Index');
zlabel('Local Field (a.u.)');
title('After PDF (RDF)');
colorbar;

% --- 変更点: subplotのハンドル(ax4)を取得 ---
ax6 = subplot(1, 3, 3);
img_slice_tomorrow_mesh = img_slice_tomorrow;
img_slice_tomorrow_mesh(img_slice_tomorrow_mesh == 0) = NaN; % 0の値を非表示に

mesh(ax6, img_slice_tomorrow_mesh);
axis tight;
daspect([1,1,1/50]);
axis on;
colormap(ax6, 'default');
xlabel('X Index');
ylabel('Y Index');
zlabel('Local Field (a.u.)');
title('homo (homo_RDF)');
colorbar;
