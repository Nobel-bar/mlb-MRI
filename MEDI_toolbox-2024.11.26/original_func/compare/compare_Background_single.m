%================================================================
% PDF適用前後・背景磁場 3項目比較プログラム
%
% 概要:
%   以下の3つのマップを2Dおよび3Dで並べて表示します。
%   1. iFreq (PDF適用前)
%   2. iFreq_bg (除去された背景磁場 = iFreq - RDF)
%   3. RDF (PDF適用後)
%
% 依存ファイル:
%   - 'phase.mat' (iFreq を含む)
%   - 'PDF.mat' (RDF を含む)
%   - 'Mask.mat' (Mask を含む)
%
% 使い方:
%   - このスクリプトを実行します。
%================================================================
clear variables;
close all;


%% --- 1. データ読み込み ---
% (パスはご自身の環境に合わせて設定してください)
image_file_1 = 'F:/hamaguchi/copy/20241205_RawData_H/Volunteer_Rotate_H/2DGE_0deg_H/total_slice';
image_file_2 = 'output_data';

save_path = fullfile(image_file_1, image_file_2);

fprintf('データを読み込んでいます...\n');
try
    load(fullfile(save_path, 'phase.mat'), 'iFreq');
    load(fullfile(save_path, 'PDF.mat'), 'RDF');
    load(fullfile(save_path, 'Mask.mat'), 'Mask');
catch ME
    fprintf('ファイルの読み込みに失敗しました。\n');
    fprintf('phase.mat, PDF.mat, Mask.mat が %s に存在するか確認してください。\n', save_path);
    rethrow(ME);
end

fprintf('データの読み込みが完了しました。\n');


%% 2. マスクの適用
% ★変更点: iFreq_to_show を追加
if ~exist('Mask', 'var')
    warning('変数 "Mask" が見つかりません。');
    fprintf('マスクなしで表示を試みますが、背景ノイズも表示されます。\n');
    iFreq_to_show = iFreq;
    RDF_to_show = RDF;
    iFreq_bg_to_show = iFreq - RDF;
else
    fprintf('脳マスクを適用しています...\n');
    iFreq_to_show = iFreq .* Mask; % ★追加
    RDF_to_show = RDF .* Mask;
    iFreq_bg_to_show = (iFreq - RDF) .* Mask;
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

% 2. 比較用の新しいFigureを作成 (最大化)
figure('Name', '2D PDF 3-Way Comparison (Maximized)', 'WindowState', 'maximized');
sgtitle(sprintf('2D PDF 3-Way Comparison - Slice %d', slice_to_display), 'FontWeight', 'bold');

% ★変更点: 1行3列レイアウトに変更

% 3. 左側: PDF適用前 (iFreq)
ax1 = subplot(1, 3, 1);
img_slice_ifreq = iFreq_to_show(:, :, slice_to_display); % ★変数追加
imshow(img_slice_ifreq, []);
colormap(ax1, 'gray');
axis on;
daspect([1,1,1]);
title('Before PDF (iFreq)');
xlabel('X Index'); ylabel('Y Index');
colorbar;

% 4. 中央: Background (iFreq_bg_to_show)
ax2 = subplot(1, 3, 2); % ★インデックス変更 (1,2,1) -> (1,3,2)
img_slice_bg = iFreq_bg_to_show(:, :, slice_to_display); % ★変数名変更
imshow(img_slice_bg, []);
colormap(ax2, 'gray');
axis on;
daspect([1,1,1]);
title('Background (iFreq_bg)');
xlabel('X Index');
colorbar;

% 5. 右側: PDF適用後 (RDF)
ax3 = subplot(1, 3, 3); % ★インデックス変更 (1,2,2) -> (1,3,3)
img_slice_rdf = RDF_to_show(:, :, slice_to_display); % ★変数名変更
imshow(img_slice_rdf, []);
colormap(ax3, 'gray');
axis on;
daspect([1,1,1]);
title('After PDF (RDF)');
xlabel('X Index');
colorbar;

% ★変更点: 3プロット用の 'Position' に調整
ax1.Position = [0.05, 0.05, 0.28, 0.85];
ax2.Position = [0.36, 0.05, 0.28, 0.85];
ax3.Position = [0.67, 0.05, 0.28, 0.85];


%% --- 4. mesh (3D) での比較 ---

fprintf('スライス %d の 3D (mesh) 比較を表示します。\n', slice_to_display);

% 1. メッシュ比較用の新しいFigureを作成 (最大化)
figure('Name', '3D Mesh 3-Way Comparison (Maximized)', 'WindowState', 'maximized');
sgtitle(sprintf('3D Mesh 3-Way Comparison - Slice %d', slice_to_display), 'FontWeight', 'bold');

% ★変更点: 1行3列レイアウトに変更

% 2. 左側: PDF適用前 (iFreq)
ax_m1 = subplot(1, 3, 1);
img_slice_ifreq_mesh = img_slice_ifreq; % 2Dプロットのデータを利用
img_slice_ifreq_mesh(img_slice_ifreq_mesh == 0) = NaN; % 0の値を非表示に

mesh(ax_m1, img_slice_ifreq_mesh);
axis tight;
daspect([1,1,1/50]); % daspect([50 50 1]) と同じ
axis on;
colormap(ax_m1, 'default');
xlabel('X Index');
ylabel('Y Index');
zlabel('Field Map (a.u.)');
title('Before PDF (iFreq)');
colorbar;

% 3. 中央: Background (iFreq_bg_to_show)
ax_m2 = subplot(1, 3, 2); % ★インデックス変更 (1,2,1) -> (1,3,2)
img_slice_bg_mesh = img_slice_bg; % 2Dプロットのデータを利用
img_slice_bg_mesh(img_slice_bg_mesh == 0) = NaN; % 0の値を非表示に

mesh(ax_m2, img_slice_bg_mesh);
axis tight;
daspect([1,1,1/50]);
axis on;
colormap(ax_m2, 'default');
xlabel('X Index');
ylabel('Y Index');
zlabel('Field Map (a.u.)');
title('Background (iFreq_bg)');
colorbar;

% 4. 右側: PDF適用後 (RDF)
ax_m3 = subplot(1, 3, 3); % ★インデックス変更 (1,2,2) -> (1,3,3)
img_slice_rdf_mesh = img_slice_rdf; % 2Dプロットのデータを利用
img_slice_rdf_mesh(img_slice_rdf_mesh == 0) = NaN; % 0の値を非表示に

mesh(ax_m3, img_slice_rdf_mesh);
axis tight;
daspect([1,1,1/50]);
axis on;
colormap(ax_m3, 'default');
xlabel('X Index');
ylabel('Y Index');
zlabel('Local Field (a.u.)');
title('After PDF (RDF)');
colorbar;
% (オプション) RDFはZ軸のスケールが他と大きく異なる場合、zlimで固定
% zlim(ax_m3, [-0.5 0.5]);


% ★変更点: 3プロット用の 'Position' に調整
ax_m1.Position = [0.05, 0.05, 0.28, 0.85];
ax_m2.Position = [0.36, 0.05, 0.28, 0.85];
ax_m3.Position = [0.67, 0.05, 0.28, 0.85];

fprintf('表示が完了しました。\n');

% (コメントアウトされていたローカル関数部分は削除しました)