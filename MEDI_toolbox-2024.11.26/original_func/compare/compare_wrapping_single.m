%================================================================
% PDF (背景磁場除去) 結果 確認用プログラム
% 
% 概要:
%   QSM_processing.m で計算された、wrap適用前の磁場マップ (iFreq_raw) と
%   適用後の局所磁場マップ (iFreq) を2Dおよび3Dで並べて表示します。
%
% 依存ファイル:
%   - QSM_processing.m を実行した際に生成される '.mat' ファイル
%
% 使い方:
%   1. QSM_processing.m を実行します。
%================================================================

%% --- 1. データ読み込み ---
% QSM_processing.m で保存された RDF.mat ファイルを読み込みます。

% --- 1. 初期設定 ---
image_file_1 = 'F:/hamaguchi/copy/20241205_RawData_H/Volunteer_Rotate_H/2DGE_0deg_H/total_slice';
image_file_2 = 'output_data';

save_path = fullfile(image_file_1, image_file_2);

load(fullfile(save_path, 'phase.mat'));
load(fullfile(save_path, 'Mask.mat'));

fprintf('データの読み込みが完了しました。\n');


%% 2. マスクの適用
if ~exist('Mask', 'var')
    warning('変数 "Mask" が RDF.mat 内に見つかりません。');
    fprintf('マスクなしで表示を試みますが、背景ノイズも表示されます。\n');
else
    fprintf('脳マスクを適用しています...\n');
    iFreq_raw = iFreq_raw .* Mask;
    iFreq = iFreq .* Mask;
end


%% --- 3. imshow (2D) での比較 ---

% 1. スライス番号の指定
total_slices = size(iFreq_raw, 3);
slice_to_display = round(total_slices / 2); % 中央スライス

if slice_to_display < 1 || slice_to_display > total_slices
    fprintf('エラー: スライス番号 %d は無効です。1から %d の間で指定してください。\n', ...
            slice_to_display, total_slices);
    return;
end

fprintf('スライス %d の 2D (imshow) 比較を表示します。\n', slice_to_display);

% 2. 比較用の新しいFigureを作成
% --- 変更点: 'WindowState', 'maximized' を追加 ---
figure('Name', '2D Wrapping Comparison (Maximized)', 'WindowState', 'maximized');
sgtitle(sprintf('2D Wrapping Comparison - Slice %d', slice_to_display), 'FontWeight', 'bold');

% 3. 左側: Wrap適用前 (iFreq_raw)
% --- 変更点: subplotのハンドル(ax1)を取得 ---
ax1 = subplot(1, 2, 1);
img_slice_before = iFreq_raw(:, :, slice_to_display);
imshow(img_slice_before, []);
colormap(ax1, 'gray');
axis on;
daspect([1,1,1]);
title('Before Wrapping (iFreq_raw)');
xlabel('X Index'); ylabel('Y Index');
colorbar;

% 4. 右側: Wrappig適用後 (iFreq)
% --- 変更点: subplotのハンドル(ax2)を取得 ---
ax2 = subplot(1, 2, 2);
img_slice_after = iFreq(:, :, slice_to_display);
imshow(img_slice_after, []);
colormap(ax2, 'gray');
axis on;
daspect([1,1,1]);
title('After Wrapping (iFreq)');
xlabel('X Index');
colorbar;

% --- 変更点: サブプロットの位置を調整して余白を削減 ---
% [left, bottom, width, height] (0から1の正規化座標)
% 上部のsgtitleの分を考慮しつつ、左右と下の余白を詰めます
ax1.Position = [0.05, 0.05, 0.43, 0.85];
ax2.Position = [0.52, 0.05, 0.43, 0.85];


%% --- 4. mesh (3D) での比較 ---

fprintf('スライス %d の 3D (mesh) 比較を表示します。\n', slice_to_display);

% 1. メッシュ比較用の新しいFigureを作成
% --- 変更点: 'WindowState', 'maximized' を追加 ---
figure('Name', '3D Mesh Comparison (Maximized)', 'WindowState', 'maximized');
sgtitle(sprintf('3D Mesh Comparison - Slice %d', slice_to_display), 'FontWeight', 'bold');

% 2. 左側: Wrap適用前 (iFreq_raw)
% --- 変更点: subplotのハンドル(ax3)を取得 ---
ax3 = subplot(1, 2, 1);
img_slice_before_mesh = img_slice_before;
img_slice_before_mesh(img_slice_before_mesh == 0) = NaN; % 0の値を非表示に

mesh(ax3, img_slice_before_mesh);
axis tight;
daspect([1,1,1/50]);
axis on;
colormap(ax3, 'default');
xlabel('X Index');
ylabel('Y Index');
zlabel('Field Map (a.u.)');
title('Before Wrapping (iFreq_raw)');
colorbar;

% 3. 右側: Wrap適用後 (iFreq)
% --- 変更点: subplotのハンドル(ax4)を取得 ---
ax4 = subplot(1, 2, 2);
img_slice_after_mesh = img_slice_after;
img_slice_after_mesh(img_slice_after_mesh == 0) = NaN; % 0の値を非表示に

mesh(ax4, img_slice_after_mesh);
axis tight;
daspect([1,1,1/50]);
axis on;
colormap(ax4, 'default');
xlabel('X Index');
ylabel('Y Index');
zlabel('Local Field (a.u.)');
title('After Wrapping (iFreq)');
colorbar;

% --- 変更点: サブプロットの位置を調整して余白を削減 ---
ax3.Position = [0.05, 0.05, 0.43, 0.85];
ax4.Position = [0.52, 0.05, 0.43, 0.85];


fprintf('表示が完了しました。\n');ax4.Position = [left_margin_2, bottom_margin, plot_width, plot_height];


%{

% --- 1. 初期設定 ---
image_file_1 = 'Volunteer_Rotate_H';
image_file_2 = '2DGE_0deg_H/total_slice';

save_path = fullfile('.', image_file_1, image_file_2);
save_raw_data(fullfile(save_path, 'phase_before.raw'), img_slice_raw);
save_raw_data(fullfile(save_path, 'phase_after.raw'), img_slice_unwrapped);


% -------------------------------------------------------------------
% スクリプトの最後にローカル関数を定義します
% -------------------------------------------------------------------
function save_raw_data(filepath, data)
    fid = fopen(filepath, 'w');
    if fid == -1
        error('ファイルが開けませんでした: %s', filepath);
    end
    fwrite(fid, data, 'double');
    fclose(fid);
end
%}