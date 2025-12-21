%================================================================
% これはただの位相の比較 結果 確認用プログラム (マスク収縮版)
%================================================================

%% --- 1. データ読み込み ---
path_mat_GE_Mag   = "F:\hamaguchi\20251204\2DGE\3_output_data\Mask.mat";
path_mat_GE_Phase = "F:\hamaguchi\20251204\2DGE\3_output_data\phase.mat";
path_mat_SE_Mag   = "F:\hamaguchi\20251204\2DSE\3_output_data\Mask.mat";
path_mat_SE_Phase = "F:\hamaguchi\20251204\2DSE\3_output_data\phase.mat";

unit_str = 'Hz';      % 単位定義
erosion_size = 10;     % ★変更: マスクを削るピクセル数 (大きいほど内側まで削れます)

%% 2. データ読み込み
fprintf('データを読み込み中...\n');

[img_mag_GE, img_mask_GE, img_freq_GE] = load_mat_data(path_mat_GE_Mag, path_mat_GE_Phase);
[img_mag_SE, img_mask_SE, img_freq_SE] = load_mat_data(path_mat_SE_Mag, path_mat_SE_Phase);

% 3Dデータ処理 (中央スライス抽出)
if ndims(img_freq_GE) == 3
    sliceIdx = round(size(img_freq_GE, 3) / 2);
    
    % GEデータのスライス
    if ~isempty(img_mag_GE), img_mag_GE = img_mag_GE(:,:,sliceIdx); end
    if ~isempty(img_mask_GE), img_mask_GE = img_mask_GE(:,:,sliceIdx); end
    img_freq_GE = img_freq_GE(:,:,sliceIdx);
    
    % SEデータのスライス
    if ~isempty(img_mag_SE), img_mag_SE = img_mag_SE(:,:,sliceIdx); end
    if ~isempty(img_mask_SE), img_mask_SE = img_mask_SE(:,:,sliceIdx); end
    img_freq_SE = img_freq_SE(:,:,sliceIdx);
end


%% 2. マスクの適用 (収縮処理を追加)

% 構造要素の作成（円形に削ります）
se = strel('disk', erosion_size); 

% --- GEデータのマスク適用 ---
if isempty(img_mask_GE)
    fprintf('GE: マスクが見つからないため、そのまま表示します。\n');
else
    fprintf('GE: マスク境界を %d ピクセル削って適用しています...\n', erosion_size);
    % マスクを収縮 (論理値に変換して処理 -> doubleに戻す)
    img_mask_GE = imerode(logical(img_mask_GE), se);
    img_mask_GE = double(img_mask_GE);
    
    % 適用
    img_freq_GE = img_freq_GE .* img_mask_GE;
end

% --- SEデータのマスク適用 ---
if isempty(img_mask_SE)
    fprintf('SE: マスクが見つからないため、そのまま表示します。\n');
else
    fprintf('SE: マスク境界を %d ピクセル削って適用しています...\n', erosion_size);
    % マスクを収縮
    img_mask_SE = imerode(logical(img_mask_SE), se);
    img_mask_SE = double(img_mask_SE);
    
    % 適用
    img_freq_SE = img_freq_SE .* img_mask_SE;
end


%% --- 3. imshow (2D) での比較 ---

fprintf('スライス %d の 2D (imshow) 比較を表示します。\n', sliceIdx);

figure('Name', '2D Wrapping Comparison (Maximized)', 'WindowState', 'maximized');
sgtitle(sprintf('2D Wrapping Comparison - Slice %d (Eroded Mask)', sliceIdx), 'FontWeight', 'bold');

% --- 左側: (GE) ---
ax1 = subplot(1, 2, 1);
imshow(img_freq_GE, []);
colormap(ax1, 'gray');
axis on; daspect([1,1,1]);
title('GREのB0磁場');
xlabel('X Index'); ylabel('Y Index');
h = colorbar;
ylabel(h, unit_str, 'FontSize', 12, 'FontWeight', 'bold');

% --- 右側: (SE) ---
ax2 = subplot(1, 2, 2);
imshow(img_freq_SE, []);
colormap(ax2, 'gray');
axis on; daspect([1,1,1]);
title('SEのB0磁場');
xlabel('X Index');
h = colorbar;
ylabel(h, unit_str, 'FontSize', 12, 'FontWeight', 'bold');

% 位置調整
ax1.Position = [0.05, 0.05, 0.43, 0.85];
ax2.Position = [0.52, 0.05, 0.43, 0.85];


%% --- 4. mesh (3D) での比較 ---

fprintf('スライス %d の 3D (mesh) 比較を表示します。\n', sliceIdx);

figure('Name', '3D Mesh Comparison (Maximized)', 'WindowState', 'maximized');
sgtitle(sprintf('3D Mesh Comparison - Slice %d (Eroded Mask)', sliceIdx), 'FontWeight', 'bold');

% --- 左側: GE ---
ax3 = subplot(1, 2, 1);
img_slice_before_mesh = img_freq_GE;
img_slice_before_mesh(img_slice_before_mesh == 0) = NaN; 

mesh(ax3, img_slice_before_mesh);
axis tight; daspect([1,1,1/50]); axis on;
colormap(ax3, 'default');
xlabel('X Index'); ylabel('Y Index');
zlabel(unit_str, 'FontSize', 14, 'FontWeight', 'bold');
title('GRE'); colorbar;

% --- 右側: SE ---
ax4 = subplot(1, 2, 2);
img_slice_after_mesh = img_freq_SE;
img_slice_after_mesh(img_slice_after_mesh == 0) = NaN; 

mesh(ax4, img_slice_after_mesh);
axis tight; daspect([1,1,1/50]); axis on;
colormap(ax4, 'default');
xlabel('X Index'); ylabel('Y Index');
zlabel(unit_str, 'FontSize', 14, 'FontWeight', 'bold');
title('SE'); colorbar;

% 位置調整
ax3.Position = [0.05, 0.05, 0.43, 0.85];
ax4.Position = [0.52, 0.05, 0.43, 0.85];

fprintf('完了しました。\n');


%% 関数群
function [m, mask, f] = load_mat_data(fm, fp)
    m = []; mask = []; f = [];
    if exist(fm, 'file')
        d = load(fm);
        if isfield(d, 'iMag'), m = double(d.iMag); end
        if isfield(d, 'Mask'), mask = double(d.Mask);
        elseif isfield(d, 'mask'), mask = double(d.mask); end
    end
    if exist(fp, 'file')
        d = load(fp);
        if isfield(d, 'iFreq'), f = double(d.iFreq);
        elseif isfield(d, 'phase'), f = double(d.phase); end
    end
end