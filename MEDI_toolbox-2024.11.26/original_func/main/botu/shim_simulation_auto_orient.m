%% Shim Simulation: Auto-Orientation (座標軸自動マッチング版)
% シム電流値からシムを再現しようとするコード
clear; clc; close all;

%% 1. 設定：ファイルパス
% -------------------------------------------------------------------------
dicomFolder_ON    = 'F:\hamaguchi\20251215\dual_echo\27\1_original_data';
path_mat_ON_Mag   = 'F:\hamaguchi\20251215\dual_echo\27\3_qsm_data\Mask.mat';
path_mat_ON_Phase = 'F:\hamaguchi\20251215\dual_echo\27\3_qsm_data\phase.mat';

dicomFolder_OFF    = 'F:\hamaguchi\20251215\dual_echo\27Z\1_original_data';
path_mat_OFF_Mag   = 'F:\hamaguchi\20251215\dual_echo\27Z\3_qsm_data\Mask.mat';
path_mat_OFF_Phase = 'F:\hamaguchi\20251215\dual_echo\27Z\3_qsm_data\phase.mat';
% -------------------------------------------------------------------------

%% 2. データ読み込み
fprintf('データを読み込み中...\n');
[img_mag_ON, img_freq_ON]   = load_mat_data(path_mat_ON_Mag, path_mat_ON_Phase);
[img_mag_OFF, img_freq_OFF] = load_mat_data(path_mat_OFF_Mag, path_mat_OFF_Phase);

% 3Dデータ処理
if ndims(img_mag_ON) == 3
    sliceIdx = round(size(img_mag_ON, 3) / 2);
    img_mag_ON = img_mag_ON(:,:,sliceIdx); img_freq_ON = img_freq_ON(:,:,sliceIdx);
    img_mag_OFF = img_mag_OFF(:,:,sliceIdx); img_freq_OFF = img_freq_OFF(:,:,sliceIdx);
end

% 実測差分
dField_measured = img_freq_ON - img_freq_OFF;
mask = (img_mag_ON > max(img_mag_ON(:))*0.1) & (img_mag_OFF > max(img_mag_OFF(:))*0.1);
dField_measured(~mask) = NaN;

%% 3. DICOM情報抽出
[shim_ON,  info_ON]  = get_dicom_info(dicomFolder_ON);
[shim_OFF, info_OFF] = get_dicom_info(dicomFolder_OFF);








% スライス位置とFOV情報
try
    if isfield(info_ON, 'PixelSpacing')
        dy = info_ON.PixelSpacing(1); dx = info_ON.PixelSpacing(2);
    else, dy=1; dx=1; end
    
    if isfield(info_ON, 'SliceLocation'), z_mm = info_ON.SliceLocation;
    else, z_mm = 0; end
    
    fprintf('  Pixel Spacing : %.2f x %.2f mm\n', dx, dy);
    fprintf('  Slice Location: %.2f mm\n', z_mm);
catch
    dx=1; dy=1; z_mm=0;
end

% シム電流差分
currents = struct();
tags = {'X','Y','Z','XY','XZ','YZ','X2Y2','Z2'}; 
fprintf('\n--- Delta Shim Currents ---\n');
for i = 1:length(tags)
    t = tags{i};
    v_on = 0; if isfield(shim_ON,t), v_on = shim_ON.(t); end
    v_off = 0; if isfield(shim_OFF,t), v_off = shim_OFF.(t); end
    currents.(t) = v_on - v_off;
    fprintf(' %s: %d\n', t, currents.(t));
end

%% 4. 正規化座標系でのマッチング (Normalized Coordinates)
[ny, nx] = size(dField_measured);

% FOVの半径を基準にする (Normalization Radius)
fov_x = nx * dx;
fov_y = ny * dy;
R0 = max(fov_x, fov_y) / 2; % 正規化半径

% 座標グリッド (-R0 ~ +R0)
x_vec = ((1:nx) - nx/2) * dx;
y_vec = ((1:ny) - ny/2) * dy;
[xx_mm, yy_mm] = meshgrid(x_vec, y_vec);

% 正規化座標 (-1 ~ +1)
xx_norm = xx_mm / R0;
yy_norm = yy_mm / R0;
zz_norm = ones(size(xx_norm)) * (z_mm / R0);

idx = find(mask);
vals_measured = dField_measured(idx);

% パターン探索
patterns = {
    'Normal',          xx_norm,       yy_norm;
    'Swap XY',         yy_norm,       xx_norm;
    'Rot 90',         -yy_norm,       xx_norm;
    'Rot 180',        -xx_norm,      -yy_norm;
    'Rot 270',         yy_norm,      -xx_norm; % 前回これが有力
    'Flip X',         -xx_norm,       yy_norm;
    'Flip Y',          xx_norm,      -yy_norm;
    'Swap & Flip',    -yy_norm,      -xx_norm; 
};

best_R = -1;
best_pattern_idx = 0;
best_shape = [];

fprintf('\n--- 座標軸マッチング探索 (Normalized) ---\n');

for p = 1:size(patterns, 1)
    name = patterns{p, 1};
    X_map = patterns{p, 2}; 
    Y_map = patterns{p, 3}; 
    
    x_r = X_map(idx);
    y_r = Y_map(idx);
    z_r = zz_norm(idx);
    
    % Basis Functions (Normalized [-1, 1])
    % これにより 1次項(x) と 2次項(x^2) のオーダーが揃います
    
    term_X  = x_r;
    term_Y  = y_r;
    term_Z  = z_r;
    
    term_XY = x_r .* y_r;
    term_XZ = x_r .* z_r;
    term_YZ = y_r .* z_r;
    term_X2Y2 = x_r.^2 - y_r.^2;
    term_Z2   = z_r.^2 - 0.5 * (x_r.^2 + y_r.^2);
    
    % Shape calculation
    shape_vec = ...
        currents.X * term_X + ...
        currents.Y * term_Y + ...
        currents.Z * term_Z + ...
        currents.XY * term_XY + ...
        currents.XZ * term_XZ + ...
        currents.YZ * term_YZ + ...
        currents.X2Y2 * term_X2Y2 + ...
        currents.Z2 * term_Z2;
    
    if std(shape_vec)==0
        corr_r = 0;
    else
        c = corrcoef(vals_measured, shape_vec);
        corr_r = c(1,2);
    end
    
    fprintf('  Pattern "%s": R = %.4f\n', name, corr_r);
    
    if abs(corr_r) > best_R
        best_R = abs(corr_r);
        best_pattern_idx = p;
        best_shape = shape_vec;
    end
end

fprintf('\n>>> 最適パターン: %s (R = %.4f)\n', patterns{best_pattern_idx, 1}, best_R);

%% 5. 最終フィッティングと表示
poly = polyfit(best_shape, vals_measured, 1);
alpha = poly(1);
beta = poly(2);

sim_map = zeros(ny, nx);
sim_map(idx) = alpha * best_shape + beta;
sim_map(~mask) = NaN;

fprintf('推定感度: %.4e (Hz/DAC)\n', alpha);

% 可視化
figure('Name', 'Shim Sim Normalized', 'Color', 'w', 'Position', [100, 100, 1200, 500]);
% オートスケール範囲（外れ値除去）
c_min = quantile(vals_measured, 0.05);
c_max = quantile(vals_measured, 0.95);
clim = [c_min, c_max];

subplot(1,3,1); imagesc(dField_measured, clim); axis image off; 
title('Measured'); colorbar; colormap(gca, 'jet');

subplot(1,3,2); imagesc(sim_map, clim); axis image off; 
title({'Simulated', patterns{best_pattern_idx, 1}}); colorbar; colormap(gca, 'jet');

subplot(1,3,3); imagesc(dField_measured - sim_map); axis image off; 
title('Residual'); colorbar; colormap(gca, 'jet');

%% --- 修正版: Mask適用 3D可視化 ---
% Mask.mat の情報を使い、背景をNaN化して脳だけを表示します

% 1. マスクの読み込みとスライス抽出
if exist(path_mat_ON_Mag, 'file')
    m_data = load(path_mat_ON_Mag);
    if isfield(m_data, 'Mask')
        % Mask変数が論理値(logical)か確認しdouble化
        raw_mask = double(m_data.Mask);
    elseif isfield(m_data, 'iMag')
        % Maskがない場合はMagから自動生成
        raw_mask = double(m_data.iMag > max(m_data.iMag(:))*0.1);
    end
    
    % データと同じスライス位置を切り出す
    if ndims(raw_mask) == 3
        % 直前の解析で sliceIdx が定義されている前提
        if ~exist('sliceIdx', 'var')
             sliceIdx = round(size(raw_mask, 3) / 2); 
        end
        display_mask = raw_mask(:,:,sliceIdx);
    else
        display_mask = raw_mask;
    end
else
    error('Maskファイルが見つかりません');
end

% 2. 表示用データの作成 (背景をNaNにする)
% NaNを入れると、MATLABのsurfはその部分を透明(描画しない)扱いにします
plot_meas = dField_measured;
plot_meas(display_mask == 0) = NaN;

plot_sim = sim_map;
plot_sim(display_mask == 0) = NaN;

plot_resid = dField_measured - sim_map;
plot_resid(display_mask == 0) = NaN;


% 3. 描画
figure('Name', 'Masked 3D Shim Field', 'Color', 'w', 'Position', [100, 100, 1200, 500]);
[h, w] = size(plot_meas);
[X, Y] = meshgrid(1:w, 1:h);

% Z軸の範囲設定 (外れ値を飛ばして見やすくする)
valid_vals = plot_meas(~isnan(plot_meas));
z_min = quantile(valid_vals, 0.01);
z_max = quantile(valid_vals, 0.99);

% --- (1) Measured ---
subplot(1, 3, 1);
surf(X, Y, plot_meas, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
title('Measured (Masked)');
view(-30, 60); % 見やすい角度
axis tight; axis off; % 軸や枠線を消してスッキリさせる
colormap(gca, 'jet');
caxis([z_min, z_max]); 
zlim([z_min, z_max]); % スパイクを除去して拡大
camlight; lighting gouraud; % 立体感

% --- (2) Simulated ---
subplot(1, 3, 2);
surf(X, Y, plot_sim, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
title({'Simulated', patterns{best_pattern_idx, 1}});
view(-30, 60);
axis tight; axis off;
colormap(gca, 'jet');
caxis([z_min, z_max]); 
zlim([z_min, z_max]);
camlight; lighting gouraud;

% --- (3) Residual ---
subplot(1, 3, 3);
surf(X, Y, plot_resid, 'EdgeColor', 'none', 'FaceAlpha', 1.0);
title('Residual (Masked)');
view(-30, 60);
axis tight; axis off;
colormap(gca, 'jet');
caxis([z_min, z_max]); 
zlim([z_min, z_max]);
camlight; lighting gouraud;

% マウス操作を有効化
rotate3d on;

%% 関数群
function [shimStruct, dicomInfo] = get_dicom_info(targetPath)
    shimStruct = struct(); dicomInfo = [];
    if ~isfolder(targetPath), error(['No folder: ' targetPath]); end
    files = dir(fullfile(targetPath, '*'));
    for i = 1:length(files)
        fname = files(i).name;
        if startsWith(fname,'.') || files(i).isdir, continue; end
        try
            info = dicominfo(fullfile(files(i).folder, fname));
            if isfield(info, 'Private_0029_1022')
                raw = info.Private_0029_1022;
                if isa(raw,'uint8')||isa(raw,'int8'), raw=char(raw'); else, raw=string(raw); end
                pts = strsplit(raw, ',');
                for k=1:length(pts)
                    it=strtrim(pts{k});
                    if contains(it,'='), kv=strsplit(it,'='); 
                        v=str2double(kv{2}); if ~isnan(v), shimStruct.(strtrim(kv{1}))=v; end
                    end
                end
                dicomInfo = info; break;
            end
        catch, continue; end
    end
end

function [m, f] = load_mat_data(fm, fp)
    d=load(fm); if isfield(d,'iMag'),m=d.iMag; elseif isfield(d,'Mask'),m=d.Mask; else, m=[]; end
    d=load(fp); if isfield(d,'iFreq'),f=d.iFreq; elseif isfield(d,'phase'),f=d.phase; else, f=[]; end
    m=double(m); f=double(f);
end