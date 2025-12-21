%% Shim Simulation: 5-Way Evaluation (精度評価付き)
clear; clc; close all;

%% 1. 設定：ファイルパス
% -------------------------------------------------------------------------
% Shim ON データ (基準)
dicomFolder_ON    = 'F:\hamaguchi\20251215\dual_echo\27\1_original_data';
path_mat_ON_Mag   = 'F:\hamaguchi\20251215\dual_echo\27\3_qsm_data\Mask.mat';
path_mat_ON_Phase = 'F:\hamaguchi\20251215\dual_echo\27\3_qsm_data\phase.mat';

% Shim OFF データ
dicomFolder_OFF    = 'F:\hamaguchi\20251215\dual_echo\27Z\1_original_data';
path_mat_OFF_Mag   = 'F:\hamaguchi\20251215\dual_echo\27Z\3_qsm_data\Mask.mat';
path_mat_OFF_Phase = 'F:\hamaguchi\20251215\dual_echo\27Z\3_qsm_data\phase.mat';
% -------------------------------------------------------------------------

%% 2. データ読み込み & 前処理
fprintf('データを読み込み中...\n');
[img_mag_ON, img_freq_ON]   = load_mat_data(path_mat_ON_Mag, path_mat_ON_Phase);
[img_mag_OFF, img_freq_OFF] = load_mat_data(path_mat_OFF_Mag, path_mat_OFF_Phase);

% 3Dデータ処理 (中央スライス抽出)
if ndims(img_mag_ON) == 3
    sliceIdx = round(size(img_mag_ON, 3) / 2);
    img_mag_ON = img_mag_ON(:,:,sliceIdx); img_freq_ON = img_freq_ON(:,:,sliceIdx);
    img_mag_OFF = img_mag_OFF(:,:,sliceIdx); img_freq_OFF = img_freq_OFF(:,:,sliceIdx);
end

% 共通マスク作成
mask = (img_mag_ON > max(img_mag_ON(:))*0.1) & (img_mag_OFF > max(img_mag_OFF(:))*0.1);

% 実測差分 (Measured Delta)
dField_measured = img_freq_ON - img_freq_OFF;
dField_measured(~mask) = NaN;

%% 3. DICOM情報とシム電流
[shim_ON,  info_ON]  = get_dicom_info(dicomFolder_ON);
[shim_OFF, info_OFF] = get_dicom_info(dicomFolder_OFF);

% スライス位置とFOV情報
try
    if isfield(info_ON, 'PixelSpacing'), dy = info_ON.PixelSpacing(1); dx = info_ON.PixelSpacing(2); else, dy=1; dx=1; end
    if isfield(info_ON, 'SliceLocation'), z_mm = info_ON.SliceLocation; else, z_mm = 0; end
catch, dx=1; dy=1; z_mm=0; end

% シム電流差分
currents = struct();
tags = {'X','Y','Z','XY','XZ','YZ','X2Y2','Z2'}; 
for i = 1:length(tags)
    t = tags{i}; v_on=0; v_off=0;
    if isfield(shim_ON,t), v_on=shim_ON.(t); end
    if isfield(shim_OFF,t), v_off=shim_OFF.(t); end
    currents.(t) = v_on - v_off;
end

%% 4. シミュレーション (Best Fit探索)
[ny, nx] = size(dField_measured);
fov_x = nx * dx; fov_y = ny * dy; R0 = max(fov_x, fov_y) / 2; 
x_vec = ((1:nx) - nx/2) * dx; y_vec = ((1:ny) - ny/2) * dy;
[xx_mm, yy_mm] = meshgrid(x_vec, y_vec);
xx_norm = xx_mm / R0; yy_norm = yy_mm / R0; zz_norm = ones(size(xx_norm)) * (z_mm / R0);

idx = find(mask);
vals_measured = dField_measured(idx);

patterns = {
    'Normal', xx_norm, yy_norm; 'Swap XY', yy_norm, xx_norm;
    'Rot 90', -yy_norm, xx_norm; 'Rot 180', -xx_norm, -yy_norm;
    'Rot 270', yy_norm, -xx_norm; 'Flip X', -xx_norm, yy_norm;
    'Flip Y', xx_norm, -yy_norm; 'Swap & Flip', -yy_norm, -xx_norm; 
};

best_R = -1; best_pattern_idx = 0; best_shape = [];
for p = 1:size(patterns, 1)
    X_map = patterns{p, 2}; Y_map = patterns{p, 3}; 
    x_r = X_map(idx); y_r = Y_map(idx); z_r = zz_norm(idx);
    
    term_X = x_r; term_Y = y_r; term_Z = z_r;
    term_XY = x_r.*y_r; term_XZ = x_r.*z_r; term_YZ = y_r.*z_r;
    term_X2Y2 = x_r.^2 - y_r.^2; term_Z2 = z_r.^2 - 0.5*(x_r.^2+y_r.^2);
    
    shape_vec = currents.X*term_X + currents.Y*term_Y + currents.Z*term_Z + ...
                currents.XY*term_XY + currents.XZ*term_XZ + currents.YZ*term_YZ + ...
                currents.X2Y2*term_X2Y2 + currents.Z2*term_Z2;
    
    if std(shape_vec)~=0
        c = corrcoef(vals_measured, shape_vec);
        if abs(c(1,2)) > best_R
            best_R = abs(c(1,2)); best_pattern_idx = p; best_shape = shape_vec;
        end
    end
end

%% 5. 5つの画像の生成
poly = polyfit(best_shape, vals_measured, 1);
alpha = poly(1); beta = poly(2);

% --- 1. Measured Delta (実測差分) ---
img_1_meas_delta = dField_measured;

% --- 2. Simulated Delta (シミュレーション) ---
sim_map = zeros(ny, nx);
sim_map(idx) = alpha * best_shape + beta;
sim_map(~mask) = NaN;
img_2_sim_delta = sim_map;

% --- 3. Predicted Shim ON (OFF画像 + シミュレーション) ---
img_3_pred_ON = img_freq_OFF + sim_map;
img_3_pred_ON(~mask) = NaN;

% --- 4. Real Shim ON (実測ON画像) ---
img_4_real_ON = img_freq_ON;
img_4_real_ON(~mask) = NaN;

% --- 5. Error (Predicted ON - Real ON) ---
% これがゼロに近いほど予測精度が高い
img_5_error = img_3_pred_ON - img_4_real_ON;

%% 6. 定量評価 (Evaluation Metrics)
% 誤差の統計量を計算
err_vals = img_5_error(idx);
meas_delta_vals = vals_measured;
sim_delta_vals = sim_map(idx);

rmse = sqrt(mean(err_vals.^2)); % Root Mean Square Error
mae  = mean(abs(err_vals));     % Mean Absolute Error
delta_range = max(meas_delta_vals) - min(meas_delta_vals);

fprintf('\n========================================\n');
fprintf('   シミュレーション精度評価レポート\n');
fprintf('========================================\n');
fprintf('1. 形状の一致度 (Correlation R): %.4f\n', best_R);
fprintf('   -> 1.0に近いほど、傾きや曲面の「形」が正しく再現されています。\n');
fprintf('\n');
fprintf('2. 予測誤差 (Error Metrics):\n');
fprintf('   - RMSE (平均的なズレ)  : %.4f Hz\n', rmse);
fprintf('   - MAE  (絶対値誤差平均): %.4f Hz\n', mae);
fprintf('\n');
fprintf('3. 評価指標 (対レンジ比):\n');
fprintf('   - 補正した総量 (Range) : %.4f Hz\n', delta_range);
fprintf('   - 残差比率 (RMSE/Range): %.2f %%\n', (rmse/delta_range)*100);
fprintf('   -> この値が小さいほど、主要な磁場変動を説明できています。\n');
fprintf('========================================\n');

%% 7. 2D可視化 (5枚並べ)
figure('Name', '5-Way Evaluation 2D', 'Color', 'w', 'Position', [50, 100, 1800, 400]);

% カラーレンジの統一
% Delta系 (1, 2)
clim_delta = [quantile(vals_measured, 0.05), quantile(vals_measured, 0.95)];
% Absolute系 (3, 4)
valid_on = img_4_real_ON(~isnan(img_4_real_ON));
clim_abs = [quantile(valid_on, 0.05), quantile(valid_on, 0.95)];
% Error系 (5) - ゼロ中心
valid_err = err_vals;
clim_err = [quantile(valid_err, 0.05), quantile(valid_err, 0.95)];

subplot(1, 5, 1); imagesc(img_1_meas_delta, clim_delta); axis image off; 
title({'1. Measured Delta', '(Shim ON - OFF)'}); colorbar; colormap(gca, 'jet');

subplot(1, 5, 2); imagesc(img_2_sim_delta, clim_delta); axis image off; 
title({'2. Simulated Delta', '(Calculated sim_map)'}); colorbar; colormap(gca, 'jet');

subplot(1, 5, 3); imagesc(img_3_pred_ON, clim_abs); axis image off; 
title({'3. Predicted Shim ON', '(Shim OFF + sim_map)'}); colorbar; colormap(gca, 'jet');

subplot(1, 5, 4); imagesc(img_4_real_ON, clim_abs); axis image off; 
title({'4. Real Shim ON', '(Actual Measurement)'}); colorbar; colormap(gca, 'jet');

subplot(1, 5, 5); imagesc(img_5_error, clim_err); axis image off; 
title({'5. Prediction Error', '(Predicted - Real)'}); colorbar; colormap(gca, 'jet');


%% 8. 3D可視化 (Mask適用・5枚並べ)
figure('Name', '5-Way Evaluation 3D', 'Color', 'w', 'Position', [50, 100, 1800, 400]);

% 背景NaN処理
p1=img_1_meas_delta; p1(~mask)=NaN;
p2=img_2_sim_delta;  p2(~mask)=NaN;
p3=img_3_pred_ON;    p3(~mask)=NaN;
p4=img_4_real_ON;    p4(~mask)=NaN;
p5=img_5_error;      p5(~mask)=NaN;

[h, w] = size(p1); [X, Y] = meshgrid(1:w, 1:h);

% 1. Measured Delta
subplot(1, 5, 1); surf(X, Y, p1, 'EdgeColor','none','FaceAlpha',1);
title('1. Measured Delta'); view(-30, 60); axis tight off;
colormap(gca,'jet'); caxis(clim_delta); zlim(clim_delta); camlight; lighting gouraud;

% 2. Simulated Delta
subplot(1, 5, 2); surf(X, Y, p2, 'EdgeColor','none','FaceAlpha',1);
title('2. Simulated Delta'); view(-30, 60); axis tight off;
colormap(gca,'jet'); caxis(clim_delta); zlim(clim_delta); camlight; lighting gouraud;

% 3. Predicted ON
subplot(1, 5, 3); surf(X, Y, p3, 'EdgeColor','none','FaceAlpha',1);
title('3. Predicted ON'); view(-30, 60); axis tight off;
colormap(gca,'jet'); caxis(clim_abs); zlim(clim_abs); camlight; lighting gouraud;

% 4. Real ON
subplot(1, 5, 4); surf(X, Y, p4, 'EdgeColor','none','FaceAlpha',1);
title('4. Real ON'); view(-30, 60); axis tight off;
colormap(gca,'jet'); caxis(clim_abs); zlim(clim_abs); camlight; lighting gouraud;

% 5. Error
subplot(1, 5, 5); surf(X, Y, p5, 'EdgeColor','none','FaceAlpha',1);
title('5. Error (Diff)'); view(-30, 60); axis tight off;
colormap(gca,'jet'); caxis(clim_err); zlim(clim_err); camlight; lighting gouraud;

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