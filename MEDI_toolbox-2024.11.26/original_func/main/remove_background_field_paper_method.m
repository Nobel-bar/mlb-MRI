%% Background Field Removal (Based on PMC6239980)
% ファントムフィッティングとシム差分補正によるEF/SF分離
clear; clc; close all;

%% 1. 設定：解析パラメータ & ファイルパス
% -------------------------------------------------------------------------
% 前回の解析で特定した値を入力してください
SENSITIVITY = 1.1137e-04;  % 感度 (Hz/DAC)
PATTERN_IDX = 2;           % 2 = 'Swap XY' (最適パターン)

% --- 人体データ (Human / Subject) ---
path_hum_mag   = 'F:\hamaguchi\20251215\dual_echo\27\3_qsm_data\Mask.mat';
path_hum_phase = 'F:\hamaguchi\20251215\dual_echo\27\3_qsm_data\phase.mat';
dicom_hum      = 'F:\hamaguchi\20251215\dual_echo\27\1_original_data';

% --- ファントムデータ (Phantom / Reference) ---
% ※シムOFFのデータ、あるいは基準となるファントムデータ
path_pha_mag   = 'F:\hamaguchi\20251215\dual_echo\27Z\3_qsm_data\Mask.mat';
path_pha_phase = 'F:\hamaguchi\20251215\dual_echo\27Z\3_qsm_data\phase.mat';
dicom_pha      = 'F:\hamaguchi\20251215\dual_echo\27Z\1_original_data';
% -------------------------------------------------------------------------

%% 2. データの読み込み
fprintf('データを読み込んでいます...\n');
[img_mag_hum, img_freq_hum] = load_mat_data(path_hum_mag, path_hum_phase);
[img_mag_pha, img_freq_pha] = load_mat_data(path_pha_mag, path_pha_phase);

% 3Dデータ処理 (中央スライス抽出)
if ndims(img_mag_hum) == 3
    sl = round(size(img_mag_hum, 3) / 2);
    img_mag_hum = img_mag_hum(:,:,sl); img_freq_hum = img_freq_hum(:,:,sl);
    img_mag_pha = img_mag_pha(:,:,sl); img_freq_pha = img_freq_pha(:,:,sl);
end

% マスク作成
mask_hum = img_mag_hum > (max(img_mag_hum(:))*0.1);
mask_pha = img_mag_pha > (max(img_mag_pha(:))*0.1);

%% 3. Detrend (デトレンド) - 0次オフセット補正
% 論文にある "Initially acquired data is used as a baseline to detrend" の簡易実装
% 平均周波数を0に合わせることで、全体的なズレを補正します。
mean_hum = mean(img_freq_hum(mask_hum));
mean_pha = mean(img_freq_pha(mask_pha));

TF_hum = img_freq_hum - mean_hum; % Total Field (Human)
B_pha  = img_freq_pha - mean_pha; % Phantom Base Field

%% 4. Phantom Polynomial Fitting (3次多項式)
% ファントムの磁場を3次関数で近似し、滑らかな「マグネット磁場」を推定します
fprintf('ファントムデータを3次多項式でフィッティング中...\n');

[ny, nx] = size(B_pha);
dx = 1; dy = 1; 
R0 = max(nx, ny) / 2;
x_vec = ((1:nx) - nx/2) / R0;
y_vec = ((1:ny) - ny/2) / R0;
[xx, yy] = meshgrid(x_vec, y_vec);

% 座標変換 (Human解析時と同じパターンを適用)
switch PATTERN_IDX
    case 2, X=yy; Y=xx; % Swap XY
    otherwise, X=xx; Y=yy; % Default
end

% フィッティング用基底関数 (3次まで: 1, x, y, x2, xy, y2, x3, x2y, xy2, y3)
idx_pha = find(mask_pha);
x_p = X(idx_pha); y_p = Y(idx_pha);
b_p = B_pha(idx_pha);

A_poly = [ones(size(x_p)), ...
          x_p, y_p, ...
          x_p.^2, x_p.*y_p, y_p.^2, ...
          x_p.^3, x_p.^2.*y_p, x_p.*y_p.^2, y_p.^3];

coeffs_poly = A_poly \ b_p;

% 全画面分の3次曲面を生成 (B_magnet_fitted)
A_full = [ones(numel(X),1), ...
          X(:), Y(:), ...
          X(:).^2, X(:).*Y(:), Y(:).^2, ...
          X(:).^3, X(:).^2.*Y(:), X(:).*Y(:).^2, Y(:).^3];
      
B_magnet_fitted = reshape(A_full * coeffs_poly, ny, nx);

%% 5. Delta Shim Calculation (シム差分の計算)
% EF = B_magnet_fitted + Delta_Shim
fprintf('シム差分を計算中...\n');

[shim_hum, ~] = get_dicom_info(dicom_hum);
[shim_pha, ~] = get_dicom_info(dicom_pha);

dShim = struct();
tags = {'X','Y','Z','XY','XZ','YZ','X2Y2','Z2'};
for i = 1:length(tags)
    t = tags{i};
    v_h=0; if isfield(shim_hum,t), v_h=shim_hum.(t); end
    v_p=0; if isfield(shim_pha,t), v_p=shim_pha.(t); end
    dShim.(t) = v_h - v_p; % Human - Phantom
end

% シム磁場マップの生成 (感度 SENSITIVITY を使用)
scale = SENSITIVITY; 
% ※注: 3次フィッティングの座標系(X,Y)と合わせるため、同じメッシュを使用
% Z項は2Dスライスでは定数または1次勾配として振る舞うため、簡易的に扱います
z_0 = 0; 

% 基底関数 (シムコイル用)
term_X = X; term_Y = Y; term_Z = z_0;
term_XY = X.*Y; term_XZ = X.*z_0; term_YZ = Y.*z_0;
term_X2Y2 = X.^2 - Y.^2; term_Z2 = z_0^2 - 0.5*(X.^2+Y.^2);

B_delta_shim = scale * (...
    dShim.X*term_X + dShim.Y*term_Y + dShim.Z*term_Z + ...
    dShim.XY*term_XY + dShim.XZ*term_XZ + dShim.YZ*term_YZ + ...
    dShim.X2Y2*term_X2Y2 + dShim.Z2*term_Z2);

%% 6. Calculate EF & SF (分離計算)
% EF (Equipment Field) = Phantom_Fit + Delta_Shim
EF = B_magnet_fitted + B_delta_shim;

% SF (Subject Field) = TF (Human) - EF
SF = TF_hum - EF;

% マスク処理
EF(~mask_hum) = NaN;
SF(~mask_hum) = NaN;
TF_hum(~mask_hum) = NaN;

%% 7. 結果の可視化
figure('Name', 'Separation Result (PMC6239980)', 'Color', 'w', 'Position', [50, 100, 1500, 500]);

% レンジ設定 (SFを見やすくするためオートスケール)
c_sf = [quantile(SF(~isnan(SF)), 0.05), quantile(SF(~isnan(SF)), 0.95)];
c_tf = [quantile(TF_hum(~isnan(TF_hum)), 0.05), quantile(TF_hum(~isnan(TF_hum)), 0.95)];

% 1. TF (Total Field: Human Raw)
subplot(1, 4, 1);
imagesc(TF_hum, c_tf); axis image off; colormap(gca, 'jet'); colorbar;
title({'1. Total Field (TF)', 'Measured Human Data'});

% 2. EF Components (Phantom Fit)
subplot(1, 4, 2);
imagesc(B_magnet_fitted .* mask_hum, c_tf); axis image off; colormap(gca, 'jet'); colorbar;
title({'2. Phantom Fit (3rd Order)', 'Base Magnet Field'});

% 3. EF (Equipment Field: Fit + Shim)
subplot(1, 4, 3);
imagesc(EF, c_tf); axis image off; colormap(gca, 'jet'); colorbar;
title({'3. Equipment Field (EF)', 'Phantom Fit + Delta Shim'});

% 4. SF (Subject Field: TF - EF)
subplot(1, 4, 4);
imagesc(SF, c_sf); axis image off; colormap(gca, 'jet'); colorbar;
title({'4. Subject Field (SF)', 'Result: TF - EF'});

%% 8. 3D可視化 (SFの抽出確認)
figure('Name', '3D Subject Field', 'Color', 'w', 'Position', [100, 100, 1200, 500]);
[h, w] = size(SF); [mX, mY] = meshgrid(1:w, 1:h);

% EF (装置磁場)
subplot(1, 2, 1);
surf(mX, mY, EF, 'EdgeColor','none','FaceAlpha',0.8);
title('Equipment Field (EF)'); view(-30, 60); axis tight off;
colormap(gca, 'jet'); camlight; lighting gouraud;

% SF (組織磁場) -> これが欲しいデータ
subplot(1, 2, 2);
surf(mX, mY, SF, 'EdgeColor','none','FaceAlpha',1.0);
title('Subject Field (SF)'); view(-30, 60); axis tight off;
colormap(gca, 'jet'); camlight; lighting gouraud;
% SFは微細構造なのでZレンジを狭める
zlim(c_sf); 

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

これで行われていることと同じですか