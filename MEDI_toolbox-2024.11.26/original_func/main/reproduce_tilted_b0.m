%% Reproduce Tilted B0 Field
% 傾いた頭部のB0磁場を、リファレンスSFの回転とEFの再計算で再現する
clear; clc; close all;

%% 1. 設定：ファイルパス
% -------------------------------------------------------------------------
% キャリブレーションファイル (JSON)
jsonPath = 'F:\hamaguchi\20251215\dual_echo\background/analysis_result.json';

% --- A. リファレンス画像 (傾く前の真っ直ぐなデータ: Source) ---
% ※ここからSF(組織磁場)を抽出します
path_ref_mag   = 'F:\hamaguchi\20251215\dual_echo\27\3_qsm_data\Mask.mat';
path_ref_phase = 'F:\hamaguchi\20251215\dual_echo\27\3_qsm_data\phase.mat';
dicom_ref      = 'F:\hamaguchi\20251215\dual_echo\27\1_original_data';

% --- B. ターゲット画像 (約15度傾いたデータ: Target) ---
% ※この角度を検出し、このシム電流でB0を再現します
% (実際の傾いたデータのパスに書き換えてください)
path_tilt_mag   = 'F:\hamaguchi\20251215\dual_echo\29\3_qsm_data\Mask.mat'; % 仮
path_tilt_phase = 'F:\hamaguchi\20251215\dual_echo\29\3_qsm_data\phase.mat'; % 仮
dicom_tilt      = 'F:\hamaguchi\20251215\dual_echo\29\1_original_data';      % 仮

% --- C. ファントムデータ (EF計算用) ---
path_pha_mag   = 'F:\hamaguchi\20251215\dual_echo\25\3_qsm_data\Mask.mat';
path_pha_phase = 'F:\hamaguchi\20251215\dual_echo\25\3_qsm_data\phase.mat';
dicom_pha      = 'F:\hamaguchi\20251215\dual_echo\25\1_original_data';

%% 2. 準備: パラメータとデータの読み込み
fprintf('初期化中...\n');
% JSONロード
if ~exist(jsonPath, 'file'), error('JSONなし'); end
fid=fopen(jsonPath); p=jsondecode(fread(fid,'*char')'); fclose(fid);
SENSITIVITY = p.SENSITIVITY; PATTERN_IDX = p.PATTERN_IDX;
B_magnet_fitted = p.B_magnet_fitted;

% データロード (iMag, Phase, Mask を全て読み込む)
% ※ ここで3つ目の戻り値として mask を受け取ります
[mag_ref,  freq_ref,  mask_ref]  = load_and_center(path_ref_mag,  path_ref_phase);
[mag_tilt, freq_tilt, mask_tilt] = load_and_center(path_tilt_mag, path_tilt_phase);
[mag_pha,  freq_pha,  mask_pha]  = load_and_center(path_pha_mag,  path_pha_phase);

% 読み込んだマスクがある場合はそれを使用し、なければ自動生成する保険をかける
if isempty(mask_ref),  mask_ref  = mag_ref > max(mag_ref(:))*0.15; end
if isempty(mask_tilt), mask_tilt = mag_tilt > max(mag_tilt(:))*0.15; end
% ※ファントム等は別途閾値が必要な場合があるため確認
if isempty(mask_pha),  mask_pha  = mag_pha > max(mag_pha(:))*0.15; end

% マスクを論理型(logical)に変換しておく
mask_ref = logical(mask_ref);
mask_tilt = logical(mask_tilt);

%% 3. 角度検出 (Image Registration)
fprintf('画像の傾きを検出中...\n');

% 強度画像(Magnitude)を使って回転量を計算
% optimizerの設定 (剛体変換: 回転+平行移動)
[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 300;

% 変換行列の推定 (Ref -> Tilt)
tform = imregtform(mag_ref, mag_tilt, 'rigid', optimizer, metric);

% 回転角度の抽出
T = tform.T;
rotation_rad = atan2(T(2,1), T(1,1));
rotation_deg = rad2deg(rotation_rad);
fprintf('検出された回転角度: %.2f 度\n', rotation_deg);

%% 4. リファレンスSFの抽出
fprintf('リファレンスからSFを抽出中...\n');

% 4-1. リファレンスのEFを計算 (Ref電流 - Phantom電流)
[shim_ref, ~] = get_dicom_info(dicom_ref);
[shim_pha, ~] = get_dicom_info(dicom_pha);
dShim_ref = calc_delta_shim(shim_ref, shim_pha);

% シム磁場マップ作成
[ny, nx] = size(freq_ref);
[X, Y] = get_coords(nx, ny, PATTERN_IDX);
B_shim_ref = calc_shim_map(dShim_ref, X, Y, SENSITIVITY);

% EF_ref = Magnet + Shim
EF_ref = B_magnet_fitted + B_shim_ref;

% 4-2. SF_ref の抽出 (TF - EF)
% デトレンド (0次)
mean_ref = mean(freq_ref(mask_ref));
SF_ref = (freq_ref - mean_ref) - EF_ref;
SF_ref(~mask_ref) = 0; % 回転のためにNaNではなく0にする

%% 5. SFの回転 (Rotation)
fprintf('SFを回転させてターゲットに合わせ中...\n');

% 検出した変換(tform)を使ってSFを回転
% OutputViewをターゲット画像に合わせることで、位置ズレも補正
SF_tilted_sim = imwarp(SF_ref, tform, 'OutputView', imref2d(size(mag_tilt)));

% 回転後のマスクも作成 (表示用)
mask_rotated = imwarp(mask_ref, tform, 'OutputView', imref2d(size(mag_tilt)));
mask_combined = mask_rotated & mask_tilt; % 両方にある領域

%% 6. ターゲット用EFの作成 (EF Reconstruction)
fprintf('ターゲット(傾いた状態)のEFを構築中...\n');

% 6-1. ターゲットのEFを計算 (Tilt電流 - Phantom電流)
[shim_tilt, ~] = get_dicom_info(dicom_tilt);
dShim_tilt = calc_delta_shim(shim_tilt, shim_pha);

% シム磁場マップ作成
B_shim_tilt = calc_shim_map(dShim_tilt, X, Y, SENSITIVITY);

% EF_tilt = Magnet(固定) + Shim(ターゲットの電流)
EF_tilt_sim = B_magnet_fitted + B_shim_tilt;

%% 7. B0磁場の再現と合成
fprintf('B0磁場を合成中...\n');

% Simulated B0 = Rotated SF + New EF
% ※ターゲットの平均周波数オフセット(0次)を足して実測レベルに合わせる
mean_tilt = mean(freq_tilt(mask_tilt));
B0_reproduced = SF_tilted_sim + EF_tilt_sim + mean_tilt;

% マスク適用
B0_reproduced(~mask_combined) = NaN;
freq_tilt_disp = freq_tilt;
freq_tilt_disp(~mask_combined) = NaN;

%% 8. 結果の比較・可視化 (Masked & XY-Swapped)
fprintf('可視化中 (Masked & XY-Swapped)...\n');

% --- 1. 差分の計算 (転置前) ---
% ※ ここで計算することで、データ同士の整合性を保ちます
diff_map_raw = B0_reproduced - freq_tilt;

% --- 2. マスク適用 (転置前) ---
% 全てのデータに共通マスクを適用し、背景をNaNにします
meas_masked   = freq_tilt;       meas_masked(~mask_combined)   = NaN;
sf_masked     = SF_tilted_sim;   sf_masked(~mask_combined)     = NaN;
ef_masked     = EF_tilt_sim;     ef_masked(~mask_combined)     = NaN;
repro_masked  = B0_reproduced;   repro_masked(~mask_combined)  = NaN;
diff_masked   = diff_map_raw;    diff_masked(~mask_combined)   = NaN;

% --- 3. 表示用に転置 (XY入れ替え) ---
% .' を使って行列を転置します
plot_meas  = meas_masked.';
plot_sf    = sf_masked.';
plot_ef    = ef_masked.';
plot_repro = repro_masked.';
plot_diff  = diff_masked.';

% --- 4. レンジ設定 (NaNを除外して計算) ---
valid_vals = plot_meas(~isnan(plot_meas));
if isempty(valid_vals), c_range=[-1 1]; else
    c_range = [quantile(valid_vals, 0.01), quantile(valid_vals, 0.99)];
end

valid_diff = plot_diff(~isnan(plot_diff));
if isempty(valid_diff), err_range=[-0.5 0.5]; else
    err_range = [quantile(valid_diff, 0.05), quantile(valid_diff, 0.95)];
end

figure('Name', 'Reproduction of Tilted B0 (Swapped View)', 'Color', 'w', 'Position', [50, 50, 1400, 800]);

% 1. 実測 (Measured)
subplot(2, 3, 1);
imagesc(plot_meas, c_range); axis image off; colormap(gca, 'jet'); colorbar;
title({'1. Measured B0', ['Angle: ' num2str(rotation_deg, '%.1f') ' deg']});

% 2. 抽出・回転したSF
subplot(2, 3, 4);
imagesc(plot_sf, c_range/2); axis image off; colormap(gca, 'jet'); colorbar;
title('2. Rotated Subject Field (SF)');

% 3. 新しいEF (Magnet + Tilt Shim)
subplot(2, 3, 5);
imagesc(plot_ef, c_range); axis image off; colormap(gca, 'jet'); colorbar;
title('3. New Equipment Field (EF)');

% 4. 再現 (Reproduced)
subplot(2, 3, 2);
imagesc(plot_repro, c_range); axis image off; colormap(gca, 'jet'); colorbar;
title({'4. Reproduced B0', '(Rotated SF + New EF)'});

% 5. 誤差 (Difference)
subplot(2, 3, 3);
% 誤差が見やすいように、レンジを少し狭めに固定または調整します
% 必要に応じて err_range を手動で [-0.5, 0.5] などに設定しても良いです
imagesc(plot_diff, err_range); axis image off; colormap(gca, 'jet'); colorbar;
title({'5. Difference', '(Reproduced - Measured)'});

%% --- 定量評価 (Advanced Metrics) ---

% 1. NaNを除外してベクトル化 (一列に並べる)
vec_meas  = plot_meas(~isnan(plot_meas) & ~isnan(plot_repro));
vec_repro = plot_repro(~isnan(plot_meas) & ~isnan(plot_repro));

% 2. 各種指標の計算
% (1) 相関係数 (Correlation)
R_mat = corrcoef(vec_meas, vec_repro);
r_val = R_mat(1,2);

% (2) RMSE (Root Mean Square Error) - 平均的なズレ(Hz)
rmse_val = sqrt(mean((vec_repro - vec_meas).^2));

% (3) MAE (Mean Absolute Error) - 外れ値に強い誤差(Hz)
mae_val  = mean(abs(vec_repro - vec_meas));

% (4) 回帰分析 (傾きと切片) -> 感度やオフセットの検証用
p = polyfit(vec_meas, vec_repro, 1);
slope_val = p(1);     % 理想は 1.0
intercept_val = p(2); % 理想は 0.0

% 3. コンソール出力
fprintf('\n========================================\n');
fprintf('   精度評価レポート (Quantitative Metrics)\n');
fprintf('========================================\n');
fprintf('1. 形状の一致度 (Correlation R): %.4f\n', r_val);
fprintf('   (1.0に近いほど、分布の形が似ている)\n\n');

fprintf('2. 絶対値の誤差 (Absolute Error):\n');
fprintf('   - RMSE (平均誤差) : %.2f Hz\n', rmse_val);
fprintf('   - MAE  (絶対誤差) : %.2f Hz\n', mae_val);
fprintf('   (0に近いほど、Hz単位での数値が正確)\n\n');

fprintf('3. 回帰分析 (Regression):\n');
fprintf('   - Slope (傾き)    : %.4f (理想: 1.0)\n', slope_val);
fprintf('   - Bias  (切片)    : %.2f Hz (理想: 0.0)\n', intercept_val);
fprintf('========================================\n');

% 4. 散布図 (Scatter Plot) の作成
% 実測 vs 再現 の対応関係を可視化します
figure('Name', 'Accuracy Scatter Plot', 'Color', 'w');
scatter(vec_meas(1:100:end), vec_repro(1:100:end), 10, 'filled', 'MarkerFaceAlpha', 0.5);
% ※ データ点が多すぎるので 1/100 に間引いて表示しています
hold on;

% 理想線 (y=x)
min_val = min([vec_meas; vec_repro]);
max_val = max([vec_meas; vec_repro]);
plot([min_val max_val], [min_val max_val], 'r--', 'LineWidth', 2);

grid on; axis square;
xlabel('Measured B0 (Hz)');
ylabel('Reproduced B0 (Hz)');
title({'Accuracy Check: Measured vs Reproduced', sprintf('R=%.3f, RMSE=%.1f Hz', r_val, rmse_val)});
legend('Data Points', 'Ideal (y=x)', 'Location', 'best');

%% ---------------------------------------------------------
% 補助関数群
% ---------------------------------------------------------
function [mag, freq, mask] = load_and_center(fM, fP)
    % MagnitudeとMaskの読み込み
    d = load(fM);
    
    % iMag (強度画像)
    if isfield(d, 'iMag')
        m = d.iMag;
    elseif isfield(d, 'Mask') 
        % iMagがなくMaskしかない場合のフォールバック(稀)
        m = d.Mask; 
    else
        m = [];
    end
    
    % Mask (マスク画像)
    if isfield(d, 'Mask')
        msk = d.Mask;
    else
        msk = [];
    end

    % Phase/Freqの読み込み
    d = load(fP);
    if isfield(d, 'iFreq')
        f = d.iFreq;
    elseif isfield(d, 'phase')
        f = d.phase;
    else
        f = [];
    end

    % 3Dデータの場合、中心スライスを抽出
    if ndims(m) == 3
        sl = round(size(m, 3) / 2);
        if ~isempty(m), m = m(:, :, sl); end
        if ~isempty(f), f = f(:, :, sl); end
        if ~isempty(msk), msk = msk(:, :, sl); end
    end

    % double型に変換して出力
    mag = double(m);
    freq = double(f);
    mask = double(msk);
end

function [dShim] = calc_delta_shim(s1, s2)
    tags={'X','Y','Z','XY','XZ','YZ','X2Y2','Z2'};
    for i=1:length(tags)
        t=tags{i}; v1=0; v2=0;
        if isfield(s1,t), v1=s1.(t); end
        if isfield(s2,t), v2=s2.(t); end
        dShim.(t) = v1 - v2;
    end
end

function [X, Y] = get_coords(nx, ny, p_idx)
    dx=1; dy=1; R0=max(nx,ny)/2;
    xv=((1:nx)-nx/2)/R0; yv=((1:ny)-ny/2)/R0;
    [xx,yy]=meshgrid(xv,yv);
    switch p_idx
        case 2, X=yy; Y=xx; % Swap XY
        otherwise, X=xx; Y=yy;
    end
end

function B = calc_shim_map(dS, X, Y, sens)
    z=0;
    B = sens * (dS.X*X + dS.Y*Y + dS.Z*z + ...
                dS.XY*X.*Y + dS.XZ*X.*z + dS.YZ*Y.*z + ...
                dS.X2Y2*(X.^2-Y.^2) + dS.Z2*(z^2-0.5*(X.^2+Y.^2)));
end

function [shimStruct, dicomInfo] = get_dicom_info(targetPath)
    shimStruct = struct(); dicomInfo = [];
    if ~isfolder(targetPath), return; end
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