%% B0 Inhomogeneity Induced Distortion Visualizer
%% Reproduce Tilted B0 Fieldに続いて行わないと失敗する．
% DICOMのバンド幅情報を読み込み、B0不均一による位置ズレ(mm)を可視化します
clearvars -except B0_reproduced freq_tilt mask_tilt dicom_tilt; close all;

%% 1. 設定：データソース
% -------------------------------------------------------------------------
% 解析対象のB0マップ (Hz)
% ※ 前回の再現データがあればそれを使いますが、なければ実測を使います
if exist('B0_reproduced', 'var')
    B0_map = B0_reproduced;
    disp('Using: Reproduced B0 Map');
elseif exist('freq_tilt', 'var')
    B0_map = freq_tilt;
    disp('Using: Measured B0 Map');
else
    % データがない場合はパスを指定してロード (例)
    path_tilt_phase = 'F:\hamaguchi\20251215\dual_echo\28_tilted\3_qsm_data\phase.mat';
    d = load(path_tilt_phase); B0_map = double(d.iFreq(:,:,round(size(d.iFreq,3)/2)));
    disp('Using: Loaded B0 Map from file');
end

% DICOMフォルダ (バンド幅情報の取得元)
dicomPath = 'F:\hamaguchi\20251215\dual_echo\29\1_original_data';
% -------------------------------------------------------------------------

%% 2. DICOM情報の取得 (バンド幅・分解能)
fprintf('DICOMヘッダーを解析中...\n');
info = get_first_dicom_info(dicomPath);

% (1) Pixel Bandwidth (Hz/pixel)
if isfield(info, 'PixelBandwidth')
    bw_per_pixel = info.PixelBandwidth;
    fprintf('  - Bandwidth     : %.2f Hz/pixel\n', bw_per_pixel);
else
    bw_per_pixel = 260; % デフォルト値 (見つからない場合)
    fprintf('  - Bandwidth     : Not found. Assuming %.2f Hz/pixel\n', bw_per_pixel);
end

% (2) Pixel Spacing (mm/pixel)
if isfield(info, 'PixelSpacing')
    py = info.PixelSpacing(1); % Row spacing
    px = info.PixelSpacing(2); % Col spacing
    fprintf('  - Pixel Spacing : %.4f x %.4f mm\n', py, px);
else
    px=1; py=1;
end

% (3) Readout Direction (周波数エンコード方向)
% "ROW"ならPhaseはRow方向なので、FreqはCol方向(横)
% "COL"ならPhaseはCol方向なので、FreqはRow方向(縦)
readout_axis = 'Unknown';
if isfield(info, 'InPlanePhaseEncodingDirection')
    phase_dir = info.InPlanePhaseEncodingDirection;
    if strcmp(phase_dir, 'ROW')
        readout_axis = 'X (Columns)'; % 横方向にズレる
        pixel_res = px;
    elseif strcmp(phase_dir, 'COL')
        readout_axis = 'Y (Rows)';    % 縦方向にズレる
        pixel_res = py;
    end
    fprintf('  - Phase Encoding: %s\n', phase_dir);
    fprintf('  - Readout Axis  : %s (Distortion Direction)\n', readout_axis);
else
    % 不明な場合はX(横)と仮定
    readout_axis = 'X (Assumed)';
    pixel_res = px;
    fprintf('  - Readout Axis  : Not found. Assuming X.\n');
end

%% 3. 位置ズレ量の計算
% Shift(pixel) = B0(Hz) / BW(Hz/pixel)
shift_pixel = B0_map / bw_per_pixel;

% Shift(mm) = Shift(pixel) * Spacing(mm/pixel)
shift_mm = shift_pixel * pixel_res;

% マスク処理
if exist('mask_tilt', 'var'), shift_mm(~mask_tilt) = NaN; end

%% 4. 可視化
figure('Name', 'Distortion Map', 'Color', 'w', 'Position', [100, 100, 1200, 500]);

% レンジ設定 (99%タイル)
v = shift_mm(~isnan(shift_mm));
clim_mm = [-max(abs(v)), max(abs(v))]; % 対称レンジ

% 1. B0 Map (Hz)
subplot(1, 3, 1);
imagesc(B0_map); axis image off; colormap(gca, 'jet'); colorbar;
title('1. B0 Inhomogeneity (Hz)');

% 2. Distortion Map (mm)
subplot(1, 3, 2);
imagesc(shift_mm, clim_mm); axis image off; colormap(gca, 'jet'); h=colorbar;
ylabel(h, 'Shift (mm)');
title({'2. Positional Shift (mm)', ['Dir: ' readout_axis]});

% 3. グリッド歪みシミュレーション (Grid Distortion)
subplot(1, 3, 3);
plot_grid_distortion(shift_mm, readout_axis, 20); % 20画素ごとのグリッド
title({'3. Grid Distortion Visualization', '(Exaggerated x10 for visibility)'});
axis image off;

%% 5. 定量レポート
fprintf('\n========================================\n');
fprintf('   位置ズレ解析レポート\n');
fprintf('========================================\n');
fprintf('最大ズレ量 (Max Shift): %.4f mm\n', max(abs(v)));
fprintf('平均ズレ量 (Mean Shift): %.4f mm\n', mean(abs(v)));
fprintf('※ %.2f Hz ずれるごとに 1 pixel ずれます\n', bw_per_pixel);
fprintf('========================================\n');

%% 補助関数
function info = get_first_dicom_info(path)
    files = dir(fullfile(path, '*'));
    for i=1:length(files)
        if files(i).isdir || startsWith(files(i).name, '.'), continue; end
        try
            info = dicominfo(fullfile(path, files(i).name));
            return;
        catch, continue; end
    end
    error('DICOM file not found');
end

function plot_grid_distortion(shift_map, axis_dir, step)
    [h, w] = size(shift_map);
    [X, Y] = meshgrid(1:w, 1:h);
    
    % 歪みを適用 (視認性のため10倍に強調)
    scale_factor = 10; 
    
    % shift_map は mm 単位なので pixel に戻す必要がありますが、
    % 簡易的にここでは「マップの値」をそのままズレとして扱います
    % (厳密には pixel_shift = shift_map / pixel_res)
    
    if contains(axis_dir, 'X')
        X_dist = X + shift_map * scale_factor;
        Y_dist = Y;
    else
        X_dist = X;
        Y_dist = Y + shift_map * scale_factor;
    end
    
    % グリッド描画
    hold on;
    % 縦線
    for i = 1:step:w
        plot(X_dist(:, i), Y_dist(:, i), 'k-', 'Color', [0 0 0 0.3]);
    end
    % 横線
    for i = 1:step:h
        plot(X_dist(i, :), Y_dist(i, :), 'k-', 'Color', [0 0 0 0.3]);
    end
    set(gca, 'YDir', 'reverse');
    xlim([1 w]); ylim([1 h]);
end