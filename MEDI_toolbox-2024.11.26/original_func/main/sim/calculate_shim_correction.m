%% Shim Correction Calculator (シム補正値の算出)
% 現在のB0マップを打ち消すために必要なシム電流値(DAC)を計算します。
clear; clc; close all;

%% 1. 設定：解析結果の入力
% -------------------------------------------------------------------------
% 解析対象のデータ（Shim ONの状態、または最新のB0マップ）
path_mag   = 'F:\hamaguchi\20251215\dual_echo\25\3_qsm_data\Mask.mat';
path_phase = 'F:\hamaguchi\20251215\dual_echo\25\3_qsm_data\phase.mat';

% 直前の解析で得られた「推定感度 (Hz/DAC)」
SENSITIVITY = 1.1137e-04; 

% 直前の解析で得られた「最適パターン」
% 1:Normal, 2:SwapXY, 3:Rot90, 4:Rot180, 5:Rot270, 6:FlipX, 7:FlipY, 8:SwapFlip
PATTERN_IDX = 2; % 'Swap XY' が最適だったため 2 を指定
% -------------------------------------------------------------------------

%% 2. データの読み込み
fprintf('データを読み込んでいます...\n');
[img_mag, img_freq] = load_mat_data(path_mag, path_phase);

% 3Dデータの場合、中央スライスを使用
if ndims(img_mag) == 3
    sliceIdx = round(size(img_mag, 3) / 2);
    img_mag = img_mag(:,:,sliceIdx);
    img_freq = img_freq(:,:,sliceIdx);
end

% マスク作成（信号強度の高い部分のみ）
mask = img_mag > (max(img_mag(:)) * 0.1);

%% 3. 球面調和関数分解 (Decomposition)
[ny, nx] = size(img_freq);
dx = 1; dy = 1; % ピクセル単位で正規化するため1でOK

% 正規化座標 (-1 ~ +1) の作成
R0 = max(nx, ny) / 2;
x_vec = ((1:nx) - nx/2) / R0;
y_vec = ((1:ny) - ny/2) / R0;
[xx, yy] = meshgrid(x_vec, y_vec);

% 座標変換パターンの適用 (解析結果に合わせる)
switch PATTERN_IDX
    case 1, X=xx; Y=yy; name='Normal';
    case 2, X=yy; Y=xx; name='Swap XY';
    case 3, X=-yy; Y=xx; name='Rot 90';
    case 4, X=-xx; Y=-yy; name='Rot 180';
    case 5, X=yy; Y=-xx; name='Rot 270';
    case 6, X=-xx; Y=yy; name='Flip X';
    case 7, X=xx; Y=-yy; name='Flip Y';
    case 8, X=-yy; Y=-xx; name='Swap & Flip';
end
fprintf('適用パターン: %s\n', name);

% ROI抽出
idx = find(mask);
b_vals = img_freq(idx);
x_r = X(idx);
y_r = Y(idx);
z_r = zeros(size(x_r)); % 2D補正のためZ=0面と仮定

% 基底関数行列 A の作成 [N_pixels x N_coils]
% ※感度(Sensitivity)で割ることで、係数が直接「DAC値」になるように調整します
scale = 1.0 / SENSITIVITY; % Hz -> DAC変換係数

A = zeros(length(idx), 9);
A(:,1) = x_r       * scale; % X
A(:,2) = y_r       * scale; % Y
A(:,3) = z_r       * scale; % Z
A(:,4) = x_r.*y_r  * scale; % XY
A(:,5) = x_r.*z_r  * scale; % XZ
A(:,6) = y_r.*z_r  * scale; % YZ
A(:,7) = (x_r.^2 - y_r.^2) * scale; % X2-Y2
A(:,8) = (z_r.^2 - 0.5*(x_r.^2+y_r.^2)) * scale; % Z2
A(:,9) = ones(size(x_r));   % Constant (B0 offset)

%% 4. フィッティング & 補正値算出
% 最小二乗法: B_map = A * Coeffs
% Coeffs は「現在の磁場を作り出しているDAC相当量」です。
% したがって、打ち消すには「-Coeffs」を今のシム値に足せばOKです。

fprintf('最適化計算中...\n');
coeffs = A \ b_vals; 

% 結果の表示
tags = {'X','Y','Z','XY','XZ','YZ','X2Y2','Z2','CenterFreq'};
correction = -coeffs; % 逆符号を加算する

fprintf('\n========================================\n');
fprintf('   推奨シム補正値 (Add this to your Shim)\n');
fprintf('========================================\n');
fprintf('項目   |  現在の推定成分(DAC) |  補正値(Delta DAC)\n');
fprintf('----------------------------------------\n');
for i = 1:8
    fprintf('%-6s | %15.1f      | %+15.0f\n', ...
        tags{i}, coeffs(i), correction(i));
end
fprintf('----------------------------------------\n');
fprintf('※補正値(Delta DAC)を、現在のシム電流値に\n');
fprintf('  加算(Add)してください。\n');

%% 5. 予測される補正後の磁場
b_fitted = A * coeffs;
img_sim = zeros(ny, nx);
img_sim(idx) = b_fitted;
img_sim(~mask) = NaN;

img_residual = img_freq - img_sim; % 補正後の予測残留磁場

figure('Name', 'Shim Optimization Result', 'Color', 'w', 'Position', [100, 100, 1200, 400]);
subplot(1, 3, 1); imagesc(img_freq); axis image off; colorbar; title('Before (Current B0)'); colormap(gca,'jet');
subplot(1, 3, 2); imagesc(img_sim); axis image off; colorbar; title('Fitted Shim Field'); colormap(gca,'jet');
subplot(1, 3, 3); imagesc(img_residual); axis image off; colorbar; title('Predicted After Correction'); colormap(gca,'jet');


%% 関数
function [m, f] = load_mat_data(fm, fp)
    d=load(fm); if isfield(d,'iMag'),m=d.iMag; elseif isfield(d,'Mask'),m=d.Mask; else, m=[]; end
    d=load(fp); if isfield(d,'iFreq'),f=d.iFreq; elseif isfield(d,'phase'),f=d.phase; else, f=[]; end
    m=double(m); f=double(f);
end