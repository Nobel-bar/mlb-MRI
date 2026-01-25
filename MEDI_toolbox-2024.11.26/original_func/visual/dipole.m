%% ========================================================================
%  Dipole Response Model Visualization (QSM Kernel)
%  k空間上のダイポール定義と、実空間上の点広がり関数(PSF)を表示します
% ========================================================================
clear variables; close all; clc;

% --- 1. パラメータ設定 ---
N = 128;              % マトリクスサイズ (偶数推奨)
fov = 256;            % FOV (mm)
voxel_size = fov / N; % ボクセルサイズ

% グリッド生成 (-N/2 から N/2-1 まで)
[ky, kx, kz] = meshgrid(-N/2:N/2-1, -N/2:N/2-1, -N/2:N/2-1);

% --- 2. k空間 ダイポールカーネルの計算 ---
% 公式: D(k) = 1/3 - (kz^2 / k^2)
% ※ B0方向は z軸 (0,0,1) と仮定しています

k2 = kx.^2 + ky.^2 + kz.^2;
% ゼロ除算を防ぐため eps を加算、または中心を強制的に0にする
D_k = 1/3 - (kz.^2) ./ (k2 + eps);
D_k(k2 == 0) = 0; % 特異点 (DC成分) の処理

% --- 3. 実空間 ダイポール応答 (Point Spread Function) の計算 ---
% k空間のカーネルを逆フーリエ変換すると、実空間での「1点に対する磁場応答」になります
D_img = fftshift(ifftn(ifftshift(D_k)));

% --- 4. 可視化 ---
fig = figure('Name', 'Dipole Response Model', 'Color', 'w', 'Position', [100, 100, 1200, 500]);

% 中心スライスのインデックス
center_idx = N/2 + 1;

% === 左: k空間 (周波数領域) ===
subplot(1, 2, 1);
% y-z 平面 (B0方向を含む断面) を表示
imagesc(squeeze(D_k(:, center_idx, :))); 
axis image; colormap jet; colorbar;
title('1. k-space Dipole Kernel (D(k))');
xlabel('kz (B0 direction)'); ylabel('ky');
clim([-2/3, 1/3]);
% マジックアングル (約54.7度) のコーンラインを描画
hold on;
line([1, N], [1, N], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 1.5);
line([1, N], [N, 1], 'Color', 'w', 'LineStyle', '--', 'LineWidth', 1.5);
text(10, 10, 'Magic Angle Cone', 'Color', 'w', 'FontSize', 10);
hold off;

% === 右: 実空間 (画像領域) ===
subplot(1, 2, 2);
% 実空間での応答 (実部を表示)
img_slice = real(squeeze(D_img(:, center_idx, :)));

% 表示用にコントラスト強調 (中心のピークが強すぎるため対数表示または制限)
imagesc(img_slice);
axis image; colormap jet; colorbar;
title('2. Image-space Response (Dipole Field)');
xlabel('z (B0 direction)'); ylabel('y');
% コントラスト調整 (中心の極端な値を避けて表示)
c_limit = max(abs(img_slice(:))) * 0.1; 
clim([-c_limit, c_limit]);

% 解説テキスト
annotation('textbox', [0.15, 0.02, 0.7, 0.05], 'String', ...
    'B0方向(横軸)に対して、赤道方向(縦軸)はプラス、極方向(横軸)はマイナスの影響を与えます。', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 12);

fprintf('表示完了: ダイポール応答モデル\n');