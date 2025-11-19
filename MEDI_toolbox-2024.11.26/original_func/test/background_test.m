%% 位相エンコード波（縞模様）の可視化シミュレーション
%
% 目的：
% Kスペースの行（kyの値）を変えることが、
% 実空間の画像に「どのくらい細かい縞模様を掛けているか」
% を可視化する。

clear; close all; clc;

%% 1. ファントム（元の画像）の作成
N = 128; % 画像サイズ
[X, Y] = meshgrid(linspace(-1, 1, N)); % 座標グリッド (-1から1まで)

% シェップ・ローガン ファントム（脳の模型）を使用
Phantom = phantom('Modified Shepp-Logan', N);

figure('Position', [100 100 1200 600]);

% --- 元の画像を表示 ---
subplot(2, 4, 1);
imagesc(Phantom); colormap gray; axis image;
title('1. 元の画像 (\rho(x,y))');
xlabel('x'); ylabel('y');

%% 2. Kスペースの各行に対応する「波」をシミュレーション

% Kスペースの座標（空間周波数）
% -N/2 から N/2-1 までのインデックスに対応
ky_vector = linspace(-N/2, N/2-1, N);

% --- A. ky = 0 （Kスペース中心） ---
% ky = 0 は、波の周波数がゼロ＝「波がない（一定）」状態
ky_val_A = 0;
% 波の式: exp(-i * 2π * ky * y)
PhaseWave_A = exp(-1i * 2 * pi * (ky_val_A/N) * Y); % kyをNで規格化

% 波の実部（縞模様）を表示
subplot(2, 4, 2);
imagesc(real(PhaseWave_A)); colormap gray; axis image;
title('A: 波 (ky = 0 : 中心)');
xlabel('x'); ylabel('y');

% 画像に波を掛けた結果
WeightedImage_A = Phantom .* real(PhaseWave_A);
subplot(2, 4, 3);
imagesc(WeightedImage_A); colormap gray; axis image;
title('画像 × 波 A');
xlabel('x'); ylabel('y');

% --- B. ky = 16 （中間） ---
ky_val_B = 16;
PhaseWave_B = exp(-1i * 2 * pi * (ky_val_B/N) * Y);

% 波の実部（縞模様）を表示
subplot(2, 4, 5);
imagesc(real(PhaseWave_B)); colormap gray; axis image;
title('B: 波 (ky = 16 : 中間)');
xlabel('x'); ylabel('y');

% 画像に波を掛けた結果
WeightedImage_B = Phantom .* real(PhaseWave_B);
subplot(2, 4, 6);
imagesc(WeightedImage_B); colormap gray; axis image;
title('画像 × 波 B');
xlabel('x'); ylabel('y');

% --- C. ky = 64 （外側, 最大周波数付近） ---
ky_val_C = 64; % N/2
PhaseWave_C = exp(-1i * 2 * pi * (ky_val_C/N) * Y);

% 波の実部（縞模様）を表示
subplot(2, 4, 7);
imagesc(real(PhaseWave_C)); colormap gray; axis image;
title('C: 波 (ky = 64 : 外側)');
xlabel('x'); ylabel('y');

% 画像に波を掛けた結果
WeightedImage_C = Phantom .* real(PhaseWave_C);
subplot(2, 4, 8);
imagesc(WeightedImage_C); colormap gray; axis image;
title('画像 × 波 C');
xlabel('x'); ylabel('y');

%% 補足：積分（Kスペースの1点）
% 実際にKスペースの1点に記録されるのは、
% この「画像×波」を「全部足し合わせた（積分した）」値です。
% （今回は周波数エンコード(x)を無視して、y方向のみで示しています）

ky_0_Signal = sum(Phantom(:) .* real(PhaseWave_A(:)));
ky_16_Signal = sum(Phantom(:) .* real(PhaseWave_B(:)));
ky_64_Signal = sum(Phantom(:) .* real(PhaseWave_C(:)));

fprintf('--- Kスペースの1点の信号強度 (積分値) ---\n');
fprintf('S(ky=0) の信号強度 (マッチ度): %.1f\n', ky_0_Signal);
fprintf('S(ky=16) の信号強度 (マッチ度): %.1f\n', ky_16_Signal);
fprintf('S(ky=64) の信号強度 (マッチ度): %.1f\n', ky_64_Signal);

% 積分結果をプロット
subplot(2, 4, 4);
bar([0, 1, 2], [ky_0_Signal, ky_16_Signal, ky_64_Signal]);
set(gca, 'XTickLabel', {'ky=0', 'ky=16', 'ky=64'});
title('Kスペースの点 (積分値)');
ylabel('信号強度');