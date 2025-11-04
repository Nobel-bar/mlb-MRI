fprintf('k空間の1点が「何」を表しているかを可視化します。\n');
clear variables;
close all;

%% --- 1. パラメータ設定 ---
matrix_size = 64; % デモ用の画像サイズ

% k空間の1点に設定する値
A = 1.0; % 振幅 (Amplitude)
phase = pi / 4; % 位相 (Phase) (ラジアン単位、例: 45度)
value = A * exp(1i * phase);

% 点を置くk空間の座標 (kx, ky)
% (0, 0) 以外の任意の場所
kx = 8;
ky = 5;

%% --- 2. k空間の作成 ---
% k空間をゼロで初期化 (複素数行列)
k_space = complex(zeros(matrix_size, matrix_size));

% k空間の中心インデックスを計算
% (fftshift した後の中心)
center_row = matrix_size / 2 + 1; % ky = 0 に相当
center_col = matrix_size / 2 + 1; % kx = 0 に相当

% k空間の (kx, ky) = (8, 5) の位置に値を設定
% (行列のインデックス [row, col] は [center_row - ky, center_col + kx] となります)
fprintf('k空間の点 (kx=%d, ky=%d) に A=%.1f, phase=%.2f の値を設定します。\n', ...
        kx, ky, A, phase);
k_space(center_row - ky, center_col + kx) = value;

%% --- 3. k空間から実空間へ変換 ---
% k空間の原点(DC)を (1,1) に戻してから逆フーリエ変換
% (fftshift の逆が ifftshift)
real_image = ifft2(ifftshift(k_space));

%% --- 4. 結果の表示 ---
figure('Name', 'k空間の1点と実空間画像', 'WindowState', 'maximized');

% k空間の表示
subplot(1, 3, 1);
imshow(abs(k_space), []);
title(sprintf('k空間 (強度)\n(kx=%d, ky=%d) に1点', kx, ky));
xlabel('kx');
ylabel('ky');

% 実空間の「強度」画像
subplot(1, 3, 2);
imshow(abs(real_image), []);
title('実空間画像 (強度)');
xlabel('x');
ylabel('y');

% 実空間の「位相」画像
subplot(1, 3, 3);
imshow(angle(real_image), [-pi, pi]);
colormap(subplot(1, 3, 3), 'hsv'); % 位相は 'hsv' が見やすい
colorbar;
title('実空間画像 (位相)');
xlabel('x');
ylabel('y');

fprintf('完了。\n');