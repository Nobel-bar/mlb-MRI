% --- 1. セットアップ: お手本の画像 ---
N = 256; % 画像サイズ (256x256)

% P = phantom(N); % <- Image Processing Toolbox が必要なため削除

% === 修正箇所: Toolboxなしで「お手本の画像」を自作 ===
% (中心に円を配置する)
[X, Y] = meshgrid(linspace(-1, 1, N), linspace(-1, 1, N));
R = sqrt(X.^2 + Y.^2); % 中心からの距離
P = zeros(N, N);
P(R < 0.5) = 1; % 半径0.5の円
% === 修正ここまで ===

% --- 2. 理想のk空間（答え）を計算 ---
% fft2: 2次元フーリエ変換
% fftshift: k空間の中心（DC成分）を画像中心に移動
kspace_ideal = fftshift(fft2(P));

% --- 3. 位相エンコードによるk空間収集のシミュレーション ---

% 収集用の空のk空間を準備
kspace_acquired = zeros(N, N);

fprintf('k空間の収集中 (位相エンコードステップ)...\n');

% このforループが、TRごとに行われる位相エンコードステップをシミュレート
% (jが1からNまで変わるのが、G_y,jの強度をN段階で変えることに相当)
for j = 1:N
    % j番目のTRで、k空間の「j行目」のデータを収集する
    kspace_acquired(j, :) = kspace_ideal(j, :);
end

fprintf('k空間の収集完了。\n');

% --- 4. 画像の再構成 ---

% 収集したk空間データを逆フーリエ変換して画像を再構成
% ifftshift: k空間の中心を元に戻す
% ifft2: 2次元逆フーリエ変換
image_reconstructed = ifft2(ifftshift(kspace_acquired));

% --- 5. 結果の表示 ---

figure;

% 1. 元のファントム画像
subplot(1, 3, 1);
imshow(P, []);
title('Original Image (自作の円)');

% 2. 収集したk空間データ
% (対数を取って表示すると見やすい)
subplot(1, 3, 2);
imshow(log(abs(kspace_acquired) + 1), []);
title('Acquired k-space (収集した設計図)');

% 3. 再構成された画像
% (逆フーリエ変換の結果は複素数なので、絶対値を取って表示)
subplot(1, 3, 3);
imshow(abs(image_reconstructed), []);
title('Reconstructed Image (再構成画像)');
