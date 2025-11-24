%% --- 1点ずつ処理 vs 一括処理 の比較検証プログラム ---
clear; close all; clc;

% 1. 小さなテスト画像の作成 (32x32)
N = 32;
% 真ん中に四角い物体がある画像
org_img = zeros(N, N);
org_img(N/4:3*N/4, N/4:3*N/4) = 1; 

% 2. 背景磁場（位相回転）の作成
% 画像の上から下にかけて位相がぐるっと回るような磁場
[X, Y] = meshgrid(1:N, 1:N);
background_phase = exp(1i * (Y / N * 4 * pi)); % 2回転分の位相変化

% 3. 元のk空間
k_org = fftshift(fftn(org_img));

%% --- 手法A: ループ処理 (Point-wise) ---
% あなたの元のコードのロジック
fprintf('手法A (ループ処理) を実行中...\n');
k_method_A = complex(zeros(N, N));

for i = 1:N
    for j = 1:N
        % 1点だけ取り出す
        temp_k = complex(zeros(N, N));
        temp_k(i, j) = k_org(i, j);
        
        % 実空間に戻して、磁場を掛ける
        temp_img = ifftn(ifftshift(temp_k));
        temp_img_affected = temp_img .* background_phase;
        
        % k空間に戻す
        temp_k_new = fftshift(fftn(temp_img_affected));
        
        % ★ここがポイント: 「元の座標の値」だけを取っている
        k_method_A(i, j) = temp_k_new(i, j);
    end
end

%% --- 手法B: 一括処理 (Batch) ---
% 私が提案したロジック
fprintf('手法B (一括処理) を実行中...\n');

% 実空間全体で磁場を掛けてから、一括でFFT
img_affected_all = org_img .* background_phase;
k_method_B = fftshift(fftn(img_affected_all));

%% --- 結果の比較 ---
% 実空間に戻して画像を確認
img_recon_A = ifftn(ifftshift(k_method_A));
img_recon_B = ifftn(ifftshift(k_method_B));

figure('Position', [100, 100, 1200, 600]);

% 1. 元画像
subplot(1, 3, 1);
imagesc(abs(org_img)); axis image; axis off; colormap gray;
title('Original Image');

% 2. 手法A (ループ) の結果
subplot(1, 3, 2);
imagesc(abs(img_recon_A)); axis image; axis off; colormap gray;
title({'Method A (Loop)', 'Artifact NOT appeared!', '(Incorrect)'});

% 3. 手法B (一括) の結果
subplot(1, 3, 3);
imagesc(abs(img_recon_B)); axis image; axis off; colormap gray;
title({'Method B (Batch)', 'Artifact appeared', '(Correct)'});

% 位相画像（磁場の影響）の表示
figure;
subplot(1,2,1); imagesc(angle(img_recon_A)); axis image; off; title('Method A Phase');
subplot(1,2,2); imagesc(angle(img_recon_B)); axis image; off; title('Method B Phase');