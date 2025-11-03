% --- 方法B: 線形畳み込み (あなたが依頼した処理) ---
fprintf('方法B (線形畳み込み) を実行します...\n');
tic;

% --- 1. デモデータ作成 ---
Image_A = phantom(64); % 組織

% 線形畳み込みでは、Image_B は「カーネル(ぼかしフィルター)」を意味します
Image_B_Kernel = fspecial('gaussian', [7 7], 2); % ぼかしフィルター

% --- 2. 処理 (実空間で「線形畳み込み」) ---
% [★これがあなたが依頼した convn です]
Result_Image_B = convn(Image_A, Image_B_Kernel, 'same');
time_B = toc;

fprintf('方法B 完了: %.4f 秒\n', time_B);

% --- 3. 表示 ---
figure('Name', '方法B: 線形畳み込み (ぼかし処理)');
subplot(1, 2, 1);
imshow(abs(Image_A), []);
title('元の画像 (Image_A)');
subplot(1, 2, 2);
imshow(abs(Result_Image_B), []);
title('線形畳み込みした画像 (convn(A, B))');