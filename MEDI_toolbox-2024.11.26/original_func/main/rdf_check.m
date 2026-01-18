%% --- Case 27の QSM_final と RDF_ppm の確認 ---
clear variables; close all; clc;

% 設定
check_id = 34; % マスクが成功していたCase 27を指定
base_dir = 'F:\hamaguchi\data\20251215\dual_echo'; 

base_dir = 'C:\Users\yasun\Documents\b0_mapping_project\data\20251215\dual_echo'; % 親フォルダ
% パス設定
case_name = num2str(check_id);
data_path = fullfile(base_dir, case_name, '3_QSM_data', 'QSM_PDF_Final.mat');

if ~exist(data_path, 'file')
    error('データが見つかりません: %s', data_path);
end

% データ読み込み
load(data_path, 'QSM_final', 'RDF_ppm', 'Mask');

% 表示スライス位置の決定
[nx, ny, nz] = size(QSM_final);
sl = round(nz / 2); % 真ん中のスライス

% --- 描画 ---
figure('Name', ['Check Case ' case_name], 'Position', [100, 100, 1200, 500], 'Color', 'w');

% 1. Local Field (RDF_ppm) - QSM_finalの入力となる画像
subplot(1, 3, 1);
imagesc(rot90(squeeze(RDF_ppm(:,:,sl)), 1));
axis image off; colormap gray;
caxis([-0.05, 0.05]); % 局所磁場用のコントラスト
title('Input: Local Field (RDF_ppm)');
colorbar;

% 2. QSM_final - 結果画像
subplot(1, 3, 2);
imagesc(rot90(squeeze(QSM_final(:,:,sl)), 1));
axis image off; colormap gray;
caxis([-0.15, 0.15]); % QSM_final用のコントラスト
title('Output: QSM_final');
colorbar;

% 3. Mask - 範囲
subplot(1, 3, 3);
imagesc(rot90(squeeze(Mask(:,:,sl)), 1));
axis image off; colormap gray;
title('Mask Region');

fprintf('Case %d の確認画像を表示しました。\n', check_id);