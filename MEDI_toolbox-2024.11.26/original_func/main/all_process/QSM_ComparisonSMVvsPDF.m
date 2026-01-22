%% --- QSM Comparison: SMV vs PDF ---
clear variables; close all; clc;

% ==========================================
% 設定
% ==========================================
base_dir = 'C:\Users\yasun\Documents\b0_mapping_project\data\20251215\dual_echo'; 
target_id = 28; % 比較したい症例番号

% ==========================================
case_name = num2str(target_id);
qsm_dir = fullfile(base_dir, case_name, '3_qsm_data');

path_smv = fullfile(qsm_dir, 'QSM_PDF_Final.mat');     % 以前の手法 (SMV想定)
path_pdf = fullfile(qsm_dir, 'QSM_PDF_Smooth.mat'); % 今回の手法 (PDF)

% --- 1. データ読み込み ---
if ~exist(path_smv, 'file') || ~exist(path_pdf, 'file')
    error('比較するファイルが揃っていません。\n QSM_PDF_Final: %s\n QSM_PDF_Smooth: %s', path_smv, path_pdf);
end

% SMVデータのロード
data_smv = load(path_smv, 'QSM_final', 'RDF_ppm', 'Mask');
QSM_smv = data_smv.QSM_final;
RDF_smv = data_smv.RDF_ppm;
Mask    = data_smv.Mask; % マスクは共通と仮定

% PDFデータのロード
data_pdf = load(path_pdf, 'QSM_final', 'RDF_ppm');
QSM_pdf = data_pdf.QSM_final;
RDF_pdf = data_pdf.RDF_ppm;

% 差分計算 (PDF - SMV)
Diff_QSM = QSM_pdf - QSM_smv;

% --- 2. 表示設定 ---
% 脳の中心スライスを自動特定
z_indices = find(squeeze(sum(sum(Mask, 1), 2)) > 0);
if ~isempty(z_indices)
    sl = z_indices(round(length(z_indices)/2));
else
    sl = round(size(QSM_smv, 3) / 2);
end

% 表示レンジ設定
caxis_rdf = [-0.1, 0.1];
caxis_qsm = [-0.15, 0.15];
caxis_diff = [-0.05, 0.05]; % 差分は強調して表示

% --- 3. 比較表示 (QSM & Difference) ---
fig = figure('Name', ['Comparison: Case ' case_name], 'Color', 'w', 'Position', [50, 50, 1400, 900]);

% --- Row 1: QSM画像 ---
ax1 = subplot(2, 3, 1);
imagesc(rot90(QSM_smv(:,:,sl), 1));
axis image off; colormap gray; clim(caxis_qsm);
title('Method A: QSM_PDF_Final (Previous)');

ax2 = subplot(2, 3, 2);
imagesc(rot90(QSM_pdf(:,:,sl), 1));
axis image off; colormap gray; clim(caxis_qsm);
title('Method B: QSM_PDF_Smooth (New)');

ax3 = subplot(2, 3, 3);
imagesc(rot90(Diff_QSM(:,:,sl), 1));
axis image off; colormap(ax3, jet); clim(caxis_diff); % 差分はカラー(jet)で表示
title('Difference (QSM_PDF_Final - QSM_PDF_Smooth)');
colorbar;

% --- Row 2: MinIP (静脈描出能の比較) ---
% 厚みを持たせて投影
slab = max(1, sl-5):min(size(QSM_smv,3), sl+5);
MinIP_smv = min(QSM_smv(:,:,slab), [], 3);
MinIP_pdf = min(QSM_pdf(:,:,slab), [], 3);

ax4 = subplot(2, 3, 4);
imagesc(rot90(MinIP_smv, 1));
axis image off; colormap gray; clim([-0.2, 0.1]);
title('MinIP: QSM_PDF_Final');

ax5 = subplot(2, 3, 5);
imagesc(rot90(MinIP_pdf, 1));
axis image off; colormap gray; clim([-0.2, 0.1]);
title('MinIP: QSM_PDF_Smooth');

% --- ヒストグラム比較 (値の分布) ---
subplot(2, 3, 6);
edges = linspace(-0.2, 0.2, 100);
histogram(QSM_smv(Mask==1), edges, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'DisplayName', 'SMV');
hold on;
histogram(QSM_pdf(Mask==1), edges, 'FaceColor', 'r', 'FaceAlpha', 0.5, 'DisplayName', 'PDF');
hold off;
legend; grid on;
title('Histogram Comparison');
xlabel('Susceptibility (ppm)'); ylabel('Count');
xlim([-0.2, 0.2]);

% --- ズーム同期 ---
% 片方の画像を拡大すると、他方も同時に拡大されます
linkaxes([ax1, ax2, ax3, ax4, ax5], 'xy');

fprintf('表示完了。画像をズームして細部を比較してください。\n');
fprintf('差分画像(右上)の赤/青は、手法による値のズレを示します。\n');