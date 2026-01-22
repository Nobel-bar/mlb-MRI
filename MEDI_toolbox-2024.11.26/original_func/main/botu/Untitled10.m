%% --- QSM解析 最終完成版 (All Fixes Integrated) ---
clear variables; close all; clc;

% =========================================================================
% 1. 初期設定 (ユーザー環境に合わせてパスを変更してください)
% =========================================================================
fprintf('1. パラメータを設定しています...\n');

% データフォルダのパス (Fドライブの方を有効にしています)
base_dir = 'C:\Users\yasun\Documents\b0_mapping_project\data\20251215\dual_echo'; % !! 要変更 !!

base_dir ='C:\Users\hamaguchi\project\b0_mapping_project\data\20251215\dual_echo';
target_id = '27'; % 解析対象のフォルダ名

image_file_dual_echo = fullfile(base_dir, target_id);
image_file_1 = '1_original_data';
image_file_3 = '3_qsm_data'; 

% 保存先設定
save_path = fullfile(image_file_dual_echo, image_file_3);

% =========================================================================
% 2. データ読み込み
% =========================================================================
input_dicom_path = fullfile(image_file_dual_echo, image_file_1);
fprintf('DICOM読み込み中: %s\n', input_dicom_path);

% FUJIFILM/Hitachi対応
[iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, files] = Read_DICOM(input_dicom_path, 'manufacturer', 'Hitachi Medical Corporation');

% 【修正1】B0方向を強制的にZ方向[0 0 1]に固定 (計算安定化のため)
B0_dir = [0; 0; 1];
fprintf('B0方向を固定しました: [0 0 1]\n');

% =========================================================================
% 3. マスク生成 (BETの代わりに堅牢な手法を使用)
% =========================================================================
fprintf('脳マスク生成中 (Thresholding + Erosion)...\n');
iMag = sqrt(sum(abs(iField).^2, 4));

% 閾値処理 (最大輝度の10%以上を脳とする)
threshold_ratio = 0.10; 
Mask_raw = iMag > (max(iMag(:)) * threshold_ratio);

% 穴埋め処理
for i = 1:size(Mask_raw, 3)
    Mask_raw(:,:,i) = imfill(Mask_raw(:,:,i), 'holes'); 
end

% 【修正2】マスクの縮小 (脳表面の強い磁場乱れを除外するため数ピクセル削る)
se = strel('sphere', 3); % 3ピクセル分削る
Mask = imerode(Mask_raw, se); 

% =========================================================================
% 4. 位相計算とアンラッピング (手動計算で縞模様を回避)
% =========================================================================
fprintf('位相計算 & アンラッピング中...\n');

% 【修正3】手動位相差計算 (Fit_ppm_complexを使わない)
ComplexDiff = iField(:,:,:,2) .* conj(iField(:,:,:,1));
PhaseDiff = angle(ComplexDiff);

% アンラッピング
UnwrappedPhase = unwrapPhase(iMag, PhaseDiff, matrix_size);

% 一度 ppm 単位に変換 (背景磁場除去のため)
iFreq_ppm = UnwrappedPhase / (2 * pi * delta_TE * CF/1e6);
iFreq_ppm(isnan(iFreq_ppm)) = 0;

% =========================================================================
% 5. 背景磁場除去 (PDFの代わりにSMVを使用)
% =========================================================================
fprintf('背景磁場除去 (SMVフィルタ) 実行中...\n');

% 【修正4】SMVフィルタによる背景除去 (エッジのアーチファクトに強い)
smv_radius = 5; % 半径5mm
try
    Background_Field = SMV(iFreq_ppm, matrix_size, voxel_size, smv_radius);
    RDF_ppm = (iFreq_ppm - Background_Field) .* Mask;
catch
    fprintf('警告: SMV関数エラー。Gaussianフィルタで代用します。\n');
    h = fspecial('gaussian', [21 21], 5);
    Background_Field = imfilter(iFreq_ppm, h, 'replicate');
    RDF_ppm = (iFreq_ppm - Background_Field) .* Mask;
end

%% --- (前半の処理はそのまま...) ---

% =========================================================================
% 6. QSM計算 (MEDI_L1) - 修正版
% =========================================================================
fprintf('QSM計算準備 (RDF.mat 作成)...\n');

% 1. データをMEDI用に準備 (スケーリング & 変数名定義)
% MEDIは内部で iFreq または RDF という名前の変数を期待します
scale_factor = 2 * pi * delta_TE * CF * 1e-6;

iFreq = RDF_ppm * scale_factor;        % スケールアップした入力データ
RDF   = iFreq;                         % 念のため RDF という名前でも保存
N_std = (ones(size(iFreq)) * 0.002) * scale_factor; % ノイズマップ

% 2. 必要な変数をすべて RDF.mat に保存
% これで "File not found" エラーを回避します
rdf_mat_path = fullfile(save_path, 'RDF.mat');
save(rdf_mat_path, 'iFreq', 'RDF', 'N_std', 'iMag', 'Mask', ...
    'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir');

fprintf('RDF.mat を保存しました: %s\n', rdf_mat_path);
fprintf('MEDI_L1 を実行します...\n');

% 3. カレントディレクトリを一時的に変更して実行
% (MEDIツールボックスはカレントディレクトリのファイルを優先して読む癖があるため)
current_dir = pwd;
cd(save_path); 

try
    % 引数を指定せず実行（自動的に RDF.mat を読み込みます）
    % lambda だけ指定します
    QSM = MEDI_L1('lambda', 1000); 
    
    % もし output が ppm 単位に戻っていなければ戻す（通常は不要だが念のため）
    % QSM = QSM / scale_factor; 
catch ME
    cd(current_dir); % エラーが起きても必ず元の場所に戻る
    error('MEDI実行エラー: %s', ME.message);
end

cd(current_dir); % 元のディレクトリに戻る

% =========================================================================
% 7. サイズ復元と保存
% =========================================================================
% マスク生成時の元サイズを取得 (Mask_raw がなければ Mask から推定)
[nx, ny, nz_orig] = size(Mask); 
% もしErosion前のMask_rawがメモリに残っていればそれを使う
if exist('Mask_raw', 'var')
    [~, ~, nz_orig] = size(Mask_raw);
end

if size(QSM, 3) ~= nz_orig
    QSM_final = zeros(nx, ny, nz_orig);
    % マスクがある範囲を特定して埋め込む
    z_indices = find(squeeze(sum(sum(Mask, 1), 2)) > 0);
    if ~isempty(z_indices)
        start_z = min(z_indices);
        % サイズが合うように調整
        write_size = min(size(QSM,3), nz_orig - start_z + 1);
        QSM_final(:,:, start_z : start_z + write_size - 1) = QSM(:,:, 1:write_size);
    else
        QSM_final = QSM; % 復元できない場合はそのまま
    end
else
    QSM_final = QSM;
end

% 結果の保存
save(fullfile(save_path, 'QSM_Final.mat'), 'QSM_final', 'RDF_ppm', 'Mask');
fprintf('保存完了: %s\n', fullfile(save_path, 'QSM_Final.mat'));

% =========================================================================
% 8. 結果表示
% =========================================================================
figure('Name', 'Final QSM Result', 'Color', 'w', 'Position', [100, 100, 1200, 600]);
sl = round(nz_orig/2);

% Input RDF (ppm)
subplot(1, 3, 1);
imagesc(rot90(RDF_ppm(:,:,sl), 1)); 
axis image off; colormap gray; caxis([-0.1, 0.1]); 
title('Input: Local Field (ppm)');

% Output QSM
subplot(1, 3, 2);
imagesc(rot90(QSM_final(:,:,sl), 1)); 
axis image off; colormap gray; caxis([-0.15, 0.15]); 
title('Output: QSM (Susceptibility)'); colorbar;

% MinIP (静脈強調)
subplot(1, 3, 3);
slab = max(1, sl-5):min(nz_orig, sl+5);
mip_img = min(QSM_final(:,:,slab), [], 3);
imagesc(rot90(mip_img, 1)); 
axis image off; colormap gray; caxis([-0.2, 0.1]); 
title('MinIP (Veins)');

fprintf('全工程が完了しました。\n');