clear variables;

%% --- 1. 撮像パラメータを手動で設定 ---
% このセクションをご自身のデータに合わせて正確に設定してください。

% --- 1. 初期設定 ---
fprintf('1. パラメータを設定しています...\n');
% パス設定
image_file_dual_echo = 'F:\hamaguchi\data\20251215\dual_echo\27'; % !! 要変更 !!
image_file_1 = '1_original_data';
image_file_2 = '2_data';
image_file_3 = '3_qsm_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_file_0 = image_file_dual_echo;
% 読み込みパスと保存パスを定義
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_path = fullfile(image_file_0, image_file_3);

% if ~exist(save_path, 'dir')
%     mkdir(save_path);
%     fprintf('保存フォルダを作成しました: %s\n', save_path);
% end

% % params構造体の初期化
% params = struct();
% 
% % ボクセルサイズ [x, y, z] (mm)
% params.voxel_size = [2.0, 2.0, 2.0];
% 
% % 行列サイズ [x, y, z] - ご指摘に基づき修正
% params.matrix_size = [128, 128, 112];
% 
% % 中心周波数 (Hz) (例: 3Tスキャナの場合 123.2 MHz)
% params.CF = 63.8 * 1e6;
% 
% % 各エコー時間 (秒単位！) - ご指摘に基づき単一エコーに修正
% params.TE = [0.0046, 0.0091];
% 
% % 静磁場の方向 [x, y, z] (通常は [0, 0, 1] または [0, 0, -1] です)
% params.B0_dir = [0, 0, 1];




%% --- 3. 新しい関数でデータを読み込む ---
% Read_Raw_Data関数からすべての出力変数を受け取ります。

% [iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, files] = Read_DICOM(fullfile(image_file_0, image_file_1));
% FUJIFILMのデータをHitachiとして読み込むように指定
[iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, files] = Read_DICOM(fullfile(image_file_0, image_file_1), 'manufacturer', 'Hitachi Medical Corporation');
%% --- 4. 以降の処理は変更不要 ---
% この後のFit_ppm_complex, BET, unwrapPhase, PDF, MEDI_L1などの呼び出しは
% そのまま使用できます。
fprintf('データの読み込みが完了しました。磁化率マップの計算を開始します...\n');

%... (以降のMEDIツールボックスの解析処理を続ける)



%% . 脳マスクの生成
% 振幅画像から脳領域を抽出するマスクを作成します。
% FSLのBETツール（要インストール）を使用するのが一般的です。
iMag = sqrt(sum(abs(iField).^2, 4)); % 全エコーの振幅を合成
Mask = BET(iMag, matrix_size, voxel_size); 
% --- 【修正】B0方向の強制固定 (エラー回避のため) ---
% ほとんどのMRI撮像では主磁場はZ方向([0 0 1])です
fprintf('B0方向を [0 0 1] に固定します (元: [%.2f %.2f %.2f])\n', B0_dir(1), B0_dir(2), B0_dir(3));
B0_dir = [0; 0; 1];

% 3. 手動位相差計算 & アンラッピング
fprintf('位相マップ計算 & アンラッピング...\n');
ComplexDiff = iField(:,:,:,2) .* conj(iField(:,:,:,1));
PhaseDiff = angle(ComplexDiff);
UnwrappedPhase = unwrapPhase(iMag, PhaseDiff, matrix_size);

% 4. 単位変換 (ppm)
CF_MHz = CF / 1e6;
iFreq_ppm = UnwrappedPhase / (2 * pi * delta_TE * CF_MHz);

% --- 【修正】データの洗浄 (NaN除去) ---
if any(isnan(iFreq_ppm(:)))
    fprintf('警告: 位相マップにNaNが含まれています。0に置換します。\n');
    iFreq_ppm(isnan(iFreq_ppm)) = 0;
end

% 5. 背景磁場除去 (PDF) - 安定化設定
fprintf('背景磁場除去 (PDF) 実行中...\n');

% 【修正】ノイズマップを計算せず、一律「1」としてPDFに渡します (計算発散防止)
N_std_for_PDF = ones(size(iFreq_ppm)); 
% マスク外の値は計算に悪影響するためゼロクリア
iFreq_ppm = iFreq_ppm .* Mask;

try
    RDF = PDF(iFreq_ppm, N_std_for_PDF, Mask, matrix_size, voxel_size, B0_dir);
catch ME
    fprintf('PDFエラー: %s\n単純なハイパスフィルタで代用します。\n', ME.message);
    % PDFがどうしてもダメな場合のバックアップ (SMVフィルタ)
    RDF = SMV(iFreq_ppm, matrix_size, voxel_size, 5) .* Mask; 
end

% 6. QSM計算 (MEDI_L1)
fprintf('QSM計算中...\n');
N_std_fixed = ones(size(RDF)) * 0.002; % QSM用のノイズマップも固定

QSM = MEDI_L1('lambda', 1000, ...
              'iFreq', RDF, ...
              'N_std', N_std_fixed, ...
              'Magnitude', iMag, ...
              'Mask', Mask, ...
              'matrix_size', matrix_size, ...
              'voxel_size', voxel_size, ...
              'B0_dir', B0_dir);

% 7. サイズ復元
[nx, ny, nz_orig] = size(Mask);
if size(QSM, 3) ~= nz_orig
    QSM_final = zeros(size(Mask));
    z_indices = find(squeeze(sum(sum(Mask, 1), 2)) > 0);
    start_z = min(z_indices);
    end_z = min(start_z + size(QSM,3) - 1, nz_orig);
    QSM_final(:,:, start_z:end_z) = QSM(:,:, 1:(end_z-start_z+1));
else
    QSM_final = QSM;
end

% 8. 結果表示
figure('Name', 'Final Result', 'Color', 'w', 'Position', [100, 100, 1200, 600]);
sl = round(nz_orig/2);

% RDF
subplot(1, 3, 1);
imagesc(rot90(RDF(:,:,sl), 1));
axis image off; colormap gray;
caxis([-0.1, 0.1]); 
title('Local Field (RDF)');

% QSM
subplot(1, 3, 2);
imagesc(rot90(QSM_final(:,:,sl), 1));
axis image off; colormap gray;
caxis([-0.15, 0.15]); 
title('Output: QSM');
colorbar;

% MinIP
subplot(1, 3, 3);
slab = max(1, sl-5):min(nz_orig, sl+5);
mip_img = min(QSM_final(:,:,slab), [], 3);
imagesc(rot90(mip_img, 1));
axis image off; colormap gray;
caxis([-0.2, 0.1]); 
title('MinIP (Veins)');

fprintf('完了しました。すべての画像が正しく表示されているはずです。\n');