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
save(fullfile(save_path, 'Mask.mat'), 'Mask', 'iMag');

%% 2. 磁場マップとノイズの推定
% マルチエコーの複素数データ(iField)から、位相変化をフィッティングして
% ラップされた磁場マップ(iFreq_raw)とノイズ標準偏差(N_std)を計算します。
%% 3. 位相アンラッピング
% ラップされた磁場マップの不連続性を解消します。
% When using bipolar acquisition, set unipolar to false
% 3. 【重要変更】 位相フィッティング (バイポーラとして処理)
fprintf('位相マップ計算中 (Bipolarモード)... \n');
% unipolar = false として扱います
try
    [iFreq_raw, N_std] = Fit_ppm_complex_bipolar(iField);
    iFreq = unwrapPhase(iMag, iFreq_raw*2, matrix_size)/2; % Bipolar用のアンラップ補正
    fprintf(' -> Bipolar関数で計算しました。\n');
catch
    fprintf('警告: Bipolar関数がないため、Unipolarで続行しますが、delta_TE補正を試みます。\n');
    [iFreq_raw, N_std] = Fit_ppm_complex(iField);
    iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);
end

% 4. 背景磁場除去 (PDF)
fprintf('背景磁場除去 (PDF)... \n');
RDF = PDF(iFreq, N_std, Mask, matrix_size, voxel_size, B0_dir);
 
% save(fullfile(save_path, 'PDF.mat'), 'RDF');

%% 5. QSM再構成 (MEDI_L1)
% これがMEDIアルゴリズムの中核です。局所磁場マップ(RDF)から
% 形態情報（振幅画像）を利用して磁化率マップ(QSM)を計算します。
% MEDI+0（CSFを基準とする手法）を使用する例です。
% 
% --- Mask_CSF を使わない設定に変更して実行 ---
% 'lambda_CSF' オプションを削除しました
% 
% CSFマスクの生成（R2*マップを利用）
R2s = arlo(TE, abs(iField));
Mask_CSF = extract_CSF(R2s, Mask, voxel_size);




% 5. QSM計算 (単位変換なし、まずはMeritなしで確実に表示)
fprintf('QSM計算中... \n');
% ノイズマップの固定は、Bipolar処理が正しければ不要なはずですが、
% 念のため最初は固定値で画像が出るか確認します (安全策)
N_std_fixed = ones(size(RDF)) * 0.002; 

QSM = MEDI_L1('lambda', 1000, ...
              'iFreq', RDF, ...          % 単位変換なし (そのままppmとして扱う)
              'N_std', N_std_fixed, ...  % ノイズ固定 (縞模様対策)
              'Magnitude', iMag, ...
              'Mask', Mask, ...
              'matrix_size', matrix_size, ...
              'voxel_size', voxel_size, ...
              'B0_dir', B0_dir);

% 6. サイズ復元
[nx, ny, nz_orig] = size(Mask);
if size(QSM, 3) ~= nz_orig
    QSM_final = zeros(size(Mask));
    z_indices = find(squeeze(sum(sum(Mask, 1), 2)) > 0);
    if ~isempty(z_indices)
        start_z = min(z_indices);
        end_z = min(start_z + size(QSM,3) - 1, nz_orig);
        QSM_final(:,:, start_z:end_z) = QSM(:,:, 1:(end_z-start_z+1));
    else
        QSM_final = QSM; % fallback
    end
else
    QSM_final = QSM;
end

% 7. 結果表示
figure('Name', 'Bipolar Check', 'Color', 'w', 'Position', [100, 100, 1200, 600]);
sl = round(nz_orig/2);

% Input (RDF)
subplot(1, 3, 1);
imagesc(rot90(RDF(:,:,sl), 1));
axis image off; colormap gray;
caxis([-0.1, 0.1]); % ppm単位
title('Input: RDF (No conversion)');

% Output (QSM)
subplot(1, 3, 2);
imagesc(rot90(QSM_final(:,:,sl), 1));
axis image off; colormap gray;
caxis([-0.15, 0.15]); 
title('Output: QSM (Bipolar Mode)');
colorbar;

% MinIP
subplot(1, 3, 3);
slab = max(1, sl-5):min(nz_orig, sl+5);
mip_img = min(QSM_final(:,:,slab), [], 3);
imagesc(rot90(mip_img, 1));
axis image off; colormap gray;
caxis([-0.2, 0.1]); 
title('MinIP (Veins)');

fprintf('完了。縞模様が消えて血管が見えますか？\n');