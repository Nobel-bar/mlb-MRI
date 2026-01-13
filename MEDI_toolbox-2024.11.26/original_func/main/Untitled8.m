%================================================================================================
% QSM解析実行スクリプト (Read_Raw_Data.m 使用例)
%================================================================================================
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
unipolar = true;
if unipolar
    % 要: ユニポーラ（またはモノポーラ）読み出しは、全エコーで同じ極性（向き）の傾斜磁場を使います。これにより、位相の振る舞いは比較的単純になります。
    [iFreq_raw, N_std] = Fit_ppm_complex(iField);
    % Estimate the frequency offset in each of the voxel using a complex
    % fitting (uneven echo spacing) フィッティングの残差から計算されたノイズの標準偏差です。
    % [iFreq_raw N_std] = Fit_ppm_complex_TE(iField,TE);
    
    % Spatial phase unwrapping (region-growing)
    iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);
else
    % 概要: バイポーラ読み出しは、エコーごとに傾斜磁場の極性を反転させます（例: +G, -G, +G, ...）。これは撮像を高速化できますが、傾斜磁場の不完全性などにより、奇数番目のエコーと
    [iFreq_raw, N_std] = Fit_ppm_complex_bipolar(iField);
    % Spatial phase unwrapping (region-growing)
    % 偶数番目のエコーの間で余計な位相オフセットが生じるという問題があります。
    % wraps occur at pi, not 2pi. 
    iFreq = unwrapPhase(iMag, iFreq_raw*2, matrix_size)/2;    
end


% save(fullfile(save_path, 'phase.mat'), 'iFreq', 'iFreq_raw');

%% 4. 背景磁場除去
% 脳組織外に由来する背景磁場を除去し、局所磁場マップ(RDF)を生成します。
% ここではPDF (Projection onto Dipole Fields) を使用する例を示します。
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

% 1. 単位変換 (Hz -> ppm)
% 30Hz程度ある値を0.5ppm程度に直します
if max(abs(RDF(:))) > 1.0
    fprintf('RDFをHzからppmに変換します (CF: %.2f Hz)\n', CF);
    RDF_ppm = RDF / (CF / 1e6);
else
    RDF_ppm = RDF;
end

% 2. 【最重要】ノイズマップの強制固定
% 計算されたN_stdを使わず、経験的な固定値 (0.002 ppm) を使用します。
% これにより「正則化が効きすぎてゼロになる」問題を回避します。
fprintf('ノイズマップを固定値 (0.002) に置換します...\n');
N_std_fixed = ones(size(RDF_ppm)) * 0.002; 

% 3. QSM計算 (MEDI_L1)
fprintf('QSM再計算中 (設定: lambda=1000, meritなし)...\n');
QSM_cropped = MEDI_L1('lambda', 1000, ...
                      'iFreq', RDF_ppm, ...
                      'N_std', N_std_fixed, ...     % 固定ノイズマップ
                      'Magnitude', iMag, ...
                      'Mask', Mask, ...
                      'matrix_size', matrix_size, ...
                      'voxel_size', voxel_size, ...
                      'B0_dir', B0_dir);

% 4. サイズ復元 (マスク範囲外の埋め戻し)
[nx, ny, nz_orig] = size(Mask);
[nx_q, ny_q, nz_qsm] = size(QSM_cropped);

QSM_final = zeros(size(Mask));
if nz_qsm ~= nz_orig
    fprintf('スライスサイズ復元: %d -> %d\n', nz_qsm, nz_orig);
    z_indices = find(squeeze(sum(sum(Mask, 1), 2)) > 0);
    if ~isempty(z_indices)
        start_z = min(z_indices);
        end_z = min(start_z + nz_qsm - 1, nz_orig);
        write_len = end_z - start_z + 1;
        QSM_final(:,:, start_z:end_z) = QSM_cropped(:,:, 1:write_len);
    else
        start_z = floor((nz_orig - nz_qsm)/2) + 1;
        QSM_final(:,:, start_z:start_z+nz_qsm-1) = QSM_cropped;
    end
else
    QSM_final = QSM_cropped;
end

% 5. 結果表示
figure('Name', 'Final Solution', 'Color', 'w', 'Position', [100, 100, 1200, 600]);

% スライス選択 (脳の中央)
z_indices = find(squeeze(sum(sum(Mask, 1), 2)) > 0);
sl = z_indices(round(length(z_indices)/2));

% Input (RDF)
subplot(1, 3, 1);
imagesc(rot90(RDF_ppm(:,:,sl), 1));
axis image off; colormap gray;
caxis([-0.1, 0.1]);
title('Input: RDF (ppm)');

% Output (QSM)
subplot(1, 3, 2);
img_qsm = rot90(QSM_final(:,:,sl), 1);
imagesc(img_qsm);
axis image off; colormap gray;
caxis([-0.15, 0.15]); % コントラスト調整
title('Output: QSM (Fixed N\_std)');
colorbar;

% MIP (静脈チェック)
subplot(1, 3, 3);
slab = max(1, sl-5):min(nz_orig, sl+5);
mip_img = min(QSM_final(:,:,slab), [], 3);
imagesc(rot90(mip_img, 1));
axis image off; colormap gray;
caxis([-0.2, 0.1]);
title('MinIP (Veins)');

fprintf('完了。これで脳構造と血管が見えるはずです。\n');

% バッチ処理用の修正QSMを保存する場合のコマンド
% save(fullfile(qsm_dir, 'QSM.mat'), 'QSM_final');