%% ========================================================================
%  QSM 一括解析スクリプト (Batch Processing: 24 to 33) - PDF Version
%  Method: Manual Phase -> Unwrapping -> PDF -> Scale -> MEDI
% ========================================================================
clear variables; close all; clc;

% --- 1. 共通設定 ---
base_dir = 'C:\Users\yasun\Documents\b0_mapping_project\data\20251215\dual_echo'; % 親フォルダ
toolbox_func_path = "C:\Users\yasun\Documents\mlb-MRI\MEDI_toolbox-2024.11.26\functions"; 
target_ids = 25:34; % 処理するフォルダ範囲

fprintf('バッチ処理 (PDF版) を開始します。対象: %d フォルダ\n', length(target_ids));

% --- ループ開始 ---
for idx = target_ids
    
    % メモリ整理
    clearvars -except idx target_ids base_dir toolbox_func_path;
    case_name = num2str(idx);
    
    fprintf('\n==================================================\n');
    fprintf('処理中: Case %s\n', case_name);
    fprintf('==================================================\n');
    
    % パス設定
    case_dir = fullfile(base_dir, case_name);
    input_dicom_path = fullfile(case_dir, '1_original_data');
    save_path = fullfile(case_dir, '3_qsm_data');
    
    if ~exist(input_dicom_path, 'dir')
        fprintf('スキップ: フォルダが見つかりません (%s)\n', input_dicom_path);
        continue;
    end
    if ~exist(save_path, 'dir'), mkdir(save_path); end

    try
        % ---------------------------------------------------------
        % 1. データ読み込み
        % ---------------------------------------------------------
        [iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, ~] = ...
            Read_DICOM(input_dicom_path, 'manufacturer', 'Hitachi Medical Corporation');
        
        % B0方向を強制固定 [0 0 1]
        B0_dir = [0; 0; 1];

        % ---------------------------------------------------------
        % 2. マスク生成 (強めのErosionはPDFの安定化に必須)
        % ---------------------------------------------------------
        iMag = sqrt(sum(abs(iField).^2, 4));
        Mask_raw = iMag > (max(iMag(:)) * 0.08); 
        for i = 1:size(Mask_raw, 3), Mask_raw(:,:,i) = imfill(Mask_raw(:,:,i), 'holes'); end
        
        % 縁のアーチファクトを除去 (PDFは境界に敏感なため強めに削ります)
        se = strel('sphere', 6); 
        Mask = imerode(Mask_raw, se);

        % ---------------------------------------------------------
        % 3. 位相計算 & アンラッピング (手動)
        % ---------------------------------------------------------
        ComplexDiff = iField(:,:,:,2) .* conj(iField(:,:,:,1));
        UnwrappedPhase = unwrapPhase(iMag, angle(ComplexDiff), matrix_size);
        
        % ppm変換 (Total Field)
        iFreq_ppm = UnwrappedPhase / (2 * pi * delta_TE * CF/1e6);
        iFreq_ppm(isnan(iFreq_ppm)) = 0;

        % ---------------------------------------------------------
        % 4. 背景磁場除去 (PDFに変更)
        % ---------------------------------------------------------
        fprintf('背景磁場除去 (PDF) を実行中...\n');
        
        % PDF用のノイズマップ推定 (Magnitudeの逆数を使用)
        % 信号が低いところはノイズが大きいとみなす重み付け
        N_std_est = 1 ./ (iMag + eps);
        N_std_est = N_std_est .* Mask;
        N_std_est = N_std_est / mean(N_std_est(Mask>0)); % 正規化

        try
            % PDF実行 (入力: Total Field -> 出力: Local Field/RDF)
            RDF_ppm = PDF(iFreq_ppm, N_std_est, Mask, matrix_size, voxel_size, B0_dir);
            
            % PDFがNaNを出した場合のチェック
            if any(isnan(RDF_ppm(:)))
                error('PDF produced NaN values');
            end
            fprintf(' -> PDF成功\n');
            
        catch ME
            fprintf('警告: PDF計算に失敗しました (%s)。SMVで代用します。\n', ME.message);
            smv_radius = 5;
            Background_Field = SMV(iFreq_ppm, matrix_size, voxel_size, smv_radius);
            RDF_ppm = (iFreq_ppm - Background_Field) .* Mask;
        end

        % ---------------------------------------------------------
        % 5. QSM計算 (MEDI_L1) - スケーリング対応
        % ---------------------------------------------------------
        % スケーリング係数 (ppm -> Phase scale)
        scale_factor = 2 * pi * delta_TE * CF * 1e-6;
        
        % 変数を準備して保存 (MEDI用)
        iFreq = RDF_ppm * scale_factor;
        RDF   = iFreq; 
        N_std = (ones(size(iFreq)) * 0.002) * scale_factor; 
        
        temp_rdf_path = fullfile(save_path, 'RDF.mat');
        save(temp_rdf_path, 'iFreq', 'RDF', 'N_std', 'iMag', 'Mask', ...
             'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir');
        
        % カレントディレクトリ移動して実行
        current_dir = pwd;
        cd(save_path);
        
        % MEDI実行
        fprintf('QSM計算中...\n');
        QSM = MEDI_L1('lambda', 1000);
        
        cd(current_dir);

        % ---------------------------------------------------------
        % 6. サイズ復元と保存
        % ---------------------------------------------------------
        [~, ~, nz_orig] = size(Mask_raw);
        if size(QSM, 3) ~= nz_orig
            QSM_final = zeros(size(Mask_raw));
            z_indices = find(squeeze(sum(sum(Mask, 1), 2)) > 0);
            if ~isempty(z_indices)
                start_z = min(z_indices);
                end_z = min(start_z + size(QSM,3) - 1, nz_orig);
                QSM_final(:,:, start_z:end_z) = QSM(:,:, 1:(end_z-start_z+1));
            else
                QSM_final = QSM;
            end
        else
            QSM_final = QSM;
        end
        
        % 最終保存
        save(fullfile(save_path, 'QSM_PDF_Final.mat'), 'QSM_final', 'RDF_ppm', 'Mask');
        fprintf('完了: Case %s\n', case_name);

        % 1. Input: Local Field (RDF by PDF)
        subplot(1, 3, 1);
        imagesc(rot90(RDF_ppm(:,:,sl), 1));
        axis image off; colormap gray;
        caxis([-0.1, 0.1]); 
        title('Input: RDF (via PDF method)');
        
        % 2. Output: QSM
        subplot(1, 3, 2);
        imagesc(rot90(QSM_final(:,:,sl), 1));
        axis image off; colormap gray;
        % QSMの表示レンジ (-0.15 ～ 0.15 ppm が一般的)
        caxis([-0.15, 0.15]); 
        title('Output: QSM');
        h = colorbar; h.Label.String = 'ppm';
        
        % 3. Check: MinIP (Veins)
        % スラブ厚: 10mm程度 (前後5スライス)
        slab_start = max(1, sl - 5);
        slab_end = min(size(QSM_final,3), sl + 5);
        mip_img = min(QSM_final(:,:,slab_start:slab_end), [], 3);
        
        subplot(1, 3, 3);
        imagesc(rot90(mip_img, 1));
        axis image off; colormap gray;
        caxis([-0.2, 0.1]); % 静脈を見やすくするため少し暗めに
        title('MinIP (Veins Check)');
        
        fprintf('表示完了 (Slice %d)\n', sl);

    catch ME
        cd(base_dir); % エラー時も戻る
        fprintf('エラー発生 (Case %s): %s\n', case_name, ME.message);
    end
end

fprintf('\n全データの処理が完了しました！\n');