%% ========================================================================
%  QSM High Quality Batch Script (Noise Reduction)
%  Features: PDF Method + Strong Regularization (Lambda=2500)
% ========================================================================
clear variables; close all; clc;

% --- 1. 共通設定 ---
base_dir = 'C:\Users\yasun\Documents\b0_mapping_project\data\20251215\dual_echo'; 
w

% 処理対象 (まずはCase 24でテストし、良ければ 24:33 に広げてください)
target_ids = 27:27; 

fprintf('高画質化バッチ処理を開始します。Lambda=2500, 対象: %d フォルダ\n', length(target_ids));

% --- ループ開始 ---
for idx = target_ids
    
    clearvars -except idx target_ids base_dir toolbox_func_path;
    case_name = num2str(idx);
    
    fprintf('\n--------------------------------------------------\n');
    fprintf('Processing Case %s ...\n', case_name);
    
    % パス設定
    case_dir = fullfile(base_dir, case_name);
    input_dicom_path = fullfile(case_dir, '1_original_data');
    save_path = fullfile(case_dir, '3_qsm_data');
    
    if ~exist(input_dicom_path, 'dir'), continue; end
    if ~exist(save_path, 'dir'), mkdir(save_path); end

    try
        % 1. データ読み込み & B0固定
        [iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, ~] = ...
            Read_DICOM(input_dicom_path, 'manufacturer', 'Hitachi Medical Corporation');
        B0_dir = [0; 0; 1];

        % 2. マスク生成 (Erosion=6)
        iMag = sqrt(sum(abs(iField).^2, 4));
        Mask_raw = iMag > (max(iMag(:)) * 0.08); 
        for i = 1:size(Mask_raw, 3), Mask_raw(:,:,i) = imfill(Mask_raw(:,:,i), 'holes'); end
        se = strel('sphere', 6); 
        Mask = imerode(Mask_raw, se);

        % 3. 位相計算 & アンラッピング
        ComplexDiff = iField(:,:,:,2) .* conj(iField(:,:,:,1));
        UnwrappedPhase = unwrapPhase(iMag, angle(ComplexDiff), matrix_size);
        
        iFreq_ppm = UnwrappedPhase / (2 * pi * delta_TE * CF/1e6);
        iFreq_ppm(isnan(iFreq_ppm)) = 0;

        % 4. 背景磁場除去 (PDF)
        fprintf('  Background Removal (PDF)...\n');
        N_std_est = 1 ./ (iMag + eps);
        N_std_est = N_std_est .* Mask;
        N_std_est = N_std_est / mean(N_std_est(Mask>0));

        try
            RDF_ppm = PDF(iFreq_ppm, N_std_est, Mask, matrix_size, voxel_size, B0_dir);
            if any(isnan(RDF_ppm(:))), error('NaN in PDF'); end
        catch
            fprintf('  PDF failed, fallback to SMV.\n');
            RDF_ppm = (iFreq_ppm - SMV(iFreq_ppm, matrix_size, voxel_size, 5)) .* Mask;
        end

        % 5. QSM計算 (High Regularization)
        % スケーリング
        scale_factor = 2 * pi * delta_TE * CF * 1e-6;
        
        iFreq = RDF_ppm * scale_factor;
        RDF   = iFreq; 
        N_std = (ones(size(iFreq)) * 0.002) * scale_factor; 
        
        % RDF.mat 保存
        save(fullfile(save_path, 'RDF.mat'), 'iFreq', 'RDF', 'N_std', 'iMag', 'Mask', ...
             'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir');
        
        % MEDI実行
        current_dir = pwd;
        cd(save_path);
        
        fprintf('  Running MEDI_L1 (Lambda=2500)...\n');
        % 【重要変更点】Lambdaを 1000 -> 2500 に強化してノイズを抑制
        QSM = MEDI_L1('lambda', 2500);
        
        cd(current_dir);

        % 6. サイズ復元 & 保存
        [~, ~, nz_orig] = size(Mask_raw);
        QSM_final = zeros(size(Mask_raw));
        if size(QSM, 3) ~= nz_orig
            z_indices = find(squeeze(sum(sum(Mask, 1), 2)) > 0);
            if ~isempty(z_indices)
                sz = min(size(QSM,3), nz_orig - min(z_indices) + 1);
                QSM_final(:,:, min(z_indices):min(z_indices)+sz-1) = QSM(:,:, 1:sz);
            else
                QSM_final = QSM;
            end
        else
            QSM_final = QSM;
        end
        
        % 高画質版として保存
        save(fullfile(save_path, 'QSM_HQ.mat'), 'QSM_final', 'RDF_ppm', 'Mask');
        fprintf('  Done! Saved as QSM_HQ.mat\n');

    catch ME
        cd(base_dir);
        fprintf('  Error: %s\n', ME.message);
    end
end
fprintf('\nAll processes finished.\n');