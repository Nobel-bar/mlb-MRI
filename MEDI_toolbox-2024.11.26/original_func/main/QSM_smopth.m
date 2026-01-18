%% ========================================================================
%  QSM 最終修正版 (Robust Noise Weighting for PDF)
%  Target: Keep signal-weighted noise map, but prevent NaN explosion
% ========================================================================
clear variables; close all; clc;

% --- 設定 ---
base_dir = 'C:\Users\yasun\Documents\b0_mapping_project\data\20251215\dual_echo'; 
toolbox_func_path = "C:\Users\yasun\Documents\mlb-MRI\MEDI_toolbox-2024.11.26\functions"; 

if exist(toolbox_func_path, 'dir')
    addpath(genpath(toolbox_func_path));
end

% エラーが出ていた28番を含む範囲を指定
target_ids = 28:28; 

fprintf('ロバスト重み付け処理を開始します (対象: %d件)...\n', length(target_ids));

for idx = target_ids
    clearvars -except idx target_ids base_dir toolbox_func_path;
    case_name = num2str(idx);
    fprintf('\n--------------------------------------------------\n');
    fprintf('Processing Case %s ...\n', case_name);
    




    
    case_dir = fullfile(base_dir, case_name);
    input_dicom_path = fullfile(case_dir, '1_original_data');
    save_path = fullfile(case_dir, '3_qsm_data');
    
    if ~exist(input_dicom_path, 'dir'), continue; end
    if ~exist(save_path, 'dir'), mkdir(save_path); end

    try
        % 1. 読み込み & マスク (Erosion=6)
        [iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, ~] = ...
            Read_DICOM(input_dicom_path, 'manufacturer', 'Hitachi Medical Corporation');
        B0_dir = [0; 0; 1];

        iMag = sqrt(sum(abs(iField).^2, 4));
        Mask_raw = iMag > (max(iMag(:)) * 0.08); 
        for i = 1:size(Mask_raw, 3), Mask_raw(:,:,i) = imfill(Mask_raw(:,:,i), 'holes'); end
        se = strel('sphere', 6); 
        Mask = imerode(Mask_raw, se);

        % 2. 位相計算 & アンラッピング
        ComplexDiff = iField(:,:,:,2) .* conj(iField(:,:,:,1));
        UnwrappedPhase = unwrapPhase(iMag, angle(ComplexDiff), matrix_size);
        iFreq_ppm = UnwrappedPhase / (2 * pi * delta_TE * CF/1e6);
        iFreq_ppm(isnan(iFreq_ppm)) = 0;

        % 3. 背景磁場除去 (PDF) - 【改良版ノイズ推定】
        fprintf('  PDF running (Robust Weighting)...\n');
        
        % Step A: 信号の底上げ（ゼロ除算防止）
        % 信号値の平均の5%を最小値として保証する
        safe_margin = 0.05 * mean(iMag(Mask==1)); 
        N_std_est = 1 ./ (iMag + safe_margin);
        
        % Step B: マスク適用
        N_std_est = N_std_est .* Mask;
        
        % Step C: 異常値のクリッピング（上限カット）
        % マスク内平均の20倍を超えるノイズ値は、計算不安定の原因になるためカット
        avg_noise = mean(N_std_est(Mask==1));
        clip_val = avg_noise * 20; 
        N_std_est(N_std_est > clip_val) = clip_val;
        
        % Step D: 正規化
        N_std_est = N_std_est / avg_noise;

        % PDF実行
        try
            RDF_ppm = PDF(iFreq_ppm, N_std_est, Mask, matrix_size, voxel_size, B0_dir);
            
            if any(isnan(RDF_ppm(:)))
                error('PDF returned NaN.');
            end
            fprintf('  -> PDF Success.\n');
            
        catch ME
            fprintf('  !!! PDF Failed (%s). Fallback to SMV. !!!\n', ME.message);
            RDF_ppm = (iFreq_ppm - SMV(iFreq_ppm, matrix_size, voxel_size, 5)) .* Mask;
        end

        % 4. QSM計算 (スムージング強化 N_std=0.02)
        % ※ここは前回の決定通り、強力なスムージングを維持します
        scale_factor = 2 * pi * delta_TE * CF * 1e-6;
        
        iFreq = RDF_ppm * scale_factor;
        RDF   = iFreq; 
        N_std = (ones(size(iFreq)) * 0.02) * scale_factor; 
        
        save(fullfile(save_path, 'RDF.mat'), 'iFreq', 'RDF', 'N_std', 'iMag', 'Mask', ...
             'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir');
        
        current_dir = pwd;
        cd(save_path);
        
        fprintf('  MEDI_L1 running...\n');
        QSM = MEDI_L1('lambda', 1000);
        
        cd(current_dir);

        % 5. 保存
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
        
        save(fullfile(save_path, 'QSM_PDF_Smooth.mat'), 'QSM_final', 'RDF_ppm', 'Mask');
        fprintf('  Done! Saved as QSM_PDF_Smooth.mat\n');

    catch ME
        cd(base_dir);
        fprintf('  Error in Case %s: %s\n', case_name, ME.message);
    end
end
fprintf('All Done.\n');