%% ========================================================================
%  QSM Ultra Strong Smoothing Script
%  これは、位相を保存するためだけのプログラム
% ========================================================================
clear variables; close all; clc;

% --- 設定 ---
base_dir = 'C:\Users\yasun\Documents\b0_mapping_project\data\20251215\dual_echo'; 

base_dir ='C:\Users\hamaguchi\project\b0_mapping_project\data\20251215\dual_echo';
toolbox_root_path = 'C:\Users\hamaguchi\project\mlb-MRI\MEDI_toolbox-2024.11.26'; 

if exist(toolbox_root_path, 'dir')
    addpath(genpath(toolbox_root_path));
    fprintf('Toolboxパスを追加しました: %s\n', toolbox_root_path);
else
    % パスが見つからない場合はエラーで止める
    error('Toolboxのフォルダが見つかりません。パスを確認してください:\n%s', toolbox_root_path);
end
target_ids = 24:33; 

fprintf('【超】強力スムージング処理を開始します (N_std=0.1, Lambda=5000)...\n');

for idx = target_ids
    clearvars -except idx target_ids base_dir toolbox_func_path;
    case_name = num2str(idx);
    fprintf('\nProcessing Case %s ...\n', case_name);
    
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

        % 2. 位相計算
        ComplexDiff = iField(:,:,:,2) .* conj(iField(:,:,:,1));
       
        % 変更後：変数を一度保存する
        WrappedPhase = angle(ComplexDiff); % これがアンラップ前の位相 (-pi から +pi)
        UnwrappedPhase = unwrapPhase(iMag, WrappedPhase, matrix_size);
        
        % 表示してみる（スライス56などを指定）
        figure; imagesc(rot90(WrappedPhase(:,:,56))); axis image off; colormap gray;
        title('Wrapped Phase'); colorbar;

        iFreq_ppm = UnwrappedPhase / (2 * pi * delta_TE * CF/1e6);
        iFreq_ppm(isnan(iFreq_ppm)) = 0;

        % 3. 背景磁場除去 (PDF) - Robust Weighting
        fprintf('  PDF running...\n');
        safe_margin = 0.05 * mean(iMag(Mask==1)); 
        N_std_est = 1 ./ (iMag + safe_margin);
        N_std_est = N_std_est .* Mask;
        avg_noise = mean(N_std_est(Mask==1));
        clip_val = avg_noise * 20; 
        N_std_est(N_std_est > clip_val) = clip_val;
        N_std_est = N_std_est / avg_noise;

        try
            RDF_ppm = PDF(iFreq_ppm, N_std_est, Mask, matrix_size, voxel_size, B0_dir);
            if any(isnan(RDF_ppm(:))), error('NaN'); end
        catch
            RDF_ppm = (iFreq_ppm - SMV(iFreq_ppm, matrix_size, voxel_size, 5)) .* Mask;
        end

        % 4. QSM計算 (超強力設定)
        scale_factor = 2 * pi * delta_TE * CF * 1e-6;
        
        iFreq = RDF_ppm * scale_factor;
        RDF   = iFreq; 
        
        % 【設定変更 A】ノイズ見積もりを0.1まで上げる (元は0.002)
        % これによりデータの信頼度を下げ、平滑化を受け入れさせます
        N_std = (ones(size(iFreq)) * 0.1) * scale_factor; 
        
        save(fullfile(save_path, 'RDF.mat'), 'WrappedPhase', 'iFreq', 'RDF', 'N_std', 'iMag', 'Mask', ...
             'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir');
    % 
    %     current_dir = pwd;
    %     cd(save_path);
    % 
    %     fprintf('  MEDI_L1 running (Aggressive Lambda=5000)...\n');
    %     % 【設定変更 B】Lambdaを5000まで上げる (元は1000)
    %     QSM = MEDI_L1('lambda', 5000); 
    % 
    %     cd(current_dir);
    % 
    %     % 5. 保存
    %     [~, ~, nz_orig] = size(Mask_raw);
    %     QSM_final = zeros(size(Mask_raw));
    %     if size(QSM, 3) ~= nz_orig
    %         z_indices = find(squeeze(sum(sum(Mask, 1), 2)) > 0);
    %         if ~isempty(z_indices)
    %             sz = min(size(QSM,3), nz_orig - min(z_indices) + 1);
    %             QSM_final(:,:, min(z_indices):min(z_indices)+sz-1) = QSM(:,:, 1:sz);
    %         else
    %             QSM_final = QSM;
    %         end
    %     else
    %         QSM_final = QSM;
    %     end
    % 
    %     % ファイル名を UltraSmooth にして区別
    %     save(fullfile(save_path, 'QSM_UltraSmooth.mat'), 'QSM_final', 'RDF_ppm', 'Mask');
    %     fprintf('  Done! Saved as QSM_UltraSmooth.mat\n');
    % 
    catch ME
        cd(base_dir);
        fprintf('  Error in Case %s: %s\n', case_name, ME.message);
    end
end
fprintf('All Done.\n');