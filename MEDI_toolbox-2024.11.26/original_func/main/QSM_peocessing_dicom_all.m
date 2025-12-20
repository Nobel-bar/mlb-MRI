%================================================================================================
% QSM解析実行スクリプト (Batch Processing: 24 to 33)
%================================================================================================
clear variables; close all; clc;

% --- 共通設定 (ループ外で定義) ---
base_dir = 'F:\hamaguchi\20251215\dual_echo'; % データの親フォルダ
toolbox_func_path = "F:\hamaguchi\MEDI_toolbox-2024.11.26\functions"; % MEDI関数のパス(変数名pathは予約語のため変更)

% 処理するフォルダ番号の範囲
target_ids = 33:33;

fprintf('バッチ処理を開始します。対象: %d フォルダ\n', length(target_ids));

%% --- ループ開始 ---
for idx = target_ids
    
    % ループ内で使用する変数以外をクリア（メモリ節約のため）
    % ただし、ループ制御変数(idx, target_ids)や共通設定(base_dir, toolbox_func_path)は残す
    clearvars -except idx target_ids base_dir toolbox_func_path;
    
    % フォルダ名の生成 ('24', '25' ...)
    case_name = num2str(idx);
    
    % --- 1. 撮像パラメータ・パスの設定 ---
    fprintf('\n==================================================\n');
    fprintf('処理中: %s\n', case_name);
    fprintf('==================================================\n');
    
    % 個別のパス設定
    image_file_dual_echo = fullfile(base_dir, case_name);
    
    % フォルダが存在するか確認
    if ~exist(image_file_dual_echo, 'dir')
        fprintf('警告: フォルダが見つかりません。スキップします: %s\n', image_file_dual_echo);
        continue;
    end
    
    % エラーが起きても次のループに進めるように try-catch を使用
    try
        image_file_1 = '1_original_data';
        image_file_2 = '2_data';
        image_file_3 = '3_qsm_data'; 
        image_file_4 = '4_rolate_output_data'; 
        image_file_5 = '5_fitting_output_data'; 
        
        image_file_0 = image_file_dual_echo;
        save_path = fullfile(image_file_0, image_file_3);
        
        % 保存用フォルダがなければ作成
        if ~exist(save_path, 'dir')
             mkdir(save_path);
             fprintf('保存フォルダを作成しました: %s\n', save_path);
        end

        % パラメータ設定（voxel_size等はRead_DICOMで上書きされるため、初期値が必要な場合のみ記述）
        % params.TE = [0.0046, 0.0091]; % 必要なら保持
        
        %% --- 3. 新しい関数でデータを読み込む ---
        % Read_Raw_Data関数からすべての出力変数を受け取ります。
        
        % FUJIFILMのデータをHitachiとして読み込むように指定
        input_dicom_path = fullfile(image_file_0, image_file_1);
        [iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, files] = Read_DICOM(input_dicom_path, 'manufacturer', 'Hitachi Medical Corporation');
        
        fprintf('データの読み込み完了。解析を開始します...\n');
        
        %% . 脳マスクの生成
        iMag = sqrt(sum(abs(iField).^2, 4)); % 全エコーの振幅を合成
        Mask = BET(iMag, matrix_size, voxel_size); 
        save(fullfile(save_path, 'Mask.mat'), 'Mask', 'iMag');
        
        %% 2. 磁場マップとノイズの推定 & 3. 位相アンラッピング
        % When using bipolar acquisition, set unipolar to false
        unipolar = true;
        
        if unipolar
            [iFreq_raw, N_std] = Fit_ppm_complex(iField);
            iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);
        else
            [iFreq_raw, N_std] = Fit_ppm_complex_bipolar(iField);
            iFreq = unwrapPhase(iMag, iFreq_raw*2, matrix_size)/2;     
        end
        
        save(fullfile(save_path, 'phase.mat'), 'iFreq', 'iFreq_raw');
        
        %% 4. 背景磁場除去
        RDF = PDF(iFreq, N_std, Mask, matrix_size, voxel_size, B0_dir);
        save(fullfile(save_path, 'PDF.mat'), 'RDF');
        
%         %% 5. QSM再構成 (MEDI_L1)
%         % CSFマスクの生成（R2*マップを利用）
%         R2s = arlo(TE, abs(iField));
%         
%         % 修正したextract_CSFを使用 (前の会話で修正したものをご利用ください)
%         Mask_CSF = extract_CSF(R2s, Mask, voxel_size);
%         
%         save(fullfile(save_path, 'other.mat'), 'N_std', 'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir', 'Mask_CSF');
%         
%         % MEDIツールボックスの関数パスにRDF.matを保存 (上書きになりますが、ループ内なので順次処理されます)
%         save(fullfile(toolbox_func_path, 'RDF.mat'),'iFreq', 'iFreq_raw', 'iMag', 'N_std', 'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir', 'Mask_CSF');
%         
%         % MEDI_L1関数を呼び出し
%         QSM = MEDI_L1('lambda', 1000, 'lambda_CSF', 100, 'merit');
%         
%         fprintf('QSM計算完了 (%s)\n', case_name);
%         
%         save(fullfile(save_path, 'QSM.mat'), 'QSM');
%         
%         % DICOM書き出し
%         % 出力先を明確にするため fullfile でパスを指定します
%         output_dicom_dir = fullfile(save_path, 'QSM_DICOM');
%         if ~exist(output_dicom_dir, 'dir'), mkdir(output_dicom_dir); end
%         Write_DICOM(QSM, files, output_dicom_dir);
%         
%         %% Source Separation (QSMss)
%         % Initialization
%         [chi_p_init,chi_n_init,R2p,alpha,beta] = MEDI_L1ss_init(Mask,CF,R2s,QSM,delta_TE);
%         
%         % 作業用ファイルも各患者フォルダに保存するように変更
%         temp_rdfss_path = fullfile(save_path, 'RDFss.mat');
%         save(temp_rdfss_path, 'iFreq', 'RDF', 'N_std', 'iMag', 'Mask', 'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir', 'alpha', 'beta', 'R2p', 'Mask_CSF', 'chi_p_init', 'chi_n_init');
%         
%         % Source Separation
%         % filename引数にはフルパスを渡すか、カレントディレクトリを意識する必要があります
%         X = MEDI_L1ss('lambda',1000, 'smv',5, 'filename', temp_rdfss_path, 'lambda_CSF',10);
%         
%         save(fullfile(save_path, 'SourceSep.mat'), 'X');
        
        fprintf('処理完了: %s\n', case_name);
        
    catch ME
        % エラーが発生した場合の処理
        fprintf('エラー発生! (%s): %s\n', case_name, ME.message);
        fprintf('次のデータに進みます...\n');
    end
end

fprintf('\n全データの処理が終了しました。\n');