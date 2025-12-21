%================================================================================================
% B0実行スクリプト (Batch Processing: 24 to 33)
%================================================================================================
clear variables; close all; clc;

% --- 共通設定 (ループ外で定義) ---
base_dir = 'F:\hamaguchi\20251215\dual_echo'; % データの親フォルダ
toolbox_func_path = "F:\hamaguchi\MEDI_toolbox-2024.11.26\functions"; % MEDI関数のパス(変数名pathは予約語のため変更)

% 処理するフォルダ番号の範囲
target_ids = 24:33;

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
        
        %% 3. マスク作成 (BET)
        fprintf('脳マスクを作成しています...\n');
        % 複数のエコーの強度を合成してマスクを作ります
        iMag = sqrt(sum(abs(iField).^2, 4)); 
        Mask = BET(iMag, matrix_size, voxel_size); 

        %% 4. B0マップの計算 (Fitting & Unwrapping)
        fprintf('B0マップを計算しています...\n');

        % [ステップA] 複素データから周波数のズレ(Raw B0)を計算
        % iFreq_raw: 位相の折り返し（Wrap）を含んだ状態のマップ
        [iFreq_raw, N_std] = Fit_ppm_complex(iField);

        % [ステップB] 位相アンラップ（折り返しを除去）
        % iFreq: きれいに繋がったB0マップ (単位: ppm)
        iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);

        % (オプション) 背景磁場を除去したい場合（局所磁場マップにする場合）
        % RDF = PDF(iFreq, N_std, Mask, matrix_size, voxel_size, B0_dir);

        %% 5. 結果の表示と保存
        fprintf('完了しました。画像を表示します。\n');

        % 真ん中のスライスを表示
        slice = round(matrix_size(3) / 2);

        figure('Name', 'B0 Map check');
        subplot(1,3,1); imshow(iMag(:,:,slice), []); title('Magnitude (強度)');
        subplot(1,3,2); imshow(iFreq_raw(:,:,slice), []); title('Wrapped Phase (シマウマ)');
        subplot(1,3,3); imshow(iFreq(:,:,slice), []); title('B0 Map (Unwrapped)');
        colorbar;

        % 保存
        save(fullfile(save_path, 'B0_Map_Result.mat'), 'iFreq');
%         save(fullfile(save_path, 'B0_Map_Result.mat'), 'iFreq', 'iFreq_raw', 'Mask', 'voxel_size', 'matrix_size', 'CF');
        fprintf('結果を B0_Map_Result.mat に保存しました。\n');
        
        fprintf('処理完了: %s\n', case_name);
        
    catch ME
        % エラーが発生した場合の処理
        fprintf('エラー発生! (%s): %s\n', case_name, ME.message);
        fprintf('次のデータに進みます...\n');
    end
end

fprintf('\n全データの処理が終了しました。\n');