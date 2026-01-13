%================================================================================================
% QSM解析実行スクリプト (Read_Raw_Data.m 使用例)
%================================================================================================
clear variables;

%% --- 1. 撮像パラメータを手動で設定 ---
% このセクションをご自身のデータに合わせて正確に設定してください。

% --- 1. 初期設定 ---
fprintf('1. パラメータを設定しています...\n');

% パス設定
image_file_dual_echo = 'F:\hamaguchi\data\20251215\dual_echo\27Z'; % !! 要変更 !!

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
% save_path   = "F:\hamaguchi\20251204\2DSE\3_output_data";
%% 1. 設定（ここだけ書き換えてください）
% DICOMが入っているフォルダのパス
image_file_dir = fullfile(image_file_0, image_file_1); 
% image_file_dir = 'F:\hamaguchi\20251204\2DSE\2_original_data';
%% 2. データの読み込み
fprintf('データを読み込んでいます...\n');
% メーカー指定オプションを入れて読み込みます
[iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, files] = ...
    Read_DICOM(image_file_dir, 'manufacturer', 'Hitachi Medical Corporation');

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
save(fullfile(save_path, 'B0_Map_Result.mat'), 'iFreq', 'iFreq_raw', 'Mask', 'voxel_size', 'matrix_size', 'CF');
fprintf('結果を B0_Map_Result.mat に保存しました。\n');