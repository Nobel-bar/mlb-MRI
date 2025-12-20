%================================================================================================
% QSM解析実行スクリプト (Read_Raw_Data.m 使用例)
%================================================================================================
clear all;
clear variables;

%% --- 1. 撮像パラメータを手動で設定 ---
% このセクションをご自身のデータに合わせて正確に設定してください。

% --- 1. 初期設定 ---
fprintf('1. パラメータを設定しています...\n');

% パス設定
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!

image_file_2DGE_1_2_Rotate_H_local = 'C:\Users\hamaguchi\Downloads\matlab\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_000 = "C:\Users\hamaguchi\Downloads\matlab\2DGE_0deg_H'";
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_file_temp = "C:\Users\hamaguchi\Downloads\matlab\2DGE_30sec_Rotate_H"; % !! 要変更 !!
image_file_0 = image_file_temp;
% 読み込みパスと保存パスを定義
load_base_path = fullfile(image_file_0, image_file_1);
%% --- 2. rawファイルのパスを指定 ---
% ご自身のファイルが保存されている実際のパスに書き換えてください。
mag_filepath =fullfile(load_base_path, '2DGE_30sec_Rotate_H_mag.raw');
phase_filepath =fullfile(load_base_path, '2DGE_30sec_Rotate_H_phase.raw');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_path = fullfile(image_file_0, image_file_3);

if ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf('保存フォルダを作成しました: %s\n', save_path);
end

% params構造体の初期化
params = struct();

% ボクセルサイズ [x, y, z] (mm)
params.voxel_size = [1.0, 1.0, 1.0];

% 行列サイズ [x, y, z] - ご指摘に基づき修正
params.matrix_size = [512, 512, 23];

% 中心周波数 (Hz) (例: 3Tスキャナの場合 123.2 MHz)
params.CF = 123.2 * 1e6;

% 各エコー時間 (秒単位！) - ご指摘に基づき単一エコーに修正
params.TE = [0.015];

% 静磁場の方向 [x, y, z] (通常は [0, 0, 1] または [0, 0, -1] です)
params.B0_dir = [0, 0, 1];




%% --- 3. 新しい関数でデータを読み込む ---
% Read_Raw_Data関数からすべての出力変数を受け取ります。
[iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, files] = ...
    Read_Raw_Data(mag_filepath, phase_filepath, params);

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


save(fullfile(save_path, 'phase.mat'), 'iFreq', 'iFreq_raw');

%% 4. 背景磁場除去
% 脳組織外に由来する背景磁場を除去し、局所磁場マップ(RDF)を生成します。
% ここではPDF (Projection onto Dipole Fields) を使用する例を示します。
RDF = PDF(iFreq, N_std, Mask, matrix_size, voxel_size, B0_dir);

 
save(fullfile(save_path, 'PDF.mat'), 'RDF');

%% 5. QSM再構成 (MEDI_L1)
% これがMEDIアルゴリズムの中核です。局所磁場マップ(RDF)から
% 形態情報（振幅画像）を利用して磁化率マップ(QSM)を計算します。
% MEDI+0（CSFを基準とする手法）を使用する例です。

% CSFマスクの生成（R2*マップを利用）
R2s = arlo(TE, abs(iField));
Mask_CSF = extract_CSF(R2s, Mask, voxel_size);

save(fullfile(save_path, 'other.mat'), 'N_std', 'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir', 'Mask_CSF');

path ="F:\hamaguchi\MEDI_toolbox-2024.11.26\functions";
save(fullfile(path, 'RDF.mat'),'iFreq', 'iFreq_raw', 'iMag', 'N_std', 'matrix_size', 'voxel_size', 'delta_TE', 'CF', 'B0_dir', 'Mask_CSF');

% MEDI_L1関数を呼び出し
QSM = MEDI_L1('lambda', 1000, 'lambda_CSF', 100, 'merit');

fprintf('スクリプトのこの部分までの処理が完了しました。\n');

save(fullfile(save_path, 'QSM.mat'), 'QSM');
% export QSM variable as dicom files in the 'QSM' directory
Write_DICOM(QSM, files, 'QSM')

% Initialization for Source Separation
[chi_p_init,chi_n_init,R2p,alpha,beta] = MEDI_L1ss_init(Mask,CF,R2s,QSM,delta_TE);
save RDFss.mat iFreq RDF N_std iMag Mask matrix_size voxel_size delta_TE CF B0_dir alpha beta R2p Mask_CSF chi_p_init chi_n_init

% Source Separation        
X = MEDI_L1ss('lambda',1000, 'smv',5, 'filename', 'RDFss.mat','lambda_CSF',10);
save SourceSep.mat X
