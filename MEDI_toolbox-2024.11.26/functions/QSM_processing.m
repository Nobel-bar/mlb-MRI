%================================================================================================
% QSM解析実行スクリプト (Read_Raw_Data.m 使用例)
%================================================================================================
%% --- 1. 撮像パラメータを手動で設定 ---
% このセクションをご自身のデータに合わせて正確に設定してください。

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

%% --- 2. rawファイルのパスを指定 ---
% ご自身のファイルが保存されている実際のパスに書き換えてください。
mag_filepath = '/Users/nori/Downloads/matlab/data/1st_224x352_to_512x512_magnitude.raw';
phase_filepath = '/Users/nori/Downloads/matlab/data/1st_224x352_to_512x512_phase.raw'; % 位相画像のファイル名に修正してください

%% --- 3. 新しい関数でデータを読み込む ---
% Read_Raw_Data関数からすべての出力変数を受け取ります。
[iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, files] = ...
    Read_Raw_Data(mag_filepath, phase_filepath, params);

%% --- 4. 以降の処理は変更不要 ---
% この後のFit_ppm_complex, BET, unwrapPhase, PDF, MEDI_L1などの呼び出しは
% そのまま使用できます。
fprintf('データの読み込みが完了しました。磁化率マップの計算を開始します...\n');


% iFieldが4D配列であることを確認します（単一エコーでも4次元目は1になります）
% Fit_ppm_complex は4D入力を期待するため、この処理はそのままです
[iFreq_raw, N_std] = Fit_ppm_complex(iField);

%... (以降のMEDIツールボックスの解析処理を続ける)

fprintf('スクリプトのこの部分までの処理が完了しました。\n');

fprintf('スクリプトのこの部分までの処理が完了しました。\n');
