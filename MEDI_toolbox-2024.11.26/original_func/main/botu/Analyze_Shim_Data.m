%% 1. 初期設定とパスの指定
clear; clc; close all;

% パス設定
image_file_dual_echo = 'F:\hamaguchi\20251215\dual_echo\27Z'; % !! 要変更 !!
image_file_1 = '1_original_data';
image_file_2 = '2_data';
image_file_3 = '3_qsm_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 
path_Phantom_ShimON  = 'F:\hamaguchi\20251215\dual_echo\25\1_original_data';
path_Phantom_ShimOFF = 'F:\hamaguchi\20251215\dual_echo\24\1_original_data';
path_Human_ShimON    = 'F:\hamaguchi\20251215\dual_echo\27\1_original_data';
path_Human_ShimOFF   = 'F:\hamaguchi\20251215\dual_echo\27Z\1_original_data';

% 結果を格納する構造体を作成
Data = struct();

%% 2. データの読み込みとB0マップ計算関数（下部に定義）の実行

fprintf('--- ファントムデータの処理中 ---\n');
Data.Phantom.ON  = Process_DICOM_to_B0(path_Phantom_ShimON, 'Phantom Shim ON');
Data.Phantom.OFF = Process_DICOM_to_B0(path_Phantom_ShimOFF, 'Phantom Shim OFF');

fprintf('--- ヒトデータの処理中 ---\n');
Data.Human.ON    = Process_DICOM_to_B0(path_Human_ShimON, 'Human Shim ON');
Data.Human.OFF   = Process_DICOM_to_B0(path_Human_ShimOFF, 'Human Shim OFF');

%% 3. シムデータ（シム成分）の抽出とまとめ
% シム成分 = (シムありのB0マップ) - (シムなしのB0マップ)
% ※注意: ヒトの場合、撮影間に動きがあると単純な引き算ではアーチファクトが出ます。
% 必要に応じてCoregistration（位置合わせ）を行ってください。

% ファントムのシム成分計算
Data.Phantom.ShimField = Data.Phantom.ON.B0Map - Data.Phantom.OFF.B0Map;

% ヒトのシム成分計算
Data.Human.ShimField   = Data.Human.ON.B0Map - Data.Human.OFF.B0Map;


%% 4. 結果の可視化（中央スライスを表示）

figure('Name', 'Shim Data Analysis', 'Position', [100, 100, 1200, 800]);

% 表示するスライス番号（中央）
slice_idx = round(size(Data.Phantom.ON.B0Map, 3) / 2);

% --- ファントムの表示 ---
subplot(2, 3, 1);
imagesc(Data.Phantom.OFF.B0Map(:, :, slice_idx)); axis image off; colorbar;
title('Phantom: Shim OFF (B0 Hz)');
caxis([-100 100]); % clim -> caxis に変更

subplot(2, 3, 2);
imagesc(Data.Phantom.ON.B0Map(:, :, slice_idx)); axis image off; colorbar;
title('Phantom: Shim ON (B0 Hz)');
caxis([-100 100]); % clim -> caxis に変更

subplot(2, 3, 3);
imagesc(Data.Phantom.ShimField(:, :, slice_idx)); axis image off; colorbar;
title('Phantom: Shim Contribution (Diff)');
colormap(gca, 'jet');
% ここは差分なのでレンジを指定しないか、必要に応じて caxis([-50 50]) 等を追加

% --- ヒトの表示 ---
% ヒトはマトリクスサイズが異なる場合があるので、再計算
h_slice = round(size(Data.Human.ON.B0Map, 3) / 2);

subplot(2, 3, 4);
imagesc(Data.Human.OFF.B0Map(:, :, h_slice)); axis image off; colorbar;
title('Human: Shim OFF (B0 Hz)');
caxis([-200 200]); % clim -> caxis に変更

subplot(2, 3, 5);
imagesc(Data.Human.ON.B0Map(:, :, h_slice)); axis image off; colorbar;
title('Human: Shim ON (B0 Hz)');
caxis([-200 200]); % clim -> caxis に変更

subplot(2, 3, 6);
imagesc(Data.Human.ShimField(:, :, h_slice)); axis image off; colorbar;
title('Human: Shim Contribution (Diff)');
colormap(gca, 'jet');

% もし sgtitle もエラーになる場合（R2018a以前）は、この行を削除または suptitle に変更してください
sgtitle('B0 Mapping and Shim Field Analysis');


%% === ローカル関数: 読み込みとB0計算を一括で行う ===
function Out = Process_DICOM_to_B0(dicom_path, label)
    disp(['Reading: ' label ' ...']);
    
    % ご提示のRead_DICOM関数を使用
    [iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, files] = Read_DICOM(dicom_path);
    
    % --- B0マップの計算 ---
    % iFieldは (x, y, z, echoes) の4次元配列と想定
    % 2つのエコー間の位相差から磁場を計算します
    
    if size(iField, 4) < 2
        warning('エコーが1つしかありません。正確なB0マップ計算には2エコー以上が必要です。');
        B0Map = zeros(size(iField, 1:3));
    else
        % 第1エコーと第2エコーを使用
        img1 = iField(:, :, :, 1);
        img2 = iField(:, :, :, 2);
        
        % 複素共役をとって掛け合わせることで位相差を抽出
        % angle( img2 * conj(img1) ) = phase2 - phase1
        PhaseDiff = angle(img2 .* conj(img1));
        
        % 位相アンラップ（Unwrapping）
        % ※簡易的な2Dアンラップの例です。本来は3D（PRELUDE等）が望ましいですが、
        % MATLAB標準関数で簡易実装します。
        UnwrappedPhase = zeros(size(PhaseDiff));
        for z = 1:size(PhaseDiff, 3)
            UnwrappedPhase(:, :, z) = unwrap(unwrap(PhaseDiff(:, :, z), [], 1), [], 2);
        end
        
        % Hzへの変換: B0(Hz) = ΔPhase / (2π * ΔTE)
        % delta_TEの単位が秒であることを確認（Read_DICOM内では秒になっているはずです）
        if delta_TE == 0
            % delta_TEが取得できなかった場合のフォールバック（例: TE差がDICOMタグから読めない場合）
            % disp('delta_TE is 0 or unknown. Creating raw Phase map.');
            B0Map = UnwrappedPhase; % Hz変換できず
        else
            B0Map = UnwrappedPhase / (2 * pi * delta_TE);
        end
    end
    
    % 出力構造体にまとめる
    Out.iField = iField;
    Out.B0Map  = B0Map;
    Out.voxel_size = voxel_size;
    Out.CF = CF; % 中心周波数 (Hz) -> ppm変換したい場合に (B0_Hz / CF * 1e6) で計算可能
end