%% 1. 設定とファイルの読み込み
clear; clc; close all;

% シム値(DAC値)のハードウェア限界 (推測値)
% 16bit符号付き整数の場合、最大は32767です。これに近いと飽和の危険があります。
SHIM_LIMIT = 30000;

% --- ファイル名の指定 (適宜変更してください) ---
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
magFileName =  'F:\hamaguchi\20251215\dual_echo\25\3_qsm_data/Mask.mat';  % 強度と位相が入ったMATファイル
phaseFileName =  'F:\hamaguchi\20251215\dual_echo\25\3_qsm_data/phase.mat';  % 強度と位相が入ったMATファイル
dicomFolder =  'F:\hamaguchi\20251215\dual_echo\25\1_original_data';  % シム情報を抽出したいDICOMファイル
% ---------------------------------------------

% 1-1. MATファイルの読み込み
if exist(magFileName, 'file')
    data = load(magFileName);
    disp(['Loaded: ' magFileName]);

    if isfield(data, 'iMag')
        img_mag = data.iMag;
    else
        error('MATファイル内に指定の変数が見つかりません。変数名を確認してください。');
    end
    
    data = load(phaseFileName);
    % 変数名の確認 (MATファイル内の変数名に合わせて変更してください)
    if isfield(data, 'iFreq')
        img_phase = data.iFreq;
    else
        error('MATファイル内に指定の変数が見つかりません。変数名を確認してください。');
    end
else
    error('MATファイルが見つかりません。');
end



%% 3. DICOMフォルダからシム情報を抽出
files = dir(fullfile(dicomFolder));
if isempty(files)
    error('指定フォルダにDICOMファイルが見つかりません。');
end

shimValues = struct(); % 結果格納用
foundShim = false;

disp(['フォルダ内のファイル数: ' num2str(length(files))]);
disp('シム情報を検索中...');

% フォルダ内のファイルを順にチェック (通常、シリーズ内のシム値は同一なので1つ見つかればOK)
for i = 1:length(files)
    filePath = fullfile(files(i).folder, files(i).name);
    try
        info = dicominfo(filePath);
        
        % Private Tag (0029,1022) を探す
        % MATLABでは通常 'Private_0029_1022' というフィールド名になります
        if isfield(info, 'Private_0029_1022')
            rawStr = info.Private_0029_1022;
            
            % uint8配列(バイト列)として読み込まれた場合は文字列に変換
            if isa(rawStr, 'uint8') || isa(rawStr, 'int8')
                shimStr = char(rawStr');
            else
                shimStr = string(rawStr);
            end
            
            % 文字列が "B0=..." の形式か確認
            if contains(shimStr, 'B0=')
                disp(['発見 (' files(i).name '): ' shimStr]);
                shimValues = parseShimString(shimStr);
                foundShim = true;
                break; % 1つ見つかったらループを抜ける
            end
        end
    catch
        % 読み込めないファイルはスキップ
        continue;
    end
end

if ~foundShim
    error('フォルダ内のDICOMファイルからシム情報 (0029,1022) を抽出できませんでした。');
end

%% 4. シム値の健全性判定 (Saturation Check)

fprintf('\n=== シム電流値(DAC値) 判定レポート ===\n');
fieldNames = fieldnames(shimValues);
isSaturated = false;

for k = 1:length(fieldNames)
    chName = fieldNames{k};
    val = shimValues.(chName);
    
    % 限界値に対する割合
    ratio = abs(val) / 32767 * 100; % 16bit max想定
    
    fprintf('  %-4s : %6d  (Output: %3.1f%%) ', chName, val, ratio);
    
    if abs(val) >= SHIM_LIMIT
        fprintf('-> [警告: 飽和の可能性あり]\n');
        isSaturated = true;
    elseif abs(val) > 25000
        fprintf('-> [注意: 高出力]\n');
    else
        fprintf('-> [OK]\n');
    end
end

if isSaturated
    fprintf('\n[総合判定: NG] 一部のシムコイルが出力限界に近いため、補正しきれていない可能性があります。\n');
else
    fprintf('\n[総合判定: OK] シム値は制御範囲内に収まっています。\n');
end

%% 5. 画像ベースのB0均一性判定 (MATファイルデータ使用)

fprintf('\n=== 画像ベース B0均一性判定 ===\n');

% マスク作成 (背景ノイズ除去)
mask = img_mag > (max(img_mag(:)) * 0.1); 
brain_phase = img_phase(mask);

% 統計量
phase_std = std(double(brain_phase(:)));
phase_range = range(double(brain_phase(:)));

fprintf('ROI内 位相標準偏差: %.4f\n', phase_std);

% 閾値判定 (実験に合わせて調整してください)
THRESHOLD_STD = 2.5; 

if phase_std < THRESHOLD_STD
    fprintf('[判定: OK] 画像上の磁場均一性は良好です。\n');
else
    fprintf('[判定: 注意] 磁場不均一がやや大きいです (StdDev > %.1f)。\n', THRESHOLD_STD);
    if isSaturated
        fprintf('           -> シム値の飽和が原因である可能性が高いです。\n');
    end
end

%% --- 補助関数: 文字列パース ---
function shimStruct = parseShimString(strIn)
    % 入力例: "B0=27017,X2Y2=-78,..." を構造体に変換
    shimStruct = struct();
    
    % カンマで分割
    parts = strsplit(strIn, ',');
    
    for i = 1:length(parts)
        item = strtrim(parts{i});
        if contains(item, '=')
            kv = strsplit(item, '=');
            key = strtrim(kv{1});
            val = str2double(kv{2});
            
            % 有効な数値なら構造体に格納
            if ~isnan(val)
                shimStruct.(key) = val;
            end
        end
    end
end