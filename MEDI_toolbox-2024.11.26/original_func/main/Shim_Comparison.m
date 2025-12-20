%% シム効果比較・評価プログラム (Shim ON vs OFF)
clear; clc; close all;

% --- 1. データパスの設定 ---
path_Phantom_ShimON  = 'F:\hamaguchi\20251215\dual_echo\25\1_original_data';
path_Phantom_ShimOFF = 'F:\hamaguchi\20251215\dual_echo\24\1_original_data';

% --- 2. データの抽出と解析 ---
fprintf('データを解析中...\n');
[shimON,  imgON,  infoON]  = analyze_dicom_folder(path_Phantom_ShimON);
[shimOFF, imgOFF, infoOFF] = analyze_dicom_folder(path_Phantom_ShimOFF);

if isempty(imgON) || isempty(imgOFF)
    error('画像データが正しく読み込めませんでした。パスを確認してください。');
end

%% --- 3. シム電流値 (DAC値) の比較評価 ---
fprintf('\n========================================\n');
fprintf('   シム電流値 (Shim Current) 比較レポート\n');
fprintf('========================================\n');

% 共通するシムチャンネルを取得
channels = fieldnames(shimON);
if isempty(channels)
    channels = fieldnames(shimOFF); % ONになくOFFにある場合（稀）
end

% テーブル作成用データ
t_names = {};
t_off = [];
t_on = [];
t_diff = [];
t_status = {};

limit_val = 30000; % 飽和判定閾値 (想定)

for i = 1:length(channels)
    ch = channels{i};
    
    val_on = 0;
    if isfield(shimON, ch), val_on = shimON.(ch); end
    
    val_off = 0;
    if isfield(shimOFF, ch), val_off = shimOFF.(ch); end
    
    diff_val = val_on - val_off;
    
    % 判定
    if abs(val_on) >= limit_val
        status = 'NG (飽和)';
    elseif abs(val_on) > 25000
        status = '注意 (高)';
    else
        status = 'OK';
    end
    
    % データ格納
    t_names{end+1,1} = ch;
    t_off(end+1,1)   = val_off;
    t_on(end+1,1)    = val_on;
    t_diff(end+1,1)  = diff_val;
    t_status{end+1,1} = status;
end

% テーブル表示
ResultTable = table(t_off, t_on, t_diff, t_status, ...
    'RowNames', t_names, ...
    'VariableNames', {'Shim_OFF', 'Shim_ON', 'Delta_Current', 'Status'});
disp(ResultTable);

% 総合判定（電流値）
if any(contains(t_status, 'NG'))
    fprintf('\n[電流値判定: NG] 一部のチャンネルが飽和しています。シムが不正確な可能性があります。\n');
else
    fprintf('\n[電流値判定: OK] すべてのチャンネルが制御範囲内で動作しています。\n');
end


%% --- 4. 画像によるシム精度の評価 (均一性 & 信号強度) ---
% シムが良い = 磁場が均一 = T2*減衰が少ない = 信号強度が高く、ばらつきが少ない

% ROI設定 (画像中心の50%領域を使用)
[rows, cols] = size(imgON);
centerMask = false(rows, cols);
centerMask(round(rows*0.25):round(rows*0.75), round(cols*0.25):round(cols*0.75)) = true;

% 信号強度の取得
roi_ON = double(imgON(centerMask));
roi_OFF = double(imgOFF(centerMask));

% 評価指標 1: 変動係数 (CV) = 標準偏差 / 平均 (小さいほど均一)
cv_ON = std(roi_ON) / mean(roi_ON);
cv_OFF = std(roi_OFF) / mean(roi_OFF);

% 評価指標 2: 平均信号強度 (シムが良いと信号が回復することが多い)
mean_ON = mean(roi_ON);
mean_OFF = mean(roi_OFF);
signal_gain = (mean_ON - mean_OFF) / mean_OFF * 100;

fprintf('\n========================================\n');
fprintf('   画像品質・均一性 (Homogeneity) 評価\n');
fprintf('========================================\n');
fprintf('指標                  | Shim OFF    | Shim ON     | 改善率/変化\n');
fprintf('------------------------------------------------------------\n');
fprintf('変動係数 (CV: 低=良)  | %.4f      | %.4f      | %+.1f%%\n', ...
    cv_OFF, cv_ON, (cv_ON - cv_OFF)/cv_OFF*100);
fprintf('平均信号強度 (Mean)   | %.1f      | %.1f      | %+.1f%% (Signal Gain)\n', ...
    mean_OFF, mean_ON, signal_gain);

fprintf('\n[シム精度判定]: ');
if cv_ON < cv_OFF
    fprintf('Shim ONにより均一性が向上しました (精度: 良)。\n');
else
    fprintf('Shim ONによる均一性の向上が見られません (精度: 低 または 元々良好)。\n');
end

%% --- 5. 結果の可視化 ---
figure('Name', 'Shim Comparison Analysis', 'Color', 'w', 'Position', [100, 100, 1000, 500]);

% 画像表示 (オートスケール統一)
clims = [0, max([imgON(:); imgOFF(:)])];

subplot(1, 3, 1);
imshow(imgOFF, clims);
title({'Shim OFF', ['CV: ' num2str(cv_OFF, '%.4f')]});
colorbar;

subplot(1, 3, 2);
imshow(imgON, clims);
title({'Shim ON', ['CV: ' num2str(cv_ON, '%.4f')]});
colorbar;

subplot(1, 3, 3);
% 差分画像 (Shimによる信号変化の可視化)
imgDiff = double(imgON) - double(imgOFF);
imshow(imgDiff, []);
title({'Difference (ON - OFF)', 'Brighter = Signal recovered'});
colormap(gca, 'jet');
colorbar;


%% === 関数: フォルダ解析 ===
function [shimStruct, midSliceImg, dicomInfo] = analyze_dicom_folder(folderPath)
    shimStruct = struct();
    midSliceImg = [];
    dicomInfo = [];
    
    files = dir(fullfile(folderPath)); % 拡張子は適宜調整
    if isempty(files)
        warning(['ファイルが見つかりません: ' folderPath]);
        return;
    end
    
    % 真ん中のスライスを取得 (画像の代表値として使用)
    midIdx = round(length(files) / 2);
    if midIdx == 0, midIdx = 1; end
    
    targetFile = fullfile(files(midIdx).folder, files(midIdx).name);
    
    try
        info = dicominfo(targetFile);
        img = dicomread(info);
        
        dicomInfo = info;
        midSliceImg = img;
        
        % シム値の抽出 (0029,1022)
        if isfield(info, 'Private_0029_1022')
            rawStr = info.Private_0029_1022;
            if isa(rawStr, 'uint8') || isa(rawStr, 'int8')
                shimStr = char(rawStr');
            else
                shimStr = string(rawStr);
            end
            shimStruct = parseShimString(shimStr);
        end
        
        % 見つからない場合、全ファイルを走査して探す
        if isempty(fieldnames(shimStruct))
            for k = 1:length(files)
                fpath = fullfile(files(k).folder, files(k).name);
                tmpInfo = dicominfo(fpath);
                if isfield(tmpInfo, 'Private_0029_1022')
                   rawStr = tmpInfo.Private_0029_1022;
                   if isa(rawStr, 'uint8') || isa(rawStr, 'int8'), rawStr = char(rawStr'); else, rawStr = string(rawStr); end
                   shimStruct = parseShimString(rawStr);
                   break; 
                end
            end
        end
        
    catch ME
        warning(['DICOM読み込みエラー: ' ME.message]);
    end
end

function shimStruct = parseShimString(strIn)
    shimStruct = struct();
    parts = strsplit(strIn, ',');
    for i = 1:length(parts)
        item = strtrim(parts{i});
        if contains(item, '=')
            kv = strsplit(item, '=');
            key = strtrim(kv{1});
            val = str2double(kv{2});
            if ~isnan(val)
                shimStruct.(key) = val;
            end
        end
    end
end