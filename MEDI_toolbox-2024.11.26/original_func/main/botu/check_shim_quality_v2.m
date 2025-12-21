%% シム精度評価プログラム (MATファイル指定版)
% 強度画像(iMag)と位相/周波数画像(iFreq)のMATファイルを読み込み、
% DICOMのシム電流値と共に評価します。
clear; clc; close all;

%% 1. ファイルパスの設定 (ここを変更してください)
% -----------------------------------------------------------
% 解析対象のMATファイル (強度と位相が別々の場合)

magFileName =  'F:\hamaguchi\20251215\dual_echo\25\3_qsm_data/Mask.mat';  % 強度と位相が入ったMATファイル
phaseFileName =  'F:\hamaguchi\20251215\dual_echo\25\3_qsm_data/phase.mat';  % 強度と位相が入ったMATファイル

% シム電流値(タグ)を参照するための代表DICOMファイル
dicomFolder =  'F:\hamaguchi\20251215\dual_echo\25\1_original_data';  % シム情報を抽出したいDICOMファイル

% シム電流の飽和判定リミット (16bit DAC想定)
SHIM_LIMIT = 30000; 
% -----------------------------------------------------------

%% 2. 画像データの読み込み (ご指定のコード)
fprintf('データを読み込んでいます...\n');

if exist(magFileName, 'file')
    data = load(magFileName);
    disp(['Loaded Magnitude: ' magFileName]);

    if isfield(data, 'iMag')
        img_mag = data.iMag;
    else
        error('MATファイル内に変数 iMag が見つかりません。');
    end
    
    if exist(phaseFileName, 'file')
        data = load(phaseFileName);
        disp(['Loaded Phase: ' phaseFileName]);
        
        if isfield(data, 'iFreq')
            img_phase = data.iFreq;
        else
            error('MATファイル内に変数 iFreq が見つかりません。');
        end
    else
        error(['位相MATファイルが見つかりません: ' phaseFileName]);
    end
else
    error(['強度MATファイルが見つかりません: ' magFileName]);
end


%% 3. DICOMヘッダからシム電流値を抽出
files = dir(fullfile(dicomFolder));
if isempty(files)
    error('指定フォルダにDICOMファイルが見つかりません。');
end

shimValues = struct(); % 結果格納用
foundShim = false;

disp(['フォルダ内のファイル数: ' num2str(length(files))]);
disp('シム情報を検索中...');
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

% if exist(dicomFileName, 'file')
%     info = dicominfo(dicomFileName);
%     shimValues = struct();
%     
%     % Private Tag (0029,1022) の抽出
%     if isfield(info, 'Private_0029_1022')
%         rawStr = info.Private_0029_1022;
%         % uint8配列なら文字列に変換
%         if isa(rawStr, 'uint8') || isa(rawStr, 'int8')
%             shimStr = char(rawStr');
%         else
%             shimStr = string(rawStr);
%         end
%         disp(['DICOM Shim String: ' shimStr]);
%         shimValues = parseShimString(shimStr);
%     else
%         warning('DICOMヘッダにシム情報(0029,1022)が見つかりません。');
%     end
% else
%     error(['DICOMファイルが見つかりません: ' dicomFileName]);
% end

%% 4. 評価と判定

% --- 4-1. シム電流値の飽和チェック ---
fprintf('\n=== シム電流値(DAC) 判定 ===\n');
isSaturated = false;
fnames = fieldnames(shimValues);

if isempty(fnames)
    fprintf('シムデータなし\n');
else
    for k = 1:length(fnames)
        ch = fnames{k};
        val = shimValues.(ch);
        
        % 状態判定
        if abs(val) >= SHIM_LIMIT
            status = 'NG (飽和)';
            isSaturated = true;
        elseif abs(val) > 25000
            status = '注意 (高出力)';
        else
            status = 'OK';
        end
        fprintf('  %-4s : %6d -> %s\n', ch, val, status);
    end
end

% --- 4-2. 画像(B0マップ)の均一性チェック ---
fprintf('\n=== B0均一性 (iFreq) 判定 ===\n');

% マスク作成 (iMagの信号がある部分のみ抽出)
% 信号強度の最大値の10%以上を有効領域とする
mask = img_mag > (max(img_mag(:)) * 0.1);

% マスク領域内の位相(周波数)データ抽出
roi_freq = double(img_phase(mask));

% 統計量の計算
freq_std = std(roi_freq);
freq_range = max(roi_freq) - min(roi_freq);

fprintf('解析対象ピクセル数: %d\n', length(roi_freq));
fprintf('周波数ムラ (StdDev): %.4f (Hz)\n', freq_std);
fprintf('周波数レンジ (Range): %.4f (Hz)\n', freq_range);

% 判定 (閾値はHz単位と仮定。目的によりますが数Hz～10Hz以内なら優秀)
THRESHOLD_STD = 15.0; % 閾値 (必要に応じて調整してください)

if freq_std < THRESHOLD_STD
    fprintf('[判定: OK] 磁場均一性は良好です。\n');
else
    fprintf('[判定: 注意] 磁場不均一が大きいです (StdDev > %.1f)\n', THRESHOLD_STD);
    if isSaturated
        fprintf('           -> 原因: シムコイルの飽和が疑われます。\n');
    end
end

%% 5. 結果の可視化 (3D対応修正版)
% 3次元データの場合、中央のスライスを抜き出して表示します
figure('Name', 'Shim Quality Check v3', 'Color', 'w', 'Position', [100, 100, 1000, 500]);

% --- スライス選択 ---
if ndims(img_mag) == 3
    sliceIdx = round(size(img_mag, 3) / 2);
    disp(['3次元データ検出: スライス ' num2str(sliceIdx) ' を表示します。']);
    
    show_mag = img_mag(:, :, sliceIdx);
    show_phase = img_phase(:, :, sliceIdx);
else
    show_mag = img_mag;
    show_phase = img_phase;
end

% --- 強度画像 ---
subplot(1, 2, 1);
imshow(show_mag, []);
title('Magnitude (Mid-Slice)');
colorbar;

% --- 位相(周波数)画像 ---
subplot(1, 2, 2);
% 表示用のマスク作成（スライス単体に対して）
slice_mask = show_mag > (max(show_mag(:)) * 0.1);

% マスク外をNaNにして背景を見やすくする
plot_phase = double(show_phase);
plot_phase(~slice_mask) = NaN; 

imagesc(plot_phase);
axis image off;
title({'Frequency Map (Mid-Slice)', ['StdDev(Vol): ' num2str(freq_std, '%.2f') ' Hz']});
colormap(gca, 'jet');
colorbar;


%% 補助関数
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