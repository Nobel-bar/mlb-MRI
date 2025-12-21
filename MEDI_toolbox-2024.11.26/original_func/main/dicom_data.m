% metadeta やDICOMの情報を抜き出す


% ファイル名の指定
% シム電流値(タグ)を参照するための代表DICOMファイル
dicomFolder =  'F:\hamaguchi\20251215\dual_echo\25\1_original_data';  % シム情報を抽出したいDICOMファイル


%% 3. DICOMヘッダからシム電流値を抽出
files = dir(fullfile(dicomFolder));
if isempty(files)
    error('指定フォルダにDICOMファイルが見つかりません。');
end

shimValues = struct(); % 結果格納用
foundShim = false;

disp(['フォルダ内のファイル数: ' num2str(length(files))]);
disp('シム情報を検索中...');

disp(['フォルダ内のファイル数: ' num2str(length(files))]);
disp('シム情報を検索中...');
for i = 1:length(files)
    filePath = fullfile(files(i).folder, files(i).name);
    try
        info = dicominfo(filePath);
        disp(info.PatientID); % IDなどを確認できます
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

