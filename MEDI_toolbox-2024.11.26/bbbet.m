function [status, result] = bet(infile, outfile, options)
% BET - Wrapper function to call FSL's bet from Windows MATLAB via WSL
%
% Usage: bet(infile, outfile, options)
% Example: bet('Mag.nii', 'Brain', '-m -f 0.5')

    % オプションがない場合のデフォルト
    if nargin < 3
        options = '-m -f 0.5';
    end

    % 1. WindowsのパスをWSLのパス形式(/mnt/c/...)に変換する関数
    function wslPath = win2wsl(winPath)
        % 絶対パスを取得
        if ~java.io.File(winPath).isAbsolute()
            winPath = fullfile(pwd, winPath);
        end
        % バックスラッシュをスラッシュに置換
        winPath = strrep(winPath, '\', '/');
        % ドライブレター (C:/...) をマウントパス (/mnt/c/...) に置換
        driveExpr = '^([a-zA-Z]):';
        tokens = regexp(winPath, driveExpr, 'tokens');
        if ~isempty(tokens)
            driveLetter = lower(tokens{1}{1});
            wslPath = regexprep(winPath, driveExpr, ['/mnt/' driveLetter]);
        else
            error('パスの変換に失敗しました: %s', winPath);
        end
    end

    % 入出力パスを変換
    wsl_infile = win2wsl(infile);
    % 出力ファイル名に拡張子がない場合への対応などはFSLが自動で行うが
    % パス変換のために一時的に処理
    wsl_outfile = win2wsl(outfile);

    % 2. WSLコマンドの組み立て
    % FSLの環境設定を読み込んでから実行する
    % ※ /usr/local/fsl は標準的なインストール先。異なる場合は変更してください。
    fsl_path = '/usr/local/fsl'; 
    
    cmd_str = sprintf('wsl bash -c "export FSLDIR=%s; . ${FSLDIR}/etc/fsl/fsl.sh; PATH=${FSLDIR}/bin:${PATH}; bet ''%s'' ''%s'' %s"', ...
                      fsl_path, wsl_infile, wsl_outfile, options);

    % 3. 実行
    disp(['Running FSL BET via WSL...']);
    disp(['Command: ' cmd_str]);
    
    [status, result] = system(cmd_str);

    if status == 0
        disp('BET finished successfully.');
    else
        disp('Error occurring in BET:');
        disp(result);
    end
end