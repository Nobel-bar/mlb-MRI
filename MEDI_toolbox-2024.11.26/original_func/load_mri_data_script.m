function [img_4d, phase_4d] = load_mri_data(base_path, mag_name, phase_name, dims, varargin)
% LOAD_MRI_DATA MRIのバイナリデータと関連するMATファイルを読み込む関数
%
% [img_4d, phase_4d] = load_mri_data(base_path, mag_name, phase_name, dims)
% [img_4d, phase_4d] = load_mri_data(..., 'DS_FACTOR', 1, 'MaskPath', '...', 'Precision', 'single')
%
% 入力:
%   base_path  : バイナリファイルがあるフォルダパス
%   mag_name   : Magnitudeファイル名 (.binなど)
%   phase_name : Phaseファイル名 (.binなど)
%   dims       : 元データの次元 [Height, Width, Depth, (Vol)]
%
% オプション (Name-Valueペア):
%   'DS_FACTOR': ダウンサンプリング係数 (デフォルト: 1)
%   'MaskPath' : phase.mat 等があるパス (デフォルト: base_pathと同じ)
%   'Precision': freadの精度 (デフォルト: 'double')

    % --- オプション引数の解析 ---
    p = inputParser;
    addRequired(p, 'base_path', @ischar);
    addRequired(p, 'mag_name', @ischar);
    addRequired(p, 'phase_name', @ischar);
    addRequired(p, 'dims', @isnumeric);
    addParameter(p, 'DS_FACTOR', 1, @isnumeric);
    addParameter(p, 'MaskPath', base_path, @ischar); % 指定なければbase_pathと同じ
    addParameter(p, 'Precision', 'double', @ischar); % デフォルトはsingle(float32)

    parse(p, base_path, mag_name, phase_name, dims, varargin{:});
    
    % 変数への展開
    ds_factor = p.Results.DS_FACTOR;
    mask_path = p.Results.MaskPath;
    precision = p.Results.Precision;
    
    % --- 1. Magnitudeデータの読み込み ---
    mag_full_path = fullfile(base_path, mag_name);
    fid = fopen(mag_full_path, 'rb');
    if fid == -1
        error('Magnitudeファイルが開けません: %s', mag_full_path);
    end
    raw_mag = fread(fid, inf, precision);
    fclose(fid);
    
    % Reshape (オリジナルサイズ)
    mag_orig = reshape(raw_mag, dims);
    
    % --- 2. Phaseデータの読み込み (バイナリ) ---
    % ※元のコードでは読み込んでいますが、最終的なphase_4dにはphase.matのiFreqを使っています。
    %   ここでは一応読み込み処理を残しますが、必要なければ削除可能です。
    phase_full_path = fullfile(base_path, phase_name);
    fid = fopen(phase_full_path, 'rb');
    if fid == -1
        warning('Phaseファイルが開けません: %s', phase_full_path);
        % バイナリがない場合は何もしない（下のtry-catchに委ねる）
    else
        % 読み込むが、元のロジック通りならここは使われない可能性がある
        fread(fid, inf, precision); 
        fclose(fid);  
    end
    
   % --- 3. Phase/Freq情報 (.mat)  の読み込みと結合 ---
    % 元コードのロジック: バイナリではなく mask_path 内の phase.mat (iFreq) を優先使用
    try
        tmp = load(fullfile(mask_path, 'phase.mat'), 'iFreq');
        iFreq = tmp.iFreq;
        % PDF.mat (RDF) も読み込んでいましたが、出力に使われていないため省略します
        % 必要ならここで読み込んでください
    catch
        warning('phase.mat (iFreq) が見つかりません。ゼロで埋めます。');
        iFreq = zeros(dims);
    end
    
    % --- 4. ダウンサンプリング処理と出力 ---
    % Magnitude
    img_4d = mag_orig(1:ds_factor:end, 1:ds_factor:end, :);
    
    % Phase (iFreqを使用)
    phase_4d = iFreq(1:ds_factor:end, 1:ds_factor:end, :);

end


% 
% % 関数呼び出し
% [iMag_4D, iPhase_4D] = load_mri_data(...
%     load_base_path, ...      % base_path
%     mag_filename, ...        % mag_filename
%     phase_filename, ...      % phase_filename
%     dims_orig, ...                   % dims (必須)
%     'DS_FACTOR', 1, ...              % オプション
%     'MaskPath', load_mask_path, ... % オプション (別フォルダの場合)
%     'Precision', precision ...       % オプション
% );
