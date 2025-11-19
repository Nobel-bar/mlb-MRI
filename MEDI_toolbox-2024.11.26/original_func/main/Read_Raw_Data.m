function [iField, voxel_size, matrix_size, CF, delta_TE, TE, B0_dir, files] = Read_Raw_Data(mag_filepath, phase_filepath, params)
%================================================================================================================================
% Read_Raw_Data: 強度・位相のrawファイルからQSMデータを読み込みます。
%
% この関数はRead_DICOMの代替として、.mat、.nii(.nii.gz)、またはバイナリの.raw形式で保存された
% 強度画像と位相画像を読み込むために設計されています。
% 元のRead_DICOM関数と同じ出力変数を生成するため、MEDIツールボックスの
% 他の関数との互換性が確保されます。
%
% [使い方]
% 1. この関数を 'Read_Raw_Data.m' という名前で保存し、MATLABのパスに配置します。
% 2. 解析スクリプト内でこの関数を呼び出します。
%
% --- 入力 ---
% mag_filepath   (string): 強度画像ファイルへのパス。データは4D配列(x, y, z, echo)である必要があります。
% phase_filepath (string): 位相画像ファイルへのパス。データは4D配列(x, y, z, echo)で、単位はラジアンである必要があります。
% params         (struct): 撮像パラメータを含む構造体。以下のフィールドが必須です:
%   .voxel_size:  ボクセルサイズ [x, y, z] (mm単位)。
%   .matrix_size: 行列サイズ [x, y, z]。
%   .CF:          ラーモア周波数 (Hz単位)。
%   .TE:          エコー時間のベクトル (秒単位)。
%   .B0_dir:      主磁場B0の方向ベクトル [x, y, z]。
%
% --- 出力 ---
% iField:        複素数データ (4D配列: x, y, z, echo)。
% voxel_size:    入力paramsからのボクセルサイズ。
% matrix_size:   行列サイズ (データと照合済み)。
% CF:            入力paramsからのラーモア周波数。
% delta_TE:      エコー間隔。
% TE:            入力paramsからのエコー時間。
% B0_dir:        入力paramsからのB0方向。
% files:         互換性のための空の構造体。
%================================================================================================================================

fprintf('強度・位相のrawデータを読み込んでいます...\n');

%% 1. 必須入力パラメータのチェック
if ~isfield(params, 'voxel_size') || ...
   ~isfield(params, 'matrix_size') || ...
   ~isfield(params, 'CF') || ...
   ~isfield(params, 'TE') || ...
   ~isfield(params, 'B0_dir')
    error('入力 "params" 構造体に必要なフィールドがありません。voxel_size, matrix_size, CF, TE, B0_dir を指定してください。');
end

%% 2. 強度画像と位相画像の読み込み
% ファイル拡張子に基づいて読み込み方法を分岐します。
[~, ~, ext] = fileparts(mag_filepath);

if strcmpi(ext, '.nii') || endsWith(mag_filepath, '.nii.gz', 'IgnoreCase', true)
    % NIfTIファイルの読み込み
    % MATLAB R2017b以降とImage Processing Toolbox、または
    % "Tools for NIfTI and ANALYZE image" のような外部ツールボックスが必要です。
    fprintf('NIfTIファイルを読み込んでいます...\n');
    iMag = single(niftiread(mag_filepath));
    iPhase = single(niftiread(phase_filepath));
elseif strcmpi(ext, '.raw')
    % バイナリ .raw ファイルの読み込み
    fprintf('.raw バイナリファイルを読み込んでいます...\n');

    % paramsから行列サイズを決定
    dims = [params.matrix_size, length(params.TE)];
    
    % --- 強度データの読み込み ---
    fid = fopen(mag_filepath, 'rb');
    if fid == -1
        error('強度ファイルを開けませんでした: %s', mag_filepath);
    end
    % データ型を単精度浮動小数点数('single')と仮定します。実際のデータ型に応じてこの値を変更する必要があるかもしれません。
    % (例: 'double', 'int16', 'uint16')
    precision = 'double=>double'; 
    iMag_vec = fread(fid, inf, precision);
    fclose(fid);
    
    % 読み込んだデータサイズの検証
    if numel(iMag_vec) ~= prod(dims)
        error('強度ファイルのデータサイズ (%d 要素) が期待される次元 [%d x %d x %d x %d] と一致しません。', numel(iMag_vec), dims(1), dims(2), dims(3), dims(4));
    end
    iMag = reshape(iMag_vec, dims);
    
    % --- 位相データの読み込み ---
    fid = fopen(phase_filepath, 'rb');
     if fid == -1
        error('位相ファイルを開けませんでした: %s', phase_filepath);
    end
    iPhase_vec = fread(fid, inf, precision);
    fclose(fid);
    
    % 読み込んだデータサイズの検証
    if numel(iPhase_vec) ~= prod(dims)
        error('位相ファイルのデータサイズ (%d 要素) が期待される次元 [%d x %d x %d x %d] と一致しません。', numel(iPhase_vec), dims(1), dims(2), dims(3), dims(4));
    end
    iPhase = reshape(iPhase_vec, dims);
else
    % .mat ファイルと仮定
    fprintf('.mat ファイルを読み込んでいます...\n');
    mag_data = load(mag_filepath);
    mag_fields = fieldnames(mag_data);
    iMag = single(mag_data.(mag_fields{1})); % .mat ファイル内の最初の変数を読み込む

    phase_data = load(phase_filepath);
    phase_fields = fieldnames(phase_data);
    iPhase = single(phase_data.(phase_fields{1})); % .mat ファイル内の最初の変数を読み込む
end

%% 3. データ寸法の検証
if ~isequal(size(iMag), size(iPhase))
    error('強度画像と位相画像の次元が同じである必要があります。');
end

matrix_size = [size(iMag, 1), size(iMag, 2), size(iMag, 3)];
if ~isequal(matrix_size, params.matrix_size)
     warning('入力 `params.matrix_size` が読み込まれたデータの次元と一致しません。データ自体の次元を使用します。');
     % matrix_size変数は既にデータから設定されているため、警告以外の対応は不要です。
end

NumEcho = size(iMag, 4);
if NumEcho ~= length(params.TE)
    error('データ内のエコー数 (%d) がTEベクトルの長さ (%d) と一致しません。', NumEcho, length(params.TE));
end

%% 4. 複素数データの生成
% 強度(iMag)と位相(iPhase)を結合して複素数データ(iField)を作成します。
iField = iMag .* exp(1i * iPhase);

%% 5. ツールボックス互換性のための出力変数設定
voxel_size = single(params.voxel_size);
CF = params.CF;
TE = single(params.TE);
B0_dir = params.B0_dir;

% エコー間隔 (delta_TE) の計算
if length(TE) == 1
    delta_TE = TE;
else
    delta_TE = TE(2) - TE(1);
end

% 互換性のためのダミー変数
files = struct;

fprintf('Rawデータの読み込みが正常に完了しました。\n');

end

