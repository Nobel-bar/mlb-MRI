%% ========================================================================
%  MAT to RAW Converter for ImageJ (Diff Analysis)
%  Target: iFreq_ppm, iFreq, and (Wrapped - iFreq) -> .raw
% ========================================================================
clear variables; close all; clc;

% --- 1. 設定 ---
% パスはご自身の環境に合わせてください
base_dir = 'C:\Users\hamaguchi\project\b0_mapping_project\data\20251215\dual_echo';
target_ids = 24:33; 

fprintf('ImageJ用 RAWデータ変換 (差分計算付き) を開始します...\n');

for idx = target_ids
    case_name = num2str(idx);
    
    % パス設定
    case_dir = fullfile(base_dir, case_name);
    mat_path = fullfile(case_dir, '3_qsm_data', 'RDF.mat');
    
    if ~exist(mat_path, 'file')
        fprintf('Skip: %s が見つかりません\n', mat_path);
        continue;
    end
    
    % --- 2. データの読み込み ---
    % 必要な変数をロード
    data = load(mat_path, 'iFreq_ppm', 'iFreq');
    
    if ~isfield(data, 'iFreq_ppm') || ~isfield(data, 'iFreq')
        fprintf('Skip: Case %s に必要な変数(iFreq_ppm/iFreq)がありません\n', case_name);
        continue;
    end
    
    % --- 3. 差分の計算 ---
    % iFreq_ppm (-pi~pi) から iFreq (背景除去済みアンラップ位相) を引く
    % これにより「除去された背景磁場 」が残ります
    Diff_Phase = data.iFreq_ppm - data.iFreq;
    
    % --- 4. RAW保存の実行 ---
    % 保存先フォルダ
    save_dir = fullfile(case_dir, '3_qsm_data');
    
    % 画像サイズ取得 (ログ表示用)
    [w, h, d] = size(data.iFreq_ppm);
    
    % (1) iFreq_ppm の保存
    save_raw_for_imagej(data.iFreq_ppm, fullfile(save_dir, sprintf('Img_UnWrapped.raw')));
    
    % (2) iFreq の保存
    save_raw_for_imagej(data.iFreq, fullfile(save_dir, sprintf('Img_iFreq.raw')));
    
    % (3) 差分 (Diff) の保存
    save_raw_for_imagej(Diff_Phase, fullfile(save_dir, sprintf('Img_Diff.raw')));
    
    fprintf('Saved Case %s: Size=[%d x %d x %d]\n', case_name, w, h, d);
end

fprintf('\n全変換完了。\n');
fprintf('--------------------------------------------------\n');
fprintf('【ImageJ 入力パラメータ (Import > Raw)】\n');
fprintf('  Image Type:      32-bit Real\n');
fprintf('  Width:           %d\n', w);
fprintf('  Height:          %d\n', h);
fprintf('  Number of images:%d\n', d);
fprintf('  Little-endian:   [Check ON]\n');
fprintf('--------------------------------------------------\n');

% =========================================================================
%  サブ関数: ImageJ用に形式を整えて保存
% =========================================================================
function save_raw_for_imagej(img_data, filepath)
    % 1. single (float32) に変換
    img_single = single(img_data);
    
    % 2. 次元の入れ替え (MATLAB:列優先 -> ImageJ:行優先 への対策)
    % これをしないとImageJで見たときに画像が90度回転・反転してしまいます
    % img_permuted = permute(img_single, [2 1 3]);
    
    % 3. 書き込み
    fid = fopen(filepath, 'w', 'l'); % Little-Endian
    if fid == -1
        error('ファイル書き込みエラー: %s', filepath);
    end
    fwrite(fid, img_single, 'single');
    fclose(fid);
end