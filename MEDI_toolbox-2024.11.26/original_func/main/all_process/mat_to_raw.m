%% ========================================================================
%  MAT to RAW Converter for ImageJ
%  Target: WrappedPhase in RDF.mat -> .raw (Float32)
% ========================================================================
clear variables; close all; clc;

% --- 1. 設定 (ご自身の環境に合わせてください) ---
base_dir = 'C:\Users\hamaguchi\project\b0_mapping_project\data\20251215\dual_echo';
target_ids = 24:34; 

fprintf('ImageJ用 RAWデータ変換を開始します...\n');

for idx = target_ids
    case_name = num2str(idx);
    
    % パス設定
    case_dir = fullfile(base_dir, case_name);
    mat_path = fullfile(case_dir, '3_qsm_data', 'RDF.mat');
    
    % 保存ファイル名 (例: WrappedPhase_Case24.raw)
    save_filename = sprintf('WrappedPhase.raw');
    save_path = fullfile(case_dir, '3_qsm_data', save_filename);
    
    if ~exist(mat_path, 'file')
        fprintf('Skip: %s が見つかりません\n', mat_path);
        continue;
    end
    
    % --- 2. 読み込み ---
    % 必要な変数だけ読み込みます
    load(mat_path, 'WrappedPhase'); 
    
    if ~exist('WrappedPhase', 'var')
        fprintf('Skip: Case %s に WrappedPhase が保存されていません\n', case_name);
        continue;
    end
    
    % --- 3. データ形式の調整 (重要) ---
    % MATLABは double(64bit) ですが、ImageJで扱いやすい single(32bit float) に変換します
    data_to_save = single(WrappedPhase);
    
    % 【重要】次元の入れ替え
    % MATLABは列優先(Column-major)、ImageJは行優先(Row-major)的な表示をするため
    % そのまま保存すると画像が90度回転して反転して見えます。
    % permute(img, [2 1 3]) でX軸とY軸を入れ替えておくと、ImageJで正立して見えます。
    % data_to_save = permute(data_to_save, [2 1 3]);
    
    % サイズ取得 (ImageJ入力用に記録)
    [width, height, depth] = size(data_to_save);
    
    % --- 4. RAW保存 ---
    fid = fopen(save_path, 'w', 'l'); % 'l' = Little-Endian (Windows標準)
    if fid == -1
        error('ファイルを保存できませんでした: %s', save_path);
    end
    fwrite(fid, data_to_save, 'single'); % float32形式で書き込み
    fclose(fid);
    
    % --- 5. ImageJ用情報の出力 ---
    fprintf('--------------------------------------------------\n');
    fprintf('Saved: %s\n', save_filename);
    fprintf('【ImageJ 入力パラメータ】\n');
    fprintf('  Image Type:      32-bit Real\n');
    fprintf('  Width (pixels):  %d\n', width);
    fprintf('  Height (pixels): %d\n', height);
    fprintf('  Number of images:%d\n', depth);
    fprintf('  Little-endian byte order: [Check ON]\n');
end
fprintf('\n全変換完了。\n');