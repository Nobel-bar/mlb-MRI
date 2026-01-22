%% WrappedPhase 4スライス表示スクリプト
clear variables; close all; clc;

% =========================================================================
% 1. 設定 (パスを確認してください)
% =========================================================================
% データがある親フォルダ (例: 'F:\...\dual_echo' や 'C:\Users\...\dual_echo')
base_dir ='C:\Users\hamaguchi\project\b0_mapping_project\data\20251215\dual_echo';
target_id = '27';

% ファイルパスの構築
mat_file_path = fullfile(base_dir, target_id, '3_qsm_data', 'RDF.mat');

% =========================================================================
% 2. データの読み込み
% =========================================================================
if ~exist(mat_file_path, 'file')
    error('ファイルが見つかりません: %s', mat_file_path);
end

fprintf('データを読み込んでいます: %s ...\n', mat_file_path);
% 変数の中身を確認してから読み込む
vars = whos('-file', mat_file_path);
if ismember('WrappedPhase', {vars.name})
    load(mat_file_path, 'WrappedPhase');
    target_data = WrappedPhase;
    var_name = 'WrappedPhase';
elseif ismember('RDF', {vars.name})
    % WrappedPhaseがない場合、RDF変数を代用
    fprintf('WrappedPhaseが見つかりません。"RDF" を読み込みます。\n');
    load(mat_file_path, 'RDF');
    target_data = RDF;
    var_name = 'RDF';
else
    error('指定されたファイルに WrappedPhase または RDF 変数が含まれていません。');
end

% =========================================================================
% 3. 4スライスの表示
% =========================================================================
[nx, ny, nz] = size(target_data);

% 表示するスライス位置を自動計算 (全体の20%, 40%, 60%, 80%の位置)
slice_indices = round(linspace(1, nz, 6)); 
slices_to_show = slice_indices(2:5); 

figure('Name', ['Check: ' var_name ' (Case ' target_id ')'], ...
       'Color', 'w', 'Position', [100, 100, 1000, 800]);

for i = 1:4
    sl = slices_to_show(i);
    
    subplot(2, 2, i);
    % 画像を90度回転させて表示 (MATLABの仕様対策)
    imagesc(rot90(target_data(:,:,sl), 1));
    
    axis image off; 
    colormap gray;
    colorbar;
    
    % コントラスト調整 (ppmデータの値域に合わせて調整)
    % 脳が見えやすい範囲: -0.1 ～ 0.1 程度 (データによっては -0.5~0.5)
    caxis([-3.2, 3.2]); 
    
    title(sprintf('Slice %d / %d', sl, nz), 'FontSize', 12);
end

sgtitle(['Variable: ' var_name ' (Case ' target_id ')'], 'FontSize', 14, 'FontWeight', 'bold');
fprintf('表示完了。\n');