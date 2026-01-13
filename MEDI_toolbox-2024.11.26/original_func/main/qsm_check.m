%% --- QSM画像表示用スクリプト (修正版) ---
clear variables; close all; clc;

% ==========================================
% 設定項目
% ==========================================
check_id = 27; % 確認したいフォルダ番号
base_dir = 'F:\hamaguchi\data\20251215\dual_echo'; % データの親フォルダ
% ==========================================

% パスの生成
case_name = num2str(check_id);
image_file_dual_echo = fullfile(base_dir, case_name);
image_file_3 = '3_qsm_data';
save_path = fullfile(image_file_dual_echo, image_file_3);
mat_file_path = fullfile(save_path, 'Reproduction_Inputs.mat');

% データの読み込み
if exist(mat_file_path, 'file')
    fprintf('データを読み込んでいます: %s\n', mat_file_path);
    data = load(mat_file_path, 'QSM');
    QSM = data.QSM;
else
    error('ファイルが見つかりません: %s', mat_file_path);
end

% 画像サイズとスライス選択
[nx, ny, nz] = size(QSM);

% スライス選択（脳がある真ん中あたりのスライスを抽出）
start_slice = round(nz * 0.25); 
end_slice   = round(nz * 0.75);
slice_indices = round(linspace(start_slice, end_slice, 9));

% 表示設定
display_range = [-0.15, 0.15]; 

figure('Name', ['QSM Result: Case ' case_name], 'Color', 'w', 'Position', [100, 100, 1000, 1000]);

for i = 1:9
    sl = slice_indices(i);
    
    subplot(3, 3, i);
    
    % 画像データの取得
    img_slice = squeeze(QSM(:, :, sl));
    
    % 【重要】向きの調整
    % 添付画像を見ると小さく横倒しになっている可能性があるため、
    % rot90の回数を変えて調整してください (1, 2, 3, -1 など)
    % ここでは「反時計回りに90度回転」させてみます
    img_show = rot90(img_slice, 1); 
    
    imagesc(img_show); 
    
    axis image; 
    axis off; 
    colormap gray;
    
    % --- 修正箇所: 古いバージョン用に caxis を使用 ---
    caxis(display_range); 
    % ------------------------------------------------
    
    title(sprintf('Slice %d', sl));
end

% カラーバー
h = colorbar;
h.Position = [0.92 0.1 0.02 0.8]; 
h.Label.String = 'Susceptibility (ppm)';