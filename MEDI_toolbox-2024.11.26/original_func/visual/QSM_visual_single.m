%================================================================
% QSM (磁化率マップ) 単一スライス 2D/3D 表示プログラム
% 
% 概要:
%   指定した1スライスのQSMを、左側に2D (imagesc)、
%   右側に3D (mesh) で並べて表示します。
%
%   - caxis を使用します (古いMATLABバージョン対応)。
%   - ウィンドウを最大化し、プロットの余白を調整します。
%
% 依存ファイル:
%   - 'QSM.mat' ('QSM' 変数を含む)
%   - 'Mask.mat' ('Mask' 変数を含む)
%================================================================
clear variables;
close all;

%% --- 1. データ読み込み ---
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H_local = 'C:\Users\hamaguchi\Downloads\matlab\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
image_file_0 = image_file_2DGE_1_2_Rotate_H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_path = fullfile(image_file_0, image_file_3);

fprintf('データを読み込んでいます...\n');
try
    load(fullfile(load_path, 'QSM.mat'), 'QSM'); 
    load(fullfile(load_path, 'Mask.mat'), 'Mask'); 
catch ME
    fprintf('ファイルの読み込みに失敗しました。\n');
    fprintf('QSM.mat と Mask.mat が %s に存在するか確認してください。\n', load_path);
    rethrow(ME);
end

if ~exist('Mask', 'var')
    warning('変数 "Mask" が見つかりません。マスクなしで続行します。');
    Mask = ones(size(QSM)); 
end

fprintf('データの読み込みが完了しました。\n');


%% --- 2. 表示設定 ---

% 1. 表示するスライス番号を指定
total_slices = size(QSM, 3);
slice_to_display = round(total_slices / 2); % (例: 中央スライス)

% スライス番号が有効かチェック
if slice_to_display < 1 || slice_to_display > total_slices
    fprintf('エラー: スライス番号 %d は無効です。1から %d の間で指定してください。\n', ...
            slice_to_display, total_slices);
    return;
end

% 2. QSM表示用のカラースケール
qsm_clim = [-0.5, 0.5]; % [ppm]単位。データに合わせて調整してください。


fprintf('スライス %d の 2D/3D 表示を生成します...\n', slice_to_display);

%% --- 3. 2D/3D 比較プロット ---

% 1. 比較用の新しいFigureを作成 (最大化)
figure('Name', ['QSM 2D/3D Comparison - Slice ' num2str(slice_to_display)], ...
       'WindowState', 'maximized');
sgtitle(sprintf('QSM 2D (left) vs 3D (right) - Slice %d', slice_to_display), 'FontWeight', 'bold');

% 2. 現在のスライスのデータを準備
qsm_slice_data = QSM(:, :, slice_to_display);
current_mask = Mask(:, :, slice_to_display);

% --- 3. 左側: 2D (imshow) 表示 ---
ax1 = subplot(1, 2, 1);
imagesc(ax1, qsm_slice_data, 'AlphaData', current_mask);

caxis(ax1, qsm_clim); % カラースケールを固定
axis equal tight;
axis off;
colormap(ax1, 'gray');
title('2D Display (imagesc)');
colorbar(ax1); % 2Dプロット用のカラーバー


% --- 4. 右側: 3D (mesh) 表示 ---
ax2 = subplot(1, 2, 2);

% 3D mesh用にデータを準備 (マスク外をNaNに)
qsm_slice_3d = qsm_slice_data;
qsm_slice_3d(current_mask == 0) = NaN;

mesh(ax2, qsm_slice_3d);
axis tight;
daspect([50 50 1]); % Z軸(高さ)を強調 (X:50, Y:50 に対し Z:1)
colormap(ax2, 'default');
title('3D Display (mesh)');
xlabel('X Index'); ylabel('Y Index'); zlabel('Susceptibility [ppm]');
colorbar(ax2); % 3Dプロット用のカラーバー

% 3Dプロットのカラースケールも2Dと合わせる
caxis(ax2, qsm_clim);


% --- 5. プロット位置の調整 (余白を削減) ---
% [left, bottom, width, height]
ax1.Position = [0.05, 0.05, 0.43, 0.85];
ax2.Position = [0.52, 0.05, 0.43, 0.85];

fprintf('表示が完了しました。\n');