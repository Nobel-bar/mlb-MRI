%================================================================
% QSM (磁化率マップ) 全スライス 2D/3D ペア表示プログラム
% 
% 概要:
%   全スライスを、1Figureあたり2スライスずつ表示します。
%   各Figureは 2x2 のレイアウトになります:
%   [Slice N (2D)]   [Slice N (3D)]
%   [Slice N+1 (2D)] [Slice N+1 (3D)]
%
%   - caxis を使用します (古いMATLABバージョン対応)。
%
% 依存ファイル:
%   - 'QSM.mat' ('QSM' 変数を含む)
%   - 'Mask.mat' ('Mask' 変数を含む)
%================================================================
clear variables;
close all;

%% --- 1. データ読み込み ---
% (パスはご自身の環境に合わせて設定してください)
image_file_1 = 'F:/hamaguchi/copy/20241205_RawData_H/Volunteer_Rotate_H/2DGE_0deg_H/total_slice';
image_file_2 = 'output_data';

save_path = fullfile(image_file_1, image_file_2);

fprintf('データを読み込んでいます...\n');
try
    load(fullfile(save_path, 'QSM.mat'), 'QSM'); 
    load(fullfile(save_path, 'Mask.mat'), 'Mask'); 
catch ME
    fprintf('ファイルの読み込みに失敗しました。\n');
    fprintf('QSM.mat と Mask.mat が %s に存在するか確認してください。\n', save_path);
    rethrow(ME);
end

if ~exist('Mask', 'var')
    warning('変数 "Mask" が見つかりません。マスクなしで続行します。');
    Mask = ones(size(QSM)); 
end

fprintf('データの読み込みが完了しました。\n');


%% --- 2. 表示設定 ---

% 1. 基本設定
total_slices = size(QSM, 3);
slices_per_figure = 2; % ★変更: 1Figureあたり2スライス
rows_per_slice = 1;  % 1スライスあたり1行
cols_per_slice_pair = 2; % 2D/3D の2列
rows = slices_per_figure * rows_per_slice; % Figureの総行数 ( = 2)
cols = cols_per_slice_pair;                % Figureの総列数 ( = 2)

% 2. QSM表示用のカラースケール
qsm_clim = [-0.5, 0.5]; % [ppm]単位。データに合わせて調整してください。

% 3. 画面サイズとプロット位置
screen_size = get(groot, 'ScreenSize');
figure_width = 1000;
figure_height = 900; % 2行用に高さを確保
pos_x = (screen_size(3) - figure_width) / 2;
pos_y = (screen_size(4) - figure_height) / 2;

% [left, bottom, width, height]
pos_top_left     = [0.05, 0.53, 0.43, 0.38];
pos_top_right    = [0.52, 0.53, 0.43, 0.38];
pos_bottom_left  = [0.05, 0.05, 0.43, 0.38];
pos_bottom_right = [0.52, 0.05, 0.43, 0.38];


fprintf('全スライスの 2D/3D ペア表示を生成します...\n');

%% --- 3. 全スライス ループプロット ---

% 1. 全スライスをループ処理
for slice_idx = 1:total_slices
    
    % 2. Figureを新規作成するタイミングか判定
    plot_idx_in_figure = mod(slice_idx - 1, slices_per_figure); % 0 か 1
    
    if plot_idx_in_figure == 0
        % (slice_idxが 1, 3, 5, ... の時に新しいFigureを開く)
        fig = figure; 
        set(fig, 'Position', [pos_x, pos_y, figure_width, figure_height]);
        set(fig, 'Color', 'w');
        
        figure_num = floor((slice_idx - 1) / slices_per_figure) + 1;
        
        sgtitle(['QSM 2D/3D Comparison - Figure ' num2str(figure_num) ...
                 ' (Slices ' num2str(slice_idx) ' & ' ...
                 num2str(min(slice_idx + 1, total_slices)) ')']);
    end
    
    % 3. 現在のスライスのデータを準備
    qsm_slice_data = QSM(:, :, slice_idx);
    current_mask = Mask(:, :, slice_idx);
    
    % 4. Figure内の行 (1:上段, 2:下段) を決定
    row_in_figure = plot_idx_in_figure + 1;

    
    % --- 5. 左側: 2D (imshow) 表示 ---
    if row_in_figure == 1
        ax_2d = subplot(rows, cols, 1); % subplot(2,2,1)
        set(ax_2d, 'Position', pos_top_left);
    else
        ax_2d = subplot(rows, cols, 3); % subplot(2,2,3)
        set(ax_2d, 'Position', pos_bottom_left);
    end
    
    imagesc(ax_2d, qsm_slice_data, 'AlphaData', current_mask);
    caxis(ax_2d, qsm_clim);
    axis equal tight;
    axis off;
    colormap(ax_2d, 'gray');
    title(ax_2d, ['Slice ' num2str(slice_idx) ' (2D)']);
    colorbar(ax_2d);


    % --- 6. 右側: 3D (mesh) 表示 ---
    if row_in_figure == 1
        ax_3d = subplot(rows, cols, 2); % subplot(2,2,2)
        set(ax_3d, 'Position', pos_top_right);
    else
        ax_3d = subplot(rows, cols, 4); % subplot(2,2,4)
        set(ax_3d, 'Position', pos_bottom_right);
    end
    
    % 3D mesh用にデータを準備 (マスク外をNaNに)
    qsm_slice_3d = qsm_slice_data;
    qsm_slice_3d(current_mask == 0) = NaN;

    mesh(ax_3d, qsm_slice_3d);
    axis tight;
    daspect([50 50 1]); % Z軸(高さ)を強調
    colormap(ax_3d, 'default');
    title(ax_3d, ['Slice ' num2str(slice_idx) ' (3D)']);
    xlabel(ax_3d, 'X Index'); ylabel(ax_3d, 'Y Index'); zlabel(ax_3d, 'Susceptibility [ppm]');
    colorbar(ax_3d);
    
    % 3Dプロットのカラースケールも2Dと合わせる
    caxis(ax_3d, qsm_clim);

end

fprintf('表示が完了しました。\n');

%{
2Dの図における灰色の部分は、計算された磁化率（QSM）の値が 0 [ppm] に近い、正常な脳組織を示しています。

## 解説
ご提示いただいた visualize_QSM_all_slices.m などのスクリプトでは、以下の設定を行っています。

Colormap (配色): colormap(ax, 'gray'); これは、配色を「グレースケール」（白黒の濃淡）に設定しています。

Color Scale (色の範囲): qsm_clim = [-0.5, 0.5]; caxis(ax, qsm_clim); これは、データの値を色にマッピングする範囲を指定しています。

-0.5 [ppm] 以下の値は、すべて黒で表示されます。

+0.5 [ppm] 以上の値は、すべて白で表示されます。

0 [ppm] の値は、ちょうど中間の灰色で表示されます。

なぜ灰色が「正常」なのか
QSM（磁化率）は、水（脳脊髄液, CSF）を基準 (0 [ppm]) として計算されます。 脳の大部分（灰白質や白質）の磁化率は、この基準に非常に近い、ごくわずかな値（例: -0.1 ～ +0.1 [ppm]）しか持ちません。

したがって、大部分の正常な脳組織は、0 [ppm] に近い「灰色」で表示されます。

白と黒は何か
ご提示の画像（image_6e29e1.png）に見られる真っ白な点や真っ黒な点は、 Fit_ppm_complex の計算エラーが増幅された結果、+0.5 [ppm] や -0.5 [ppm] の表示範囲をはるかに超えた**異常値（アーティファクト）**です。

白い点: +18 [ppm] といった巨大なプラスの異常値が、表示限界の +0.5 [ppm]（白）に張り付いて（飽和して）表示されています。

黒い点: 同様に、マイナスの異常値が -0.5 [ppm]（黒）に張り付いて表示されています。
%}