%================================================================
% PDF (背景磁場除去) および 多項式フィッティング結果 比較用プログラム
%　フィッティングの次元で比較している
% 概要:
%   PDF適用前の磁場マップ (iFreq) と
%   PDF適用後の局所磁場マップ (RDF) 、
%   さらに 2x2 から 8x8 までの多項式フィッティング結果 (homo_RDF) を
%   2D (imshow) と 3D (mesh) で並べて比較します。
%
% 依存ファイル:
%   - phase.mat (iFreq を含む)
%   - PDF.mat (RDF を含む)
%   - Mask.mat
%   - 2x2_fitting.mat, 3x3_fitting.mat, ..., 8x8_fitting.mat
%
% 使い方:
%   1. 必要な .mat ファイルがすべて揃っていることを確認します。
%   2. このスクリプトを実行します。
%================================================================
clear variables;
close all;

%% --- 1. 初期設定と共通データ読み込み ---
fprintf('1. パラメータを設定し、共通データを読み込んでいます...\n');

% パス設定
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
image_file_0 = image_file_00;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%変更あり%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load_path = fullfile(image_file_0, image_file_3);
load_fitting_path = fullfile(image_file_0, image_file_5);

try
    load(fullfile(load_path, 'phase.mat'), 'iFreq');
    load(fullfile(load_path, 'PDF.mat'), 'RDF');
    load(fullfile(load_path, 'Mask.mat'), 'Mask');
catch ME
    fprintf('共通データ (iFreq, RDF, Mask) の読み込みに失敗しました。\n');
    fprintf('パスを確認してください: %s\n', load_path);
    rethrow(ME);
end

fprintf('共通データの読み込みが完了しました。\n');

%% --- 2. 多項式フィッティング結果の読み込みと計算 ---
% 比較するフィッティングの次数
fit_degrees = 2:8; % 2x2 から 8x8 まで
num_fits = length(fit_degrees);

% 結果を格納するセル配列
% 1番目: iFreq, 2番目: RDF, 3番目以降: homo_RDF (2x2〜8x8)
total_plots = 2 + num_fits; % 合計9プロット
data_to_show = cell(1, total_plots);
titles = cell(1, total_plots);

% 1. マスク適用済みの iFreq と RDF を格納
fprintf('マスクを適用しています...\n');
data_to_show{1} = iFreq .* Mask;
titles{1} = 'Before PDF (iFreq)';

data_to_show{2} = RDF .* Mask;
titles{2} = 'After PDF (RDF)';

% 2. 各フィッティング結果を読み込み、homo_RDF を計算して格納
fprintf('多項式フィッティング結果 (%d〜%d) を読み込んでいます...\n', fit_degrees(1), fit_degrees(end));
for i = 1:num_fits
    deg = fit_degrees(i);
    filename_fitting = sprintf('%dx%d_fitting.mat', deg, deg);
    full_fit_path = fullfile(load_fitting_path, filename_fitting);
    
    if ~exist(full_fit_path, 'file')
        warning('フィッティングファイルが見つかりません: %s', full_fit_path);
        % データがない場合は NaN で埋めた配列を格納
        data_to_show{i + 2} = nan(size(iFreq));
        titles{i + 2} = sprintf('PolyFit (%dx%d) - 読込失敗', deg, deg);
    else
        fprintf('  %s を読み込み中...\n', filename_fitting);
        data_fit = load(full_fit_path, 'fitting');
        
        % iFreq から fitting 結果を引き、マスクを適用
        homo_RDF = (iFreq - data_fit.fitting) .* Mask;
        
        data_to_show{i + 2} = homo_RDF;
        titles{i + 2} = sprintf('PolyFit (%dx%d)', deg, deg);
    end
end
fprintf('全データの準備が完了しました。\n');


%% --- 3. スライス番号の指定 ---
total_slices = size(iFreq, 3);
slice_to_display = round(total_slices / 2); % 中央スライス

if slice_to_display < 1 || slice_to_display > total_slices
    fprintf('エラー: スライス番号 %d は無効です。\n', slice_to_display);
    return;
end

%% --- 4. imshow (2D) での比較 (3x3 グリッド) ---

fprintf('スライス %d の 2D (imshow) 比較を表示します。\n', slice_to_display);

figure('Name', '2D Comparison (All Fits)', 'WindowState', 'maximized');
sgtitle(sprintf('2D PDF & PolyFit Comparison - Slice %d', slice_to_display), 'FontWeight', 'bold');

for i = 1:total_plots
    ax = subplot(3, 3, i);
    
    % スライスデータを取得
    img_slice = data_to_show{i}(:, :, slice_to_display);
    
    % NaN（読み込み失敗時）でないかチェック
    if all(isnan(img_slice(:)))
        imshow(zeros(size(img_slice)), []); % 真っ黒な画像を表示
        title(titles{i});
        xlabel('データなし');
        continue;
    end
    
    imshow(img_slice, []);
    colormap(ax, 'gray');
    axis on;
    daspect([1,1,1]);
    title(titles{i});
    xlabel('X Index');
    if mod(i-1, 3) == 0 % 各行の左端のみYラベルを表示
        ylabel('Y Index');
    end
    colorbar;
end

%% --- 5. mesh (3D) での比較 (3x3 グリッド) ---

fprintf('スライス %d の 3D (mesh) 比較を表示します。\n', slice_to_display);

figure('Name', '3D Mesh Comparison (All Fits)', 'WindowState', 'maximized');
sgtitle(sprintf('3D Mesh Comparison - Slice %d', slice_to_display), 'FontWeight', 'bold');

for i = 1:total_plots
    ax = subplot(3, 3, i);
    
    % スライスデータを取得
    mesh_slice = data_to_show{i}(:, :, slice_to_display);

    % NaN（読み込み失敗時）でないかチェック
    if all(isnan(mesh_slice(:)))
        title(titles{i});
        xlabel('データなし');
        continue;
    end
    
    % 0の値を非表示に (NaNにする)
    mesh_slice(mesh_slice == 0) = NaN; 
    
    if all(isnan(mesh_slice(:)))
        % マスク内が全て0だった場合
        title(titles{i});
        xlabel('マスク内データなし');
        continue;
    end

    mesh(ax, mesh_slice);
    axis tight;
    daspect([1,1,1/50]); % Z軸のスケールを調整
    axis on;
    colormap(ax, 'default');
    xlabel('X Index');
    if mod(i-1, 3) == 0 % 各行の左端のみYラベルを表示
        ylabel('Y Index');
    end
    zlabel('Field (a.u.)');
    title(titles{i});
    colorbar;
end

fprintf('比較表示が完了しました。\n');

