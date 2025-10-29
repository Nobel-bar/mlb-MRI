%{
画像データに対して指定した領域（マスク）内の値を用いて、
 2次元多項式フィッティングを行うものです。
主に、MRIの位相画像などに見られる、緩やかな背景変動（バックグラウンドノイズ）を
除去する目的で使われることが多いです。

[修正点]
1. 3Dデータを読み込み、スライスごとにループ処理を行う構造に変更。
2. 位相画像 (angle(img)) をフィッティング対象のデータ (z_data) とするように修正。
3. 構文エラー、未定義変数、古い変数名を修正。
4. 結果の保存をループ内で行い、スライス番号を含むファイル名に変更。
5. 3Dメッシュ表示は、代表として中央スライスのみ行うように変更。
%}

fprintf('スクリプトを開始します (3D多項式フィッティング)\n');
clear variables;
close all;

%% --- 1. 初期設定 ---
fprintf('1. パラメータを設定しています...\n');

% パス設定
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

image_file_0 = image_file_00; % slab用

% 読み込みパスと保存パスを定義
load_base_path = fullfile(image_file_0, image_file_3);
load_mask_path = fullfile(image_file_0, image_file_3);
save_path = fullfile(image_file_0, image_file_5);

if ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf('保存フォルダを作成しました: %s\n', save_path);
end

% 画像サイズ
matrix_size = [512, 512];

% フィッティング多項式の最大次数
fit_degree = 5;
 
%% --- 2. データの読み込み ---
fprintf('2. アンラップデータを読み込んでいます...\n');
% (修正) phase.matから読み込まれる変数を 'data' 構造体に格納
data = load(fullfile(load_base_path, 'phase.mat'));

% (修正) iFreq変数が存在するか確認
if ~isfield(data, 'iFreq')
    error('phase.mat に iFreq という名前の変数が含まれていません。');
end
% (修正) iFreq をワークスペースに展開 (変数名がiFreqの場合)
iFreq = data.iFreq;
clear data; % メモリ節約

% マスクの読み込み
load(fullfile(load_mask_path, 'Mask.mat')); % 'Mask' という変数名で読み込まれると仮定
Mask = logical(Mask);
fprintf('データの読み込みが完了しました。\n');

%% --- 3. 座標グリッドと設計行列Vの構築 (全スライス共通) ---

% ★★★ ここから修正 ★★★
% (修正) 読み込んだ iFreq のサイズを直接取得する
img_size = size(iFreq);
if numel(img_size) < 3
    img_size(3) = 1; % 2Dデータの場合
end
matrix_size_x = img_size(1);
matrix_size_y = img_size(2);
Slice = img_size(3);

fprintf('読み込んだデータサイズ: %d x %d x %d\n', matrix_size_x, matrix_size_y, Slice);

% 座標グリッドの作成 (2D)
% (修正) params.iFreq を iFreq のサイズ (img_size) に変更
[Y, X] = meshgrid(1:matrix_size_y, 1:matrix_size_x);

% 画像全体の座標 (reshapeされたX, Yではなく、2Dのまま計算)
flat_X = X(:);
flat_Y = Y(:);
% 中央スライスのインデックス（メッシュ表示用）
center_slice_idx = floor(Slice / 2) + 1;

%% --- (修正) セクション4〜6をスライスごとのループ処理に変更 ---

fprintf('全 %d スライスのフィッティング処理を開始します...\n', Slice);
% (修正) params.iFreq を iFreq のサイズ (img_size) に変更
fitting = complex(zeros([matrix_size_x, matrix_size_y, Slice],'double'));
for slice_idx = 1:Slice
    
    fprintf('  --- スライス %d / %d を処理中 ---\n', slice_idx, Slice);
    
    % 現在のスライスを抽出
    current_slice_phase = iFreq(:,:,slice_idx);
    
    % ★★★ 修正点 1: Vの計算をループ内に移動 ★★★
    % (スライスごとの2Dマスクを取得)
    current_mask_2D = Mask(:,:,slice_idx);
    
    % マスク内のデータ点のみを抽出
    x_data = X(current_mask_2D);
    y_data = Y(current_mask_2D);
    
    % (フィッティング対象のデータも現在のマスクで抽出)
    % ★★★ 修正点 2: インデックスを Slice -> slice_idx に ★★★
    % z_data = current_slice_phase(Mask(:,:,Slice)); % <-- バグ
    z_data = current_slice_phase(current_mask_2D); % <-- 修正後

    % 設計行列Vの構築 (スライスごと)
    col_idx = 1;
    V = []; 
    for total_degree = 0:fit_degree
        for i = 0:total_degree
            j = total_degree - i;
            V(:, col_idx) = (x_data.^i) .* (y_data.^j);
            col_idx = col_idx + 1;
        end
    end
    % (Vの計算ここまで)


    % --- 4. 最小二乗法による係数の求解 (スライスごと) ---
    coefficients = V \ z_data;

    % --- 5. フィッティング曲面の再構成 (スライスごと) ---
    fitting_surface = zeros([matrix_size_x, matrix_size_y]);
    col_idx = 1;

    for total_degree = 0:fit_degree
        for i = 0:total_degree
            j = total_degree - i;
            % 全てのピクセルに対して、対応する基底と係数を掛けて足し合わせる
            term_surface = coefficients(col_idx) * (flat_X.^i) .* (flat_Y.^j);
            fitting_surface = fitting_surface + reshape(term_surface, matrix_size);
            col_idx = col_idx + 1;
        end
    end

    % --- 6. 結果の表示と保存 (スライスごと) ---
    
    % (修正) 代表スライス（例：中央）のみメッシュ表示する
    if slice_idx == center_slice_idx
        figure;
        meshc(Y, X, fitting_surface);
        title(sprintf('Fitted Surface (Slice %d, Degree: %d)', slice_idx, fit_degree));
        xlabel('Y-axis');
        ylabel('X-axis');
        zlabel('Fitted Phase');
        colorbar;
        drawnow;
    end
    fitting(:,:,slice_idx) = fitting_surface;
    % (修正) ファイル名にスライス番号を含めて保存
    % filename_base = sprintf('fitting_surface_slice%03d.raw', slice_idx);
    % save_raw_data(fullfile(save_path, filename_base), fitting_surface);
    
    % (オプション) 位相補正後の画像も保存する場合
    % corrected_phase = current_slice_phase - fitting_surface;
    % corrected_img = abs(current_slice_complex) .* exp(1i * corrected_phase);
    % save_raw_data(fullfile(save_path, sprintf('corrected_img_slice%03d_Re.raw', slice_idx)), real(corrected_img));
    % save_raw_data(fullfile(save_path, sprintf('corrected_img_slice%03d_Im.raw', slice_idx)), imag(corrected_img));

end % --- スライスループの終了 ---

save(fullfile(save_path, 'fitting.mat'), 'fitting');


fprintf('...全 %d スライスのフィッティングが完了しました。\n', Slice);
fprintf('結果は %s に保存されました。\n', save_path);

% -------------------------------------------------------------------
% スクリプトの最後にローカル関数を定義します
% -------------------------------------------------------------------
function save_raw_data(filepath, data)
    % データを 'double' ではなく 'single' (元の入力と同じ) で保存
    fid = fopen(filepath, 'w');
    if fid == -1
        error('ファイルが開けませんでした: %s', filepath);
    end
    fwrite(fid, data, 'single');
    fclose(fid);
end
