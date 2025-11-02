%{
画像データに対して指定した領域（マスク）内の値を用いて、
 2次元多項式フィッティングを行うものです。
主に、MRIの位相画像などに見られる、緩やかな背景変動（バックグラウンドノイズ）を
除去する目的で使われることが多いです。
%}
fprintf('スクリプトを開始します (3D多項式フィッティング V2)\n');
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
load_base_path = fullfile(image_file_0, image_file_1);
load_mask_path = fullfile(image_file_0, image_file_3);
save_path = fullfile(image_file_0, image_file_5);

if ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf('保存フォルダを作成しました: %s\n', save_path);
end

% フィッティング多項式の最大次数
fit_degree = 5;

% --- params構造体: QSMデータの撮像パラメータ ---
params = struct();
params.voxel_size = [1.0, 1.0, 1.0];        % !! 要変更 !! : ボクセルサイズ
params.matrix_size = [512, 512, 23];        % !! 要変更 !! : 行列サイズ
params.CF = 123.2 * 1e6;                    % !! 要変更 !! : 中心周波数 (Hz)
params.TE = [0.015];                        % !! 要変更 !! : エコー時間 (秒)
params.B0_dir = [0, 0, 1];                  % !! 要変更 !! : 静磁場方向

%% --- 2. データの読み込み ---
fprintf('2. アンラップデータを読み込んでいます...\n')
Re_filepath = fullfile(load_base_path, '1st_2DGE_0deg_Re.raw');
Im_filepath = fullfile(load_base_path, '1st_2DGE_0deg_Im.raw');

% データの次元を定義 [x, y, z, echo]
dims = [params.matrix_size, length(params.TE)];
num_elements = prod(dims); 
precision = 'double=>double'; % 'double' (64bit) を使用

% --- Reデータの読み込み ---
fid_re = fopen(Re_filepath, 'rb');
if fid_re == -1, error('Reファイルを開けませんでした: %s', Re_filepath); end
Re_vec = fread(fid_re, inf, precision);
fclose(fid_re);
if numel(Re_vec) ~= num_elements
    error('強度ファイルのデータサイズが期待される次元と一致しません。');
end
Re_4D = reshape(Re_vec, dims);
clear Re_vec; % メモリ節約

% --- Imデータの読み込み ---
Im_phase = fopen(Im_filepath, 'rb');
if Im_phase == -1, error('位相ファイルを開けませんでした: %s', Im_filepath); end
Im_vec = fread(Im_phase, inf, precision);
fclose(Im_phase);
if numel(Im_vec) ~= num_elements
    error('Imファイルのデータサイズが期待される次元と一致しません。');
end
Im_4D = reshape(Im_vec, dims);
clear Im_vec; % メモリ節約

% --- 変数定義 ---
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);
num_echos = length(params.TE); 
echo_idx = 1; % 最初のエコーのみ使用

original_img = complex(Re_4D(:,:,:, echo_idx), Im_4D(:,:,:, echo_idx));
clear Re_4D Im_4D; % メモリ節約

% --- Mask.mat の読み込み ---
mask_path = fullfile(load_mask_path, 'Mask.mat');
if exist(mask_path, 'file')
    load(mask_path); % 'Mask' という変数がワークスペースに読み込まれます
    fprintf('Mask.mat を読み込みました。\n');
    Mask = logical(Mask);
else
    error('Mask.mat が見つかりませんでした: %s', mask_path);
end

fprintf('データの読み込み完了。%d スライス、%d エコーのデータを処理します。\n', num_slices, num_echos);

% (注) 読み込んだ iFreq のサイズを直接取得する
img_size = size(original_img);
if numel(img_size) < 3
    img_size(3) = 1; % 2Dデータの場合
end
matrix_size_x = img_size(1);
matrix_size_y = img_size(2);
Slice = img_size(3);

% (注) params.matrix_size と読み込みデータのスライス数が一致するか確認
if Slice ~= params.matrix_size(3)
    warning('params.matrix_size(3) と読み込んだスライス数が異なります。');
end
if matrix_size_x ~= params.matrix_size(1) || matrix_size_y ~= params.matrix_size(2)
     warning('params.matrix_size(1,2) と読み込んだXYサイズが異なります。');
end

fprintf('読み込んだデータサイズ: %d x %d x %d\n', matrix_size_x, matrix_size_y, Slice);

%% --- 3. 座標グリッドの構築 (全スライス共通) ---

% 座標グリッドの作成 (2D)
[Y, X] = meshgrid(1:matrix_size_y, 1:matrix_size_x);

% 画像全体の座標 (reshapeされたX, Yではなく、2Dのまま計算)
flat_X = X(:);
flat_Y = Y(:);

% 中央スライスのインデックス（メッシュ表示用）
center_slice_idx = floor(Slice / 2) + 1;

%% --- 4. スライスごとのループ処理 ---

fprintf('全 %d スライスのフィッティング処理を開始します...\n', Slice);

% (修正) fitting は実数(位相)を格納するため complex 不要
fitting = complex(zeros([matrix_size_x, matrix_size_y, Slice],'double'));

for slice_idx = 1:Slice
    
    fprintf('   --- スライス %d / %d を処理中 ---\n', slice_idx, Slice);
    
    % 現在のスライスを抽出
    current_slice_complex = original_img(:,:,slice_idx);
    
    % (スライスごとの2Dマスクを取得)
    current_mask_2D = Mask(:,:,slice_idx);
    
    % マスク内のデータ点のみを抽出
    x_data = X(current_mask_2D);
    y_data = Y(current_mask_2D);
    
    % ★★★ 修正点 1: フィッティング対象を angle() に ★★★
    % (位相データを抽出)
    z_data = current_slice_complex(current_mask_2D);

    % (もし z_data が列ベクトルでない場合は列ベクトルに変換)
    if size(z_data, 2) > 1
        z_data = z_data(:);
    end

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

    % --- 4.1. 最小二乗法による係数の求解 (スライスごと) ---
    % (Vがランク落ちしている場合などに対応するため、\ (mldivide) を使用)
    coefficients = V \ z_data;

    % --- 4.2. フィッティング曲面の再構成 (スライスごと) ---
    % ★★★ 修正点 2: フィッティング曲面の再構成コードを追加 ★★★
    
    % 画像全体 (flat_X, flat_Y) を使って設計行列 V_full を作成
    V_full = zeros(numel(flat_X), size(V, 2)); % Vと同じ列数
    col_idx = 1;
    for total_degree = 0:fit_degree
        for i = 0:total_degree
            j = total_degree - i;
            V_full(:, col_idx) = (flat_X.^i) .* (flat_Y.^j);
            col_idx = col_idx + 1;
        end
    end
    
    % 求めた係数を使って画像全体のフィッティング値を計算
    fitting_values_flat = V_full * coefficients;
    
    % 2D画像に戻す
    fitting_surface = reshape(fitting_values_flat, [matrix_size_x, matrix_size_y]);

    % --- 4.3. 結果の格納 ---
    fitting(:,:,slice_idx) = fitting_surface;
    
    % (オプション) 位相補正後の画像も保存する場合
    % corrected_phase = angle(current_slice_complex) - fitting_surface;
    % corrected_img = abs(current_slice_complex) .* exp(1i * corrected_phase);
    % ★★★ 修正点 3: save_raw_data を使用 ★★★
    % save_raw_data(fullfile(save_path, sprintf('corrected_img_slice%03d_Re.raw', slice_idx)), real(corrected_img));
    % save_raw_data(fullfile(save_path, sprintf('corrected_img_slice%03d_Im.raw', slice_idx)), imag(corrected_img));

end % --- スライスループの終了 ---

fprintf('...全 %d スライスのフィッティングが完了しました。\n', Slice);


%% --- 5. (追加) 中央スライスのメッシュ表示 ---
%{
fprintf('5. 中央スライス (%d) のメッシュ表示を作成します...\n', center_slice_idx);

try
    % 元の位相データ（マスク内）
    phase_orig = angle(original_img(:,:,center_slice_idx));
    mask_center = Mask(:,:,center_slice_idx);
    phase_orig(~mask_center) = NaN; % マスク外を非表示に

    % フィッティングされた背景位相
    fitting_center = fitting(:,:,center_slice_idx);
    
    figure;
    
    % 元のデータ
    subplot(1, 2, 1);
    mesh(Y, X, phase_orig);
    title(sprintf('Original Phase (Slice %d)', center_slice_idx));
    xlabel('Y'); ylabel('X'); zlabel('Phase (rad)');
    axis tight;
    
    % フィッティング結果
    subplot(1, 2, 2);
    mesh(Y, X, fitting_center);
    title(sprintf('Fitted Surface (Slice %d, Degree %d)', center_slice_idx, fit_degree));
    xlabel('Y'); ylabel('X'); zlabel('Phase (rad)');
    axis tight;
    
    % (オプション) 補正後のデータ
    % figure;
    % corrected_phase = phase_orig - fitting_center;
    % mesh(Y, X, corrected_phase);
    % title(sprintf('Corrected Phase (Slice %d)', center_slice_idx));
    % axis tight;

catch ME
    warning('メッシュ表示の作成に失敗しました: %s', ME.message);
end
%}


%% --- 6. P0補正と保存 ---
% (注) fitting 配列は実数（背景位相マップ）です。
% この実数データに対して P0 補正を行うコードになっています。

fprintf('6. P0補正と最終データの保存を開始します...\n');
%{
% --- k空間への変換とP0補正 ---
k_space_orig = fftshift(fftn(fitting));
[max_val, max_idx] = max(abs(k_space_orig(:)));
[kk, mm, nn] = ind2sub(size(k_space_orig), max_idx);
fprintf('k空間の最大値は座標 (%d, %d, %d) にあります。\n', kk, mm, nn);

% (注) fitting が実数のため、k_space_orig(max_idx) は実数のはず
% (数値誤差により虚部を持つ可能性があるため、そのまま複素数として扱う)
p0_factor = k_space_orig(max_idx) / max_val; % max_val は abs() なので実数
k_space_p0 = k_space_orig / p0_factor;
img = ifftn(ifftshift(k_space_p0));
%}
img = fitting;
% --- 保存 ---
% ★★★ 修正点 3: save_raw_data を使用 ★★★
filename_base = sprintf('fitting_P0_corrected');
save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(img));
save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(img));
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(img));
save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(img));

% (注) P0補正前の元々のフィッティング結果（実数）も保存する場合
% save_raw_data(fullfile(save_path, 'fitting_surface_3D.raw'), fitting);

fprintf('結果は %s に保存されました。\n', save_path);

% -------------------------------------------------------------------
% スクリプトの最後にローカル関数を定義します
% -------------------------------------------------------------------
function save_raw_data(filepath, data)
    % ★★★ (軽微) コメントを double に修正 ★★★
    % データを 'double' (元の入力と同じ) で保存
    fid = fopen(filepath, 'w');
    if fid == -1
        error('ファイルが開けませんでした: %s', filepath);
    end
    fwrite(fid, data, 'double');
    fclose(fid);
end
