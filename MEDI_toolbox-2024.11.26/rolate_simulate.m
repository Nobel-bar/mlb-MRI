%==================================================================================================
% QSM RAWデータ (強度/位相) を読み込み、
% 各スライスに対してk空間で2Dモーションをシミュレートし、
% 元画像とアーティファクト画像を並べて表示するMATLABスクリプト
%
% ★ 修正: 回転操作を画像空間からk空間に変更しました ★
%==================================================================================================

fprintf('スクリプトを開始します...\n');
clear variables;

%% --- 1. 撮像・シミュレーション パラメータ設定 ---
fprintf('1. パラメータを設定しています...\n');

% パス設定
image_file_1 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_2 = '1_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 

% 読み込みパスと保存パスを定義
load_base_path = fullfile(image_file_1, image_file_2);
load_mask_path = fullfile(image_file_1, image_file_3);
save_path = fullfile(image_file_1, image_file_4);

% 保存先フォルダが存在しない場合は作成する
if ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf('保存フォルダを作成しました: %s\n', save_path);
end

% --- params構造体: QSMデータの撮像パラメータ ---
params = struct();
params.voxel_size = [1.0, 1.0, 1.0];         % !! 要変更 !! : ボクセルサイズ
params.matrix_size = [512, 512, 23];        % !! 要変更 !! : 行列サイズ
params.CF = 123.2 * 1e6;                       % !! 要変更 !! : 中心周波数 (Hz)
params.TE = [0.015];                           % !! 要変更 !! : エコー時間 (秒)
params.B0_dir = [0, 0, 1];                     % !! 要変更 !! : 静磁場方向

% --- シミュレーションパラメータ ---
theta = -18.6;         % 回転角度 (度)
cross_section = 'axial'; % ファイル名用
IdealSize = 512;       % 画像処理の基準サイズ

% --- k空間ハイブリッド化パラメータ ---
row_indices = 145:368; % k空間ROIの行 (224行)
col_indices = 81:432;  % k空間ROIの列 (352列)
width = 112;           % 置き換える行数
pix_start_row = 112;   % ROIの何行目から置き換えるか

% パラメータのチェック
if ~isfield(params, 'voxel_size') || ...
   ~isfield(params, 'matrix_size') || ...
   ~isfield(params, 'CF') || ...
   ~isfield(params, 'TE') || ...
   ~isfield(params, 'B0_dir')
    error('params 構造体に必要なフィールドがありません。');
end

%% --- 2. RAWデータの読み込み ---
fprintf('2. RAWデータを読み込んでいます...\n');

% ファイル名
mag_filename = '1st_2DGE_0deg_mag.raw';     % !! 要変更 !!
phase_filename = '1st_2DGE_0deg_phase.raw'; % !! 要変更 !!

mag_filepath = fullfile(load_base_path, mag_filename);
phase_filepath = fullfile(load_base_path, phase_filename);

% データの次元を定義 [x, y, z, echo]
dims = [params.matrix_size, length(params.TE)];
num_elements = prod(dims); 
precision = 'double=>double'; % 'double' (64bit) を使用

% --- 強度データの読み込み ---
fid_mag = fopen(mag_filepath, 'rb');
if fid_mag == -1, error('強度ファイルを開けませんでした: %s', mag_filepath); end
iMag_vec = fread(fid_mag, inf, precision);
fclose(fid_mag);
if numel(iMag_vec) ~= num_elements
    error('強度ファイルのデータサイズが期待される次元と一致しません。');
end
iMag_4D = reshape(iMag_vec, dims);
clear iMag_vec; % メモリ節約

% --- 位相データの読み込み ---
fid_phase = fopen(phase_filepath, 'rb');
if fid_phase == -1, error('位相ファイルを開けませんでした: %s', phase_filepath); end
iPhase_vec = fread(fid_phase, inf, precision);
fclose(fid_phase);
if numel(iPhase_vec) ~= num_elements
    error('位相ファイルのデータサイズが期待される次元と一致しません。');
end
iPhase_4D = reshape(iPhase_vec, dims);
clear iPhase_vec; % メモリ節約

% --- Mask.mat の読み込み ---
mask_path = fullfile(load_mask_path, 'Mask.mat');
if exist(mask_path, 'file')
    load(mask_path); % 'Mask' という変数がワークスペースに読み込まれます
    fprintf('Mask.mat を読み込みました。\n');
else
    warning('Mask.mat が見つかりませんでした: %s', mask_path);
end

% --- 変数定義 ---
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);
num_echos = length(params.TE); 
echo_idx = 1; % 最初のエコーのみ使用

fprintf('データの読み込み完了。%d スライス、%d エコーのデータを処理します。\n', num_slices, num_echos);
%% --- 3. スライスごとのループ処理 ---
fprintf('3. スライスごとのシミュレーションループを開始します (全 %d スライス)...\n', num_slices);

for slice_idx = 7:num_slices-7
    
    fprintf('  --- スライス %d / %d を処理中 ---\n', slice_idx, num_slices);

    % --- 3a. 現在のスライスの2D複素数データを生成 ---
    original_img = iMag_4D(:,:,slice_idx, echo_idx) .* exp(1i * iPhase_4D(:,:,slice_idx, echo_idx));
    
    % --- 3b. 回転前のk空間を生成 ---
    % (シミュレーションの「回転前」のデータ)
    k_space_original = fftshift(fft2(original_img)); % 回転前
    
    %% --- 4. k空間の回転（モーションのシミュレート） ---
    % フーリエ変換の回転定理に基づき、k空間を直接回転させます。
    % 補間によるアーティファクトを防ぐため、画像空間と同様にゼロパディングを行います。
    

    padded_size = 3 * IdealSize;
    % 複素数のゼロ行列（k空間）を作成
    padded_k_space = complex(zeros(padded_size, padded_size));
    
    % k空間キャンバスの中央に「回転前」のk空間データを配置
    start_idx_pad = IdealSize + 1;
    end_idx_pad = 2 * IdealSize;
    padded_k_space(start_idx_pad:end_idx_pad, start_idx_pad:end_idx_pad) = k_space_original;

    % ★ 'interp2' を使用してk空間を回転 (Toolbox不要) ★
    fprintf('    (Toolbox無) interp2 を使用して手動でk空間を回転中...\n');
    
    % 1. 元k空間の座標グリッドを作成 (1...1536)
    [X_in, Y_in] = meshgrid(1:padded_size, 1:padded_size);
    
    % 2. 回転中心を計算
    centerX = (padded_size + 1) / 2;
    centerY = (padded_size + 1) / 2;

    % 3. 出力k空間の各ピクセルが、回転前のどこにあったかを逆算
    X_out = X_in;
    Y_out = Y_in;
    X_shifted = X_out - centerX;
    Y_shifted = Y_out - centerY;

    theta_rad = theta * (pi/180);
    cosT = cos(theta_rad);
    sinT = sin(theta_rad);
    
    % 逆回転の計算
    X_orig = X_shifted * cosT + Y_shifted * sinT + centerX;
    Y_orig = -X_shifted * sinT + Y_shifted * cosT + centerY;

    % 4. k空間データを実数部と虚数部に分ける
    padded_k_Re = real(padded_k_space);
    padded_k_Im = imag(padded_k_space);

    % 5. interp2 で補間を実行 (実数部と虚数部で別々に)
    rotated_k_Re = interp2(X_in, Y_in, padded_k_Re, X_orig, Y_orig, 'linear', 0);
    rotated_k_Im = interp2(X_in, Y_in, padded_k_Im, X_orig, Y_orig, 'linear', 0);

    % 6. 再び複素数k空間データに結合
    rotated_padded_k_space = rotated_k_Re + 1i * rotated_k_Im;
    
    fprintf('    ...手動回転が完了しました。\n');

    % --- 4b. 中央の切り出し ---
    % 中央の 512x512 領域（回転後のk空間）を切り出す
    % (シミュレーションの「回転後」のデータ)
    k_space_rotated = rotated_padded_k_space(start_idx_pad:end_idx_pad, start_idx_pad:end_idx_pad);
    
    % メモリ節約
    clear padded_k_space rotated_k_Re rotated_k_Im rotated_padded_k_space;
    clear X_in Y_in X_out Y_out X_shifted Y_shifted X_orig Y_orig;

    %% --- 5. k空間データのハイブリッド化 ---
    
    % --- 5a. k空間ROIの切り出し ---
    % k_space_original と k_space_rotated は既に準備完了
    cutted_original_k = k_space_original(row_indices, col_indices);
    cutted_rotated_k = k_space_rotated(row_indices, col_indices);

    % --- 5b. k空間のハイブリッド化 ---
    hybrid_cutted_k = cutted_original_k;
    hybrid_cutted_k(pix_start_row:(pix_start_row + width - 1), :) = ...
        cutted_rotated_k(pix_start_row:(pix_start_row + width - 1), :);

    %% --- 6. アーティファクト画像の再構成 ---
    
    % 512x512のk空間にハイブリッドデータを戻す
    full_hybrid_k = complex(zeros(matrix_x, matrix_y));
    full_hybrid_k(row_indices, col_indices) = hybrid_cutted_k;
    
    % 逆フーリエ変換でアーティファクト画像を生成
    artifact_img = ifft2(ifftshift(full_hybrid_k));
    
    % ★ マスクの適用 ★
    % 最終画像 (artifact_img) に対してマスクを適用し、背景ノイズを除去します
    if exist('Mask', 'var')
        artifact_img = artifact_img .* Mask(:,:,slice_idx);
    else
        warning('変数 "Mask" が見つかりません。マスク処理をスキップします。');
    end

    % メモリ節約
    clear k_space_original k_space_rotated cutted_original_k cutted_rotated_k hybrid_cutted_k full_hybrid_k;

    %% --- 7. (改良) 1スライスごとに元画像とアーティファクト画像を並べて表示 ---
    
    % 表示用に、元の画像にもマスクを適用する
    if exist('Mask', 'var')
        original_img_masked = original_img .* Mask(:,:,slice_idx);
    else
        original_img_masked = original_img; % マスクがない場合はそのまま
    end

    % 各スライスごとに新しいFigureを作成する
    figure('Name', sprintf('Slice %d Comparison', slice_idx));
    
    % 1行2列のグリッドの1番目（左側）
    subplot(1, 2, 1);
    imshow(abs(original_img_masked), []);
    title(sprintf('Original (Slice %d)', slice_idx));
    
    % 1行2列のグリッドの2番目（右側）
    subplot(1, 2, 2);
    imshow(abs(artifact_img), []);
    title(sprintf('Artifact (Slice %d)', slice_idx));

    % Figureを強制的に今すぐ描画する
    drawnow;



end % --- スライスループ (for slice_idx) の終了 ---

fprintf('...全 %d スライスのシミュレーションが完了しました。\n', num_slices);
fprintf('結果は %s に保存されました。\n', save_path);

