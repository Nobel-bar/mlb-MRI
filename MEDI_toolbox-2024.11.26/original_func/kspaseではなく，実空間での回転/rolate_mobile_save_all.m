%==================================================================================================
% QSM RAWデータ (強度/位相) を読み込み、
% 各スライスに対して2DモーションアーティファクトをシミュレートするMATLABスクリプト
%
% 概要:
% 1. 3Dの強度・位相データ (.raw) を読み込む
% 2. 3Dデータをスライスごとにループ処理
% 3. 各スライス (2D複素数画像) に対し、回転（モーション）をシミュレート
% 4. 「回転前」と「回転後」のk空間データを生成
% 5. 2つのk空間データを部分的に合成（ハイブリッド化）し、体動を模倣
% 6. 合成したk空間から画像を再構成（アーティファクト画像の生成）
% 7. 結果をスライスごとにファイル保存
%==================================================================================================

fprintf('スクリプトを開始します...\n');
clear variables;

%% --- 1. 撮像・シミュレーション パラメータ設定 ---
% このセクションは、ご自身のデータに合わせて正確に設定する必要があります。

fprintf('1. パラメータを設定しています...\n');

% !! 要変更 !! : データの基本的なパス設定
image_file_1 = '/Users/nori/Downloads/matlab/';
image_file_2 = '1_data';
image_file_3 = '3_output_data'; % 結果の保存先サブフォルダ名
image_file_4 = '4_rolate_output_data'; % 結果の保存先サブフォルダ名

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

% !! 要変更 !! : ボクセルサイズ [x, y, z] (mm)
params.voxel_size = [1.0, 1.0, 1.0];
% !! 要変更 !! : 行列サイズ [x, y, z]
params.matrix_size = [512, 512, 23];
% !! 要変更 !! : 中心周波数 (Hz) (例: 3Tスキャナ 123.2 MHz)
params.CF = 123.2 * 1e6;
% !! 要変更 !! : エコー時間 (秒単位) (マルチエコーの場合は [0.015, 0.020, ...] のようにベクトルで指定)
params.TE = [0.015];
% !! 要変更 !! : 静磁場の方向 [x, y, z] (通常は [0, 0, 1])
params.B0_dir = [0, 0, 1];

% --- シミュレーションパラメータ ---
theta = -18.6;         % 回転角度 (度)
cross_section = 'axial'; % ファイル名用
IdealSize = 512;       % 画像処理の基準サイズ (matrix_size(1) と合わせる)

% --- k空間ハイブリッド化パラメータ ---
% どのk空間領域を「回転後」データで置き換えるか
row_indices = 145:368; % k空間ROIの行 (224行)
col_indices = 81:432;  % k空間ROIの列 (352列)
width = 112;           % 置き換える行数
pix_start_row = 112;   % ROIの何行目から置き換えるか (1-based index)

% パラメータのチェック (必須フィールドの確認)
if ~isfield(params, 'voxel_size') || ...
   ~isfield(params, 'matrix_size') || ...
   ~isfield(params, 'CF') || ...
   ~isfield(params, 'TE') || ...
   ~isfield(params, 'B0_dir')
    error('params 構造体に必要なフィールドがありません。voxel_size, matrix_size, CF, TE, B0_dir を確認してください。');
end

%% --- 2. RAWデータの読み込み ---
% 強度・位相の3Dデータを読み込み、4D配列 [x, y, z, echo] に整形します。

fprintf('2. RAWデータを読み込んでいます...\n');

% !! 要変更 !! : 読み込むファイル名
mag_filename = '1st_2DGE_0deg_mag.raw';
phase_filename = '1st_2DGE_0deg_phase.raw'; % ★注意: 元のコードではmagnitude.rawになっていました。正しい位相ファイル名を指定してください。

mag_filepath = fullfile(load_base_path, mag_filename);
phase_filepath = fullfile(load_base_path, phase_filename);

% データの次元を定義 [x, y, z, echo]
dims = [params.matrix_size, length(params.TE)];
num_elements = prod(dims); % 期待される総要素数

% !! 要変更 !! : データの精度 (型)
% 'double' (64bit), 'single' (32bit), 'int16' (符号付き16bit整数) など
% 元のコードに基づき 'double' を使用します
precision = 'double=>double'; 

% --- 強度データの読み込み ---
fid_mag = fopen(mag_filepath, 'rb');
if fid_mag == -1
    error('強度ファイルを開けませんでした: %s', mag_filepath);
end
iMag_vec = fread(fid_mag, inf, precision);
fclose(fid_mag);

% 読み込んだデータサイズの検証
if numel(iMag_vec) ~= num_elements
    error('強度ファイルのデータサイズ (%d 要素) が期待される次元 [%d x %d x %d x %d = %d 要素] と一致しません。', ...
          numel(iMag_vec), dims(1), dims(2), dims(3), dims(4), num_elements);
end
% 1次元ベクトルを4D配列に整形
iMag_4D = reshape(iMag_vec, dims);
clear iMag_vec; % メモリ節約

% --- 位相データの読み込み ---
fid_phase = fopen(phase_filepath, 'rb');
 if fid_phase == -1
    error('位相ファイルを開けませんでした: %s', phase_filepath);
end
iPhase_vec = fread(fid_phase, inf, precision);
fclose(fid_phase);

% 読み込んだデータサイズの検証
if numel(iPhase_vec) ~= num_elements
    error('位相ファイルのデータサイズ (%d 要素) が期待される次元 [%d x %d x %d x %d = %d 要素] と一致しません。', ...
          numel(iPhase_vec), dims(1), dims(2), dims(3), dims(4), num_elements);
end
% 1次元ベクトルを4D配列に整形
iPhase_4D = reshape(iPhase_vec, dims);
clear iPhase_vec; % メモリ節約

% --- (オプション) Mask.mat の読み込み ---
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
num_echos = length(params.TE); % (現在は1を想定)

fprintf('データの読み込み完了。%d スライス、%d エコーのデータを処理します。\n', num_slices, num_echos);

%% --- 3. スライスごとのループ処理 ---
% 3Dデータ (iMag_4D, iPhase_4D) から1スライスずつデータを取り出し、
% 2Dのモーションシミュレーションを実行します。

% (このシミュレーションはエコーごとではなくスライスごとに行うため、
%  最初のエコー (echo_idx = 1) のみを使用します)
echo_idx = 1;

%%%%%%%%%%%%%%一旦1スライスに%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_slices = 7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('3. スライスごとのシミュレーションループを開始します (全 %d スライス)...\n', num_slices);

for slice_idx = 7:num_slices
% for slice_idx = 1:num_slices
    
    fprintf('  --- スライス %d / %d を処理中 ---\n', slice_idx, num_slices);

    % --- 3a. 現在のスライスの2D複素数データを生成 ---
    % iMag .* exp(1i * iPhase) は iMag*cos(iPhase) + 1i*iMag*sin(iPhase) と等価で高速です
    original_img = iMag_4D(:,:,slice_idx, echo_idx) .* exp(1i * iPhase_4D(:,:,slice_idx, echo_idx));
    
    
    if exist('Mask', 'var')
        original_img = original_img .* Mask(:,:,slice_idx);
    end


%% --- 4. 画像の回転（モーションのシミュレート） ---
    % imrotateによる画像の切れを防ぐため、3倍のゼロパディングキャンバスを作成
    padded_size = 3 * IdealSize;
    % 複素数のゼロ行列を作成
    padded_img = complex(zeros(padded_size, padded_size));
    
    % キャンバスの中央に現在のスライス画像(original_img)を配置
    start_idx = IdealSize + 1;
    end_idx = 2 * IdealSize;
    padded_img(start_idx:end_idx, start_idx:end_idx) = original_img;

    % ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    % 'imrotate' の代替処理 (Image Processing Toolbox が不要な方法)
    % ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
    fprintf('    (Toolbox無) interp2 を使用して手動で回転を実行中...\n');

    theta = -18.6; % 回転角度 (度)
    
    % 1. 元画像の座標グリッドを作成 (1...1536)
    [X_in, Y_in] = meshgrid(1:padded_size, 1:padded_size);
    
    % 2. 回転中心を計算 (画像のど真ん中)
    centerX = (padded_size + 1) / 2;
    centerY = (padded_size + 1) / 2;

    % 3. 出力画像の各ピクセルが、回転前のどこにあったかを逆算
    % (X, Y) も (1...1536) のグリッド
    X_out = X_in;
    Y_out = Y_in;
    
    % 座標を回転中心からの相対位置に変換
    X_shifted = X_out - centerX;
    Y_shifted = Y_out - centerY;

    % 角度をラジアンに変換
    theta_rad = theta * (pi/180);
    cosT = cos(theta_rad);
    sinT = sin(theta_rad);
    
    % 逆回転の計算: 出力(X_out, Y_out)は、入力(X_orig, Y_orig)から来た
    % (imrotate と同じ counter-clockwise のため、逆回転は +theta を使う)
    X_orig = X_shifted * cosT + Y_shifted * sinT + centerX;
    Y_orig = -X_shifted * sinT + Y_shifted * cosT + centerY;

    % 4. 複素数データを実数部と虚数部に分ける
    padded_img_Re = real(padded_img);
    padded_img_Im = imag(padded_img);

    % 5. interp2 で補間を実行
    % 'linear' は 'bilinear' と同じ線形補間
    % 'ExtrapolationValue' (範囲外の値) を 0 に設定 (imrotate の 'crop' と同じ)
    rotated_img_Re = interp2(X_in, Y_in, padded_img_Re, X_orig, Y_orig, 'linear', 0);
    rotated_img_Im = interp2(X_in, Y_in, padded_img_Im, X_orig, Y_orig, 'linear', 0);

    % 6. 再び複素数画像に結合
    rotated_padded_img = rotated_img_Re + 1i * rotated_img_Im;
    
    fprintf('    ...手動回転が完了しました。\n');
    % ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
   



% --- 4. 中央の切り出し (この後の処理は同じ) ---
    start_idx = IdealSize + 1;
    end_idx = 2 * IdealSize;

    % 中央の 512x512 領域（回転後の画像）を切り出す
    new_moved_img = rotated_padded_img(start_idx:end_idx, start_idx:end_idx);
    
    % メモリ節約のため、これ以降使わない巨大な中間変数を削除
    clear padded_img rotated_img_Re rotated_img_Im rotated_padded_img;

    %% --- 5. k空間データの生成とハイブリッド化 ---
    
    % --- 5a. k空間データの生成 ---
    % fftshiftでk空間中心(DC成分)を行列中央に配置
    k_space_original = fftshift(fft2(original_img)); % 回転前
    k_space_rotated = fftshift(fft2(new_moved_img));   % 回転後

    % --- 5b. k空間ROIの切り出し ---
    % ループを使わず、インデックス指定で一括切り出し
    cutted_original_k = k_space_original(row_indices, col_indices);
    cutted_rotated_k = k_space_rotated(row_indices, col_indices);

    % --- 5c. k空間のハイブリッド化（アーティファクト生成） ---
    % まず回転前のk空間ROIをコピー
    hybrid_cutted_k = cutted_original_k;

    % 指定された範囲(pix_start_rowからwidth行分)を回転後のk空間データで上書き
    % これが「スキャン途中で体動が起きた」ことのシミュレーションになります
    hybrid_cutted_k(pix_start_row:(pix_start_row + width - 1), :) = ...
        cutted_rotated_k(pix_start_row:(pix_start_row + width - 1), :);

    %% --- 6. アーティファクト画像の再構成 ---
    
    % 512x512のゼロ行列（k空間）を作成
    full_hybrid_k = complex(zeros(matrix_x, matrix_y));

    % ハイブリッド化したROIデータを元の位置(row_indices, col_indices)に戻す
    full_hybrid_k(row_indices, col_indices) = hybrid_cutted_k;

    % k空間の中心を四隅に戻し (ifftshift)、逆フーリエ変換 (ifft2)
    % これが最終的なアーティファクト画像です
    artifact_img = ifft2(ifftshift(full_hybrid_k));
    
    % メモリ節約
    clear k_space_original k_space_rotated cutted_original_k cutted_rotated_k hybrid_cutted_k full_hybrid_k;

    %% --- 7. 結果のファイル保存 ---
    % 現在のスライスの結果をファイルに保存
    
    pix = pix_start_row; % ファイル名用
    
    % スライス番号を含む固有のファイル名ベースを作成
    base_filename = sprintf('slice%03d_artifact_-(theta %s)-(width %s)-(start %s)_%s.raw', ...
                            slice_idx, num2str(theta), num2str(width), num2str(pix), cross_section);

    save_raw_data(fullfile(save_path, [base_filename, '_Re.raw']), real(artifact_img));
    save_raw_data(fullfile(save_path, [base_filename, '_Im.raw']), imag(artifact_img));
    save_raw_data(fullfile(save_path, [base_filename, '_mag.raw']), abs(artifact_img));
    save_raw_data(fullfile(save_path, [base_filename, '_phase.raw']), angle(artifact_img));

end % --- スライスループ (for slice_idx) の終了 ---

fprintf('...全 %d スライスのシミュレーションが完了しました。\n', num_slices);
fprintf('結果は %s に保存されました。\n', save_path);

% -------------------------------------------------------------------
% スクリプトの最後にローカル関数を定義します
% -------------------------------------------------------------------
function save_raw_data(filepath, data)
    fid = fopen(filepath, 'w');
    if fid == -1
        error('ファイルが開けませんでした: %s', filepath);
    end
    fwrite(fid, data, 'double');
    fclose(fid);
end