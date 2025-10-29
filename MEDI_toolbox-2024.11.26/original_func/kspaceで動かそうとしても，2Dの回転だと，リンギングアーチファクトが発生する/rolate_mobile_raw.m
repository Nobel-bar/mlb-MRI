%==================================================================================================
% [V2] RAWデータ (Real_0ch...) を読み込み、
% スライスごとにk空間で2Dモーションをシミュレートし、
% 元画像とアーティファクト画像を並べて表示するMATLABスクリプト

%==================================================================================================

fprintf('スクリプトを開始します (V2: k空間回転 スライスごと処理)\n');
clear variables;
close all;

%% --- 1. 初期設定 ---
fprintf('1. パラメータを設定しています...\n');

% パス設定
image_file_1 = '/Users/nori/Downloads/matlab'; % !! 要変更 !!
image_file_2 = '2_orignal_data';
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

% 入力ファイル名 (拡張子なし)
input_Re_name = 'Imgn_0ch__1_1_1_1_1_0_0_1_23_1_1_1';
input_Im_name = 'Imgn_0ch__1_1_1_1_1_0_0_1_23_1_1_1';

% --- params構造体: QSMデータの撮像パラメータ ---
params = struct();
params.voxel_size = [1.0, 1.0, 1.0];         % !! 要変更 !!
params.matrix_size = [512, 512, 23];        % !! 要変更 !!
params.CF = 123.2 * 1e6;                       % !! 要変更 !!
params.TE = [0.015];                           % !! 要変更 !!
params.B0_dir = [0, 0, 1];                     % !! 要変更 !!

% --- サイズに関するパラメータ ---
orig_matrix_x = 512; % 元データのマトリクスサイズ
orig_matrix_y = 768;
cutted_matrix_x = 224; % 実際に収集されたk空間の有効データサイズ
cutted_matrix_y = 352;
final_matrix_x = params.matrix_size(1); % 最終的に出力する画像のサイズ
final_matrix_y = params.matrix_size(2); 
extention = 2.0/1.3;
magnification = round(orig_matrix_y * extention); % Y方向の拡大後サイズ (約1182)

% --- シミュレーションパラメータ ---
theta = -18.6;         % 回転角度 (度)
cross_section = 'axial'; % ファイル名用
IdealSize = 512;       % 画像処理の基準サイズ (パディング計算用)

% --- k空間ハイブリッド化パラメータ ---
% (メモ: extend_k_space の中心基準のインデックス)
x_center_ext = floor(orig_matrix_x / 2) + 1;
x_start_ext = x_center_ext - floor(cutted_matrix_x / 2);
x_end_ext = x_start_ext + cutted_matrix_x - 1;

y_center_ext = floor(magnification / 2) + 1;
y_start_ext = y_center_ext - floor(cutted_matrix_y / 2);
y_end_ext = y_start_ext + cutted_matrix_y - 1;

% row_indices / col_indices は k-space の *絶対インデックス* を指すように変更
% 元のコードの 145:368, 81:432 は、切り出した 224x352 の中でのインデックスではなく、
% 512xMagnification の中でのインデックスである必要があります。
% ここでは、元の 512x512 を基準にしたインデックスと仮定します。
row_indices_abs = 145:368; 
col_indices_abs = y_start_ext:y_end_ext; % Y方向は切り出した領域全体と仮定

% pix_start_row は row_indices_abs の *相対インデックス* ではなく *絶対インデックス* である必要があります
% 元のコード(pix_start_row = 112)が 145 から数えて112行目という意味なら修正が必要
% ここでは 112行目から、と解釈します
pix_start_row_abs = 112+145; 
width = 112;           
% ※ハイブリッド化のインデックスはロジックに合わせて再確認してください。
%   ここでは、元のk空間の112行目から112行分、と仮定します。
hybrid_row_indices = pix_start_row_abs:(pix_start_row_abs + width - 1);


%% --- 2. RAWデータの読み込み ---
fprintf('2. RAWデータを読み込んでいます...\n');
filename_input_Re = fullfile(load_base_path, input_Re_name); % .raw を追加
filename_input_Im = fullfile(load_base_path, input_Im_name); % .raw を追加

% ベクトルとして読み込み、3D配列に変換
fileID_Re = fopen(filename_input_Re, 'r');
if fileID_Re == -1, error('ファイルが開けませんでした: %s', filename_input_Re); end
data_vector_re = fread(fileID_Re, inf, 'single');
fclose(fileID_Re);
Slice = numel(data_vector_re) / (orig_matrix_x * orig_matrix_y);
if mod(Slice, 1) ~= 0, error('実数部のファイルサイズが不正です。'); end
original_img_Re = reshape(data_vector_re, [orig_matrix_x, orig_matrix_y, Slice]);

fileID_Im = fopen(filename_input_Im, 'r');
if fileID_Im == -1, error('ファイルが開けませんでした: %s', filename_input_Im); end
data_vector_im = fread(fileID_Im, inf, 'single');
fclose(fileID_Im);
original_img_Im = reshape(data_vector_im, [orig_matrix_x, orig_matrix_y, Slice]);

orig_img_3d = complex(original_img_Re, original_img_Im);
fprintf('%d x %d x %d の画像を正常に読み込みました。\n', orig_matrix_x, orig_matrix_y, Slice);
clear original_img_Re original_img_Im data_vector_re data_vector_im;

%% --- 3. k空間への変換とP0補正 ---
fprintf('3. k空間への変換とP0補正を行っています...\n');
k_space_orig = fftshift(fftn(orig_img_3d));

[max_val, max_idx] = max(abs(k_space_orig(:)));
[kk, mm, nn] = ind2sub(size(k_space_orig), max_idx);
fprintf('k空間の最大値は座標 (%d, %d, %d) にあります。\n', kk, mm, nn);
p0_factor = k_space_orig(max_idx) / max_val;
k_space_p0 = k_space_orig / p0_factor;
clear k_space_orig;

% 拡大したk_space (回転前データ) を用意
extend_k_space = complex(zeros([orig_matrix_x magnification Slice],'double'));

y_center_org = floor(orig_matrix_y / 2) + 1;
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2);
y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1;

extend_k_space(x_start_ext:x_end_ext, y_start_ext:y_end_ext, :) = k_space_p0(x_start_ext:x_end_ext, y_start_org_cut:y_end_org_cut, :);
clear k_space_p0;

%% --- 4. スライスごとのループ処理を開始 ---
fprintf('4. スライスごとのシミュレーションループを開始します (全 %d スライス)...\n', Slice);

% パディングサイズ (3倍) を計算
padded_size = 3 * IdealSize;
% パディングされたk空間の中央（回転前k空間を配置する場所）を計算
x_start_idx_pad = IdealSize + 1;
x_end_idx_pad = 2 * IdealSize;
y_start_idx_pad = round((padded_size - magnification) / 2) + 1;
y_end_idx_pad = y_start_idx_pad + magnification - 1; % -1 が必要

% 回転計算用の座標グリッド (ループ外で1回だけ計算)
[X_in, Y_in] = meshgrid(1:padded_size, 1:padded_size);
centerX = (padded_size + 1) / 2; % 2. 回転中心を計算
centerY = (padded_size + 1) / 2;
X_out = X_in;  % 3. 出力k空間の各ピクセルが、回転前のどこにあったかを逆算
Y_out = Y_in;
X_shifted = X_out - centerX;
Y_shifted = Y_out - centerY;
theta_rad = theta * (pi/180);
cosT = cos(theta_rad);
sinT = sin(theta_rad);
X_orig = X_shifted * cosT + Y_shifted * sinT + centerX;  % 逆回転の計算
Y_orig = -X_shifted * sinT + Y_shifted * cosT + centerY;
clear X_out Y_out X_shifted Y_shifted cosT sinT theta_rad;

% 最終的な切り出し位置の計算
y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(final_matrix_y / 2);
y_end_final = y_start_final + final_matrix_y - 1;


for slice_idx = 7:Slice - 7
    
    fprintf('  --- スライス %d / %d を処理中 ---\n', slice_idx, Slice);

    % --- 4a. 回転前のk空間スライスを取得 ---
    k_space_original_slice = extend_k_space(:,:,slice_idx);

    %% --- 4b. k空間の回転（モーションのシミュレート） ---
    % 複素数のゼロ行列（k空間パディング用 2D）を作成
    padded_k_space_slice = complex(zeros(padded_size, padded_size));
    
    % k空間キャンバスの中央に「回転前」のk空間スライスを配置
    padded_k_space_slice(x_start_idx_pad:x_end_idx_pad, y_start_idx_pad:y_end_idx_pad) = k_space_original_slice;

    % ★ 'interp2' を使用してk空間を回転 (Toolbox不要) ★
    fprintf('    (Toolbox無) interp2 を使用して手動でk空間を回転中...\n');
    
    % k空間スライスを実数部と虚数部に分ける
    padded_k_Re = real(padded_k_space_slice);
    padded_k_Im = imag(padded_k_space_slice);

    % interp2 で補間を実行 (実数部と虚数部で別々に)
    rotated_k_Re = interp2(X_in, Y_in, padded_k_Re, X_orig, Y_orig, 'linear', 0);
    rotated_k_Im = interp2(X_in, Y_in, padded_k_Im, X_orig, Y_orig, 'linear', 0);

    % 再び複素数k空間データに結合
    rotated_padded_k_space_slice = rotated_k_Re + 1i * rotated_k_Im;
    fprintf('    ...手動回転が完了しました。\n');

    % --- 4c. 中央の切り出し ---
    % 中央の 512xmagnification 領域（回転後のk空間スライス）を切り出す
    k_space_rotated_slice = rotated_padded_k_space_slice(x_start_idx_pad:x_end_idx_pad, y_start_idx_pad:y_end_idx_pad);
    
    % メモリ節約
    clear padded_k_space_slice padded_k_Re padded_k_Im rotated_k_Re rotated_k_Im rotated_padded_k_space_slice;

    %% --- 5. k空間データのハイブリッド化 ---
    
    % --- 5a. k空間のハイブリッド化 ---
    % k_space_original_slice をコピー
    hybrid_k_slice = k_space_original_slice;
    
    % 指定された行 (hybrid_row_indices) を、回転後のk空間データで上書き
    hybrid_k_slice(hybrid_row_indices, :) = k_space_rotated_slice(hybrid_row_indices, :);

    %% --- 6. アーティファクト画像の再構成 ---
    
    % (メモ: full_hybrid_k ではなく hybrid_k_slice を ifft2 します)
    % 逆フーリエ変換でアーティファクト画像を生成
    artifact_img_ext = ifft2(ifftshift(hybrid_k_slice));
    
    % 最終的な 512x512 に切り出し
    artifact_img = artifact_img_ext(:, y_start_final:y_end_final);

    % (マスクは現在読み込んでいないため、適用は省略)
    
    % メモリ節約
    clear k_space_rotated_slice hybrid_k_slice artifact_img_ext;

    %% --- 7. (改良) 1スライスごとに元画像とアーティファクト画像を並べて表示 ---
    
    % 比較用の「元の画像」もk空間から再構成する
    original_img_ext = ifft2(ifftshift(k_space_original_slice));
    original_img = original_img_ext(:, y_start_final:y_end_final);
    
    % (マスク適用は省略)

    % 各スライスごとに新しいFigureを作成する
    figure('Name', sprintf('Slice %d Comparison', slice_idx));
    
    % 1行2列のグリッドの1番目（左側）
    subplot(1, 2, 1);
    imshow(abs(original_img), []);
    title(sprintf('Original (Slice %d)', slice_idx));
    
    % 1行2列のグリッドの2番目（右側）
    subplot(1, 2, 2);
    imshow(abs(artifact_img), []);
    title(sprintf('Artifact (Slice %d)', slice_idx));

    % Figureを強制的に今すぐ描画する
    drawnow;
    
    clear original_img_ext original_img; % メモリ節約

    %% --- 8. (修正) 結果のファイル保存 (ループ内に移動) ---
    
    % 最終画像を permute (転置) する (元のスクリプト 221行目)
    final_artifact_img = permute(artifact_img, [2 1]);
    
    % % ファイル名にスライス番号を含める
    % filename_base = sprintf('slice%03d_2st_rolate', slice_idx);
    % save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(final_artifact_img));
    % save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(final_artifact_img));
    % save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(final_artifact_img));
    % save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(final_artifact_img));

end % --- スライスループ (for slice_idx) の終了 ---

fprintf('...全 %d スライスのシミュレーションが完了しました。\n', Slice);
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
