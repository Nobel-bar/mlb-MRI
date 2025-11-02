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
image_file_5 = '5_fitting_output_data'; 

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
extention = 2.0/1.3;
magnification = round( params.matrix_size(2) * extention);
% --- シミュレーションパラメータ ---
theta = -18.6;         % 回転角度 (度)
cross_section = 'axial'; % ファイル名用
IdealSize = 512;       % 画像処理の基準サイズ

% --- k空間ハイブリッド化パラメータ ---
row_indices = 145:368; % k空間ROIの行 (224行)
col_indices = 81:432;  % k空間ROIの列 (352列)
width = 112;           % 置き換える行数
pix_start_row = 112;   % ROIの何行目から置き換えるか

cutted_matrix_x = 224; % 実際に収集されたk空間の有効データサイズ
cutted_matrix_y = 352;

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
load(fullfile(load_mask_path, 'phase.mat'));
load(fullfile(load_mask_path, 'PDF.mat'));
mask_path = fullfile(load_mask_path, 'Mask.mat');

% --- 変数定義 ---
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);
num_echos = length(params.TE); 
echo_idx = 1; % 最初のエコーのみ使用

fprintf('データの読み込み完了。%d スライス、%d エコーのデータを処理します。\n', num_slices, num_echos);
original_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * iPhase_4D(:,:,:, echo_idx));

% Highpass_img = original_img ./ fitting;  % fittingを正規化したうえで
Highpass_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * (RDF));
Back_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * (iFreq - RDF));


extend_org = complex(zeros(params.matrix_size(1), magnification, (params.matrix_size(3))));


y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(matrix_y / 2);
y_end_final = y_start_final + matrix_y - 1;

extend_org(:,y_start_final:y_end_final, :) = original_img;
extend_high(:,y_start_final:y_end_final, :) = Highpass_img;
extend_back(:,y_start_final:y_end_final, :) = Back_img;


% 拡大した (回転前データ) を用意 もとの大きさの5倍の大きさのキャンバスを用意してください


fprintf('    3D空間の準備が完了しました。\n');

% k空間のハイブリッド化を行う「絶対インデックス」
hybrid_row_indices = x_start_ext+112:(x_start_ext + 112 + 112 - 1); % 112行目から112行分


%% --- 4. 3D 空間の回転（モーションのシミュレート） ---
% [修正] 'imrotate3' を使用して 3D 空間全体を回転させます
fprintf('4. 3D k空間の回転 (imrotate3) を実行中...\n');

% imrotate3 は複素数データを直接扱えないため、実数部と虚数部に分けます
large_org_Re = real(large_org);
large_org_Im = imag(large_org);

% 'linear' は 'bilinear' (2D) や 'trilinear' (3D) に相当
% 'crop' で、はみ出た部分を0にします
fprintf('    実数部を回転中...\n');
rotated_org_Re = imrotate3(large_org_Re, theta, rotation_axis, 'linear', 'crop');
fprintf('    虚数部を回転中...\n');
rotated_org_Im = imrotate3(large_org_Im, theta, rotation_axis, 'linear', 'crop');

% 再び複素数k空間データに結合
extend_org_rotated = complex(rotated_org_Re, rotated_org_Im); % これが「回転後」の3D k空間


% highpass
% imrotate3 は複素数データを直接扱えないため、実数部と虚数部に分けます
large_high_Re = real(large_high);
large_high_Im = imag(large_high);

% 'linear' は 'bilinear' (2D) や 'trilinear' (3D) に相当
% 'crop' で、はみ出た部分を0にします
fprintf('    実数部を回転中...\n');
rotated_high_Re = imrotate3(large_high_Re, theta, rotation_axis, 'linear', 'crop');
fprintf('    虚数部を回転中...\n');
rotated_high_Im = imrotate3(large_high_Im, theta, rotation_axis, 'linear', 'crop');

% 再び複素数k空間データに結合
extend_high_rotated = complex(rotated_high_Re, rotated_high_Im); % これが「回転後」の3D k空間

% background
% imrotate3 は複素数データを直接扱えないため、実数部と虚数部に分けます
large_back_Re = real(large_back);
large_back_Im = imag(large_back);

% 'linear' は 'bilinear' (2D) や 'trilinear' (3D) に相当
% 'crop' で、はみ出た部分を0にします
fprintf('    実数部を回転中...\n');
rotated_back_Re = imrotate3(large_back_Re, theta, rotation_axis, 'linear', 'crop');
fprintf('    虚数部を回転中...\n');
rotated_back_Im = imrotate3(large_back_Im, theta, rotation_axis, 'linear', 'crop');

% 再び複素数k空間データに結合
extend_back_rotated = complex(rotated_back_Re, rotated_back_Im); % これが「回転後」の3D空間
fprintf('    3D空間の回転が完了しました。\n'); 


space_rotated = extend_high_rotated .* exp(1i * (iFreq - RDF));
space_rotated_base space_rotatedからこの大きさで切り出したい(params.matrix_size(1), magnification, (params.matrix_size(3)))));

%% --- 5. 3D k空間データのハイブリッド化 ---
fprintf('5. 3D 空間のハイブリッド化を実行中...\n');
k_space_rolate = fftshift(fftn(space_rotated));
hybrid_space = fftshift(fftn(large_org));
    
% 指定された行 (hybrid_row_indices) を、回転後のk空間データで上書き
% Y方向 (:,) と Z方向 (:,) の全体にわたって置き換えます
hybrid_space(hybrid_row_indices, :, :) = k_space_rolate(hybrid_row_indices, :, :);

x_center = floor(params.matrix_size(1) / 2) + 1;
x_start_cut = x_center - floor(cutted_matrix_x / 2);
x_end_cut = x_start_cut+cutted_matrix_x - 1;

y_center_org = floor(orig_matrix_y / 2) + 1;
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2);
y_end_org_cut = y_start_org_cut+cutted_matrix_y - 1;
final_k_space = complex(zeros(params.matrix_size(1), magnification, (params.matrix_size(3))));
final_k_space(x_start_cut:x_end_cut, y_start_cut:y_end_cut, :) = hybrid_space(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut, :);
   

% [修正] 3D k空間全体を ifftn で画像空間に戻します
fprintf('    アーティファクト画像を再構成中 (ifftn)...\n');
artifact_img_ext = ifftn(ifftshift(final_k_space));

artifact_img = artifact_img_ext(:, y_start_final:y_end_final, :);



%% --- 7. スライスごとの表示と保存 ---
fprintf('7. スライスごとの表示と保存を開始します...\n');

for slice_idx = 12:12
    
    fprintf('  --- スライス %d / %d を処理中 ---\n', slice_idx, Slice);

    % 現在のスライスを3Dボリュームから抽出
    original_img_slice = iMag_4D(:,:,slice_idx, echo_idx);
    artifact_img_slice = artifact_img(:,:,slice_idx);

    % (マスク適用は省略)

    % 各スライスごとに新しいFigureを作成する
    figure('Name', sprintf('Slice %d Comparison', slice_idx));
    
    % 1行2列のグリッドの1番目（左側）
    subplot(1, 2, 1);
    imshow(original_img_slice, []);
    title(sprintf('Original (Slice %d)', slice_idx));
    
    % 1行2列のグリッドの2番目（右側）
    subplot(1, 2, 2);
    imshow(abs(artifact_img_slice), []);
    title(sprintf('Artifact (Slice %d)', slice_idx));

    % Figureを強制的に今すぐ描画する
    drawnow;
    
    %% --- 8. (修正) 結果のファイル保存 (ループ内に移動) ---
    
    % 最終画像を permute (転置) する
    final_artifact_img_slice = permute(artifact_img_slice, [2 1]);
    
    % % ファイル名にスライス番号を含める
    % filename_base = sprintf('slice%03d_3D_rotate_artifact', slice_idx);
    % save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(final_artifact_img_slice));
    % save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(final_artifact_img_slice));
    % save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(final_artifact_img_slice));
    % save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(final_artifact_img_slice));

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


fprintf('...全 %d スライスのシミュレーションが完了しました。\n', num_slices);
fprintf('結果は %s に保存されました。\n', save_path);

