%==================================================================================================
% QSM RAWデータ (3D) を読み込み、
% 3D実空間で回転 (ユーザー提供のパディング+interp2ロジック) をシミュレートし、
% k空間でハイブリッド化 (一部置換) を行うMATLABスクリプト
%
% ※ Toolbox不要、位置ズレ修正版
%==================================================================================================
fprintf('スクリプトを開始します (3D実空間回転 + k空間ハイブリッド - 位置ズレ修正版)\n');
clear variables;
close all;

%% --- 1. 撮像・シミュレーション パラメータ設定 ---
fprintf('1. パラメータを設定しています...\n');

% パス設定
image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_000 = "C:\Users\hamaguchi\Downloads\matlab\2DGE_0deg_H";
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

image_file_0 = image_file_00;
image_file_0 = image_file_000;

% 読み込みパスと保存パスを定義
load_base_path = fullfile(image_file_0, image_file_1);
load_mask_path = fullfile(image_file_0, image_file_3);
save_path = fullfile(image_file_0, image_file_4);

% 保存先フォルダが存在しない場合は作成する
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% --- params構造体: QSMデータの撮像パラメータ ---
params = struct();
params.matrix_size = [512, 512, 23]; % !! 要変更 !! : 行列サイズ
params.TE = [0.015]; % !! 要変更 !! : エコー時間 (秒)

% --- 拡張と回転のパラメータ ---
extention = 2.0/1.3; % Y方向の拡張率 (約1.54倍)
magnification = round( params.matrix_size(2) * extention); % 拡張後のYサイズ (512 -> 788)
theta = -18.6; % 回転角度 (度)

% --- k空間ハイブリッド化パラメータ ---
cutted_matrix_x = 224; % 実際に収集されたk空間の有効データサイズ
cutted_matrix_y = 352;
width = 112; % 置き換える行数
pix_start_row = 112; % k空間ROIの何行目から置き換えるか (1-based)

%% --- 2. RAWデータの読み込み ---
fprintf('2. RAWデータと処理済みデータを読み込んでいます...\n');

% ファイル名
mag_filename = '1st_2DGE_0deg_mag.raw'; % !! 要変更 !!
phase_filename = '1st_2DGE_0deg_phase.raw'; % !! 要変更 !!

mag_filepath = fullfile(load_base_path, mag_filename);
phase_filepath = fullfile(load_base_path, phase_filename);

dims = [params.matrix_size, length(params.TE)];
precision = 'double=>double'; 

% --- 強度・位相データの読み込み ---
fid_mag = fopen(mag_filepath, 'rb');
if fid_mag == -1, error('強度ファイルを開けませんでした: %s', mag_filepath); end
iMag_4D = reshape(fread(fid_mag, inf, precision), dims);
fclose(fid_mag);

fid_phase = fopen(phase_filepath, 'rb');
if fid_phase == -1, error('位相ファイルを開けませんでした: %s', phase_filepath); end
iPhase_4D = reshape(fread(fid_phase, inf, precision), dims);
fclose(fid_phase);

% --- 処理済みデータ (iFreq, RDF) の読み込み ---
try
    load(fullfile(load_mask_path, 'phase.mat'), 'iFreq');
    load(fullfile(load_mask_path, 'PDF.mat'), 'RDF');
catch ME
    fprintf('phase.mat または PDF.mat の読み込みに失敗しました。\n');
    rethrow(ME);
end

% --- 変数定義 ---
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);
echo_idx = 1; 

fprintf('データの読み込み完了。%d スライス、%d エコーのデータを処理します。\n', num_slices);

original_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * iPhase_4D(:,:,:, echo_idx));
Highpass_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * (RDF));
Back_img =  exp(1i * (iFreq - RDF)); 

%% --- 3. 実空間データの拡張 (ゼロパディング) ---
fprintf('3. 実空間データのY方向を拡張 (ゼロパディング) しています...\n');

% 拡張後のサイズの複素数ゼロ行列を初期化 (512 x 788 x 23)
extend_org = complex(zeros(matrix_x, magnification, num_slices));
extend_high = complex(zeros(matrix_x, magnification, num_slices));
extend_back = complex(zeros(matrix_x, magnification, num_slices));

% 拡張後Yサイズの中心を計算
y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(matrix_y / 2);
y_end_final = y_start_final + matrix_y - 1;

extend_org(:,y_start_final:y_end_final, :) = original_img;
extend_high(:,y_start_final:y_end_final, :) = Highpass_img; 
extend_back(:,y_start_final:y_end_final, :) = Back_img; 

fprintf('   3D実空間の準備が完了しました。\n');

%% --- 4. 3D 実空間の回転（モーションのシミュレート） ---
% ユーザー提供のロジック（3倍パディング -> interp2 -> 切り出し）を使用
fprintf('4. 3D 実空間の回転 (パディング + interp2) を実行中...\n');

% ローカル関数 perform_padded_rotation を呼び出して回転
% (各変数の回転処理を共通化)
fprintf('     (original) を回転中...\n');
extend_org_rotated = perform_padded_rotation(extend_org, theta);

fprintf('     (Highpass) を回転中...\n');
extend_high_rotated = perform_padded_rotation(extend_high, theta);

fprintf('     (Background) を回転中...\n');
extend_back_rotated = perform_padded_rotation(extend_back, theta);

fprintf('   3D 実空間の回転が完了しました。\n'); 

% --- 4.4. 回転後実空間データの合成 ---
extend_rotated = extend_high_rotated .* extend_back_rotated;


%% --- 5. 3D k空間データのハイブリッド化 ---
fprintf('5. 3D k空間のハイブリッド化を実行中...\n');

fprintf('    回転前のk空間 (base) を作成中...\n');
base_and_simulate_space = fftshift(fftn(extend_org)); 
base_and_direct_space = base_and_simulate_space;

fprintf('    回転後のk空間 (rotated) を作成中...\n');
direct_space = fftshift(fftn(extend_org_rotated)); 
simulate_space = fftshift(fftn(extend_rotated)); 

x_center = floor(matrix_x / 2) + 1; 
x_start_cut = x_center - floor(cutted_matrix_x / 2); 
x_end_cut = x_start_cut + cutted_matrix_x - 1; 

y_center_org = floor(magnification / 2) + 1; 
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2); 
y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1; 

hybrid_row_indices = x_start_cut: x_end_cut; 
hydrid_y_indices = y_start_org_cut: y_end_org_cut;
hydrid_z_indices = 13: 23;
 
fprintf('    k空間データをハイブリッド化 (置換) しています...\n');
% temp_space = complex(zeros(size(base_and_simulate_space)));
% 
% % 2. 指定した範囲（インデックス）のデータだけを、元の場所と同じ位置にコピー
% temp_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices) = ...
% base_and_simulate_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices);
% 
% % 3. 元の変数を、この「範囲外がゼロになったデータ」で上書き
% base_and_simulate_space = temp_space;
base_and_simulate_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices) = 0;
[max_val, max_idx] = max(abs(base_and_simulate_space(:)));
p0_factor = base_and_simulate_space(max_idx) / max_val;
base_and_simulate_space = base_and_simulate_space / p0_factor;

kk = complex(zeros(matrix_x, magnification, num_slices)); 
kk(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices) = simulate_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices);
[max_val, max_idx] = max(abs(kk(:)));
p0_factor = kk(max_idx) / max_val;
kk = kk / p0_factor;

base_and_simulate_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices) = kk(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices);
base_and_direct_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices) = direct_space(hybrid_row_indices, hydrid_y_indices, hydrid_z_indices);

clear simulate_space direct_space; 

%% --- 6. アーティファクト画像の再構成 ---
fprintf('6. アーティファクト画像を再構成中...\n');

[max_val, max_idx] = max(abs(base_and_simulate_space(:)));
p0_factor = base_and_simulate_space(max_idx) / max_val;
base_and_simulate_space_p0 = base_and_simulate_space / p0_factor;

[max_val, max_idx] = max(abs(base_and_direct_space(:)));
p0_factor = base_and_direct_space(max_idx) / max_val;
base_and_direct_space_p0 = base_and_direct_space / p0_factor;

artifact_img_ext = ifftn(ifftshift(base_and_simulate_space_p0));
direct_img_ext = ifftn(ifftshift(base_and_direct_space_p0));

clear base_and_simulate_space_p0 base_and_direct_space_p0; 

% 中央切り出し
artifact_img = artifact_img_ext(:, y_start_final:y_end_final, :);
direct_img = direct_img_ext(:, y_start_final:y_end_final, :);

clear artifact_img_ext direct_img_ext; 

%% --- 7. スライスごとの表示と保存 ---
fprintf('7. スライスごとの表示と保存を開始します...\n');

for slice_idx = 12:12 
    
    fprintf('   --- スライス %d / %d を処理中 ---\n', slice_idx, num_slices);
    djrect_img_slice = permute(direct_img(:,:,slice_idx), [2 1]);
    artifact_img_slice = permute(artifact_img(:,:,slice_idx), [2 1]);
    
    figure('Name', sprintf('Slice %d Comparison', slice_idx), 'WindowState', 'maximized');
     
    subplot(1, 2, 1);
    imagesc(abs(djrect_img_slice));
    colormap gray; axis image; axis off;
    title(sprintf('Original (Slice %d)', slice_idx));
     
    subplot(1, 2, 2);
    imagesc(abs(artifact_img_slice));
    colormap gray; axis image; axis off;
    title(sprintf('Artifact (Slice %d)', slice_idx));
    
    drawnow;
end 

% ファイル保存
filename_base = sprintf(['h_artifact_th%.1f'], theta);
save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(artifact_img));

fprintf('...スライス %d のシミュレーションが完了しました。\n', slice_idx);
fprintf('結果は %s に保存されました。\n', save_path);


% ===================================================================
% ローカル関数定義
% ===================================================================

function save_raw_data(filepath, data)
    fid = fopen(filepath, 'w');
    if fid == -1, error('ファイルが開けませんでした: %s', filepath); end
    fwrite(fid, data, 'double');
    fclose(fid);
end

function vol_out = perform_padded_rotation(vol_in, theta)
% PERFORM_PADDED_ROTATION
% ユーザー指定の「3倍パディング -> interp2で回転 -> 中央切り出し」ロジックを
% 3Dボリュームの各スライスに適用する関数
%
%   入力: vol_in (MxNxZ の複素数データ), theta (回転角度)
%   出力: vol_out (入力と同じサイズの回転済みデータ)

    [rows, cols, num_slices] = size(vol_in);
    vol_out = complex(zeros(rows, cols, num_slices));

    % --- 1. 回転用座標グリッドの作成 (スライス共通) ---
    % 元のコードの "padded_size = 3 * IdealSize" に相当する処理
    % 矩形画像にも対応するため、それぞれの次元を3倍にします
    pad_rows_total = rows * 3;
    pad_cols_total = cols * 3;
    
    % 座標グリッド作成 (meshgrid は [x, y] = [col, row] 順)
    [X_in, Y_in] = meshgrid(1:pad_cols_total, 1:pad_rows_total);
    
    % 回転中心 (パディング画像のど真ん中)
    centerX = (pad_cols_total + 1) / 2;
    centerY = (pad_rows_total + 1) / 2;
    
    % 座標変換の準備 (出力座標 -> 入力座標 の逆変換)
    X_out = X_in;
    Y_out = Y_in;
    
    X_shifted = X_out - centerX;
    Y_shifted = Y_out - centerY;
    
    theta_rad = theta * (pi/180);
    cosT = cos(theta_rad);
    sinT = sin(theta_rad);
    
    % 逆回転座標の計算
    % (imrotate互換のため: X(col)にcos, Y(row)にsin を適用する際の符号に注意)
    % ユーザーコード: X_orig = X_shifted * cosT + Y_shifted * sinT + centerX;
    %               Y_orig = -X_shifted * sinT + Y_shifted * cosT + centerY;
    X_orig = X_shifted * cosT + Y_shifted * sinT + centerX;
    Y_orig = -X_shifted * sinT + Y_shifted * cosT + centerY;
    
    % 画像を埋め込む位置 (パディング空間の中央)
    start_r = floor((pad_rows_total - rows)/2) + 1;
    end_r   = start_r + rows - 1;
    start_c = floor((pad_cols_total - cols)/2) + 1;
    end_c   = start_c + cols - 1;

    % --- 2. スライスごとの処理ループ ---
    for z = 1:num_slices
        slice_data = vol_in(:,:,z);
        
        % --- パディング (キャンバスへ配置) ---
        padded_img = complex(zeros(pad_rows_total, pad_cols_total));
        padded_img(start_r:end_r, start_c:end_c) = slice_data;
        
        % 実部・虚部に分離
        padded_img_Re = real(padded_img);
        padded_img_Im = imag(padded_img);
        
        % --- interp2 による回転 (補間) ---
        % X_in, Y_in はグリッド、V は値、X_orig, Y_orig は参照先座標
        rotated_img_Re = interp2(X_in, Y_in, padded_img_Re, X_orig, Y_orig, 'linear', 0);
        rotated_img_Im = interp2(X_in, Y_in, padded_img_Im, X_orig, Y_orig, 'linear', 0);
        
        rotated_padded_img = complex(rotated_img_Re, rotated_img_Im);
        
        % --- 中央切り出し ---
        vol_out(:,:,z) = rotated_padded_img(start_r:end_r, start_c:end_c);
    end
end