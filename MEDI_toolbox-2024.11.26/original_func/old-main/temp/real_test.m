git clean -fd%==================================================================================================
% QSM RAWデータ (3D) を読み込み、
% 3D実空間で回転 (imrotate3) をシミュレートし、
% k空間でハイブリッド化 (一部置換) を行うMATLABスクリプト
%
% 修正: 背景磁場不均一による幾何学的歪み(Readout distortion)シミュレーションを追加
%==================================================================================================

fprintf('スクリプトを開始します (3D実空間回転 + k空間ハイブリッド + Readout Distortion)\n');
clear variables;
close all;

%% --- 1. 撮像・シミュレーション パラメータ設定 ---
fprintf('1. パラメータを設定しています...\n');

% (パス設定は元のまま維持)
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
% image_file_0 = image_file_000;

% 読み込みパスと保存パスを定義
load_base_path = fullfile(image_file_0, image_file_1);
load_mask_path = fullfile(image_file_0, image_file_3);
save_path = fullfile(image_file_0, image_file_4);

% ファイル名
mag_filename = '1st_2DGE_0deg_mag.raw'; % !! 要変更 !!
phase_filename = '1st_2DGE_0deg_phase.raw'; % !! 要変更 !!

% 保存先フォルダが存在しない場合は作成する
if ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf('保存フォルダを作成しました: %s\n', save_path);
end

% --- params構造体: QSMデータの撮像パラメータ ---
params = struct();
params.matrix_size = [512, 512, 23]; % !! 要変更 !! : 行列サイズ
params.TE = [0.015]; % !! 要変更 !! : エコー時間 (秒)

% --- ★追加: Readoutシミュレーション用パラメータ ---
params.dwell_time = 2e-5; % 20マイクロ秒 (Bandwidth = 50kHz相当) ※必要に応じて調整してください

% --- 拡張と回転のパラメータ ---
extention = 2.0/1.3; % Y方向の拡張率 (約1.54倍)
magnification = round( params.matrix_size(2) * extention); % 拡張後のYサイズ (512 -> 788)
theta = -18.6; % 回転角度 (度)
rotation_axis = [0 0 1]; % 回転軸 (Z軸)
IdealSize = 512; % 画像処理の基準サイズ (現在未使用)

% --- k空間ハイブリッド化パラメータ ---
cutted_matrix_x = 224; % 実際に収集されたk空間の有効データサイズ
cutted_matrix_y = 352;
width = 112; % 置き換える行数
pix_start_row = 112; % k空間ROIの何行目から置き換えるか (1-based)

%% --- 2. RAWデータの読み込み ---
fprintf('2. RAWデータと処理済みデータを読み込んでいます...\n');

mag_filepath = fullfile(load_base_path, mag_filename);
phase_filepath = fullfile(load_base_path, phase_filename);
precision = 'double=>double'; 

% --- 強度・位相データの読み込み ---
fid_mag = fopen(mag_filepath, 'rb');
if fid_mag == -1, error('強度ファイルを開けませんでした: %s', mag_filepath); end
data_vector_mag = fread(fid_mag, inf, precision);
fclose(fid_mag);

% 修正ポイント: エコー数 (length(params.TE)) も考慮してスライス枚数を計算する
num_echoes = length(params.TE);
pixels_per_slice = params.matrix_size(1) * params.matrix_size(2);

% (総データ数) / (1枚の画素数 * エコー数) = スライス枚数
params.matrix_size(3) = numel(data_vector_mag) / (pixels_per_slice * num_echoes);

if mod(params.matrix_size(3), 1) ~= 0
    error('ファイルサイズが不正です。設定した行列サイズやエコー数がファイルと一致しません。'); 
end

% データの次元を定義 [x, y, z, echo]
dims = [params.matrix_size, num_echoes];
% ベクトルを4D配列にリシェイプ
iMag_4D = reshape(data_vector_mag, dims);
% 不要になったベクトルを削除してメモリを解放 (推奨)
clear data_vector_mag; 

% --- 位相データの読み込み ---
fid_phase = fopen(phase_filepath, 'rb');
if fid_phase == -1, error('位相ファイルを開けませんでした: %s', phase_filepath); end
% freadの結果を直接 reshape に渡してメモリ節約
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
echo_idx = 1; % 最初のエコーのみ使用

fprintf('データの読み込み完了。%d スライス、%d エコーのデータを処理します。\n', num_slices);
original_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * iPhase_4D(:,:,:, echo_idx));

Highpass_img = iMag_4D(:,:,:, echo_idx) .* exp(1i * (RDF));
Back_img =  exp(1i * (iFreq - RDF)); 

%% --- 3. 実空間データの拡張 (ゼロパディング) ---
fprintf('3. 実空間データのY方向を拡張 (ゼロパディング) しています...\n');

extend_org = complex(zeros(matrix_x, magnification, num_slices));
extend_high = complex(zeros(matrix_x, magnification, num_slices));
extend_back = complex(zeros(matrix_x, magnification, num_slices));

y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(matrix_y / 2);
y_end_final = y_start_final + matrix_y - 1;

extend_org(:,y_start_final:y_end_final, :) = original_img;
extend_high(:,y_start_final:y_end_final, :) = Highpass_img; 
extend_back(:,y_start_final:y_end_final, :) = Back_img; 

fprintf('   3D実空間の準備が完了しました。\n');

%% --- 4. 3D 実空間の回転（モーションのシミュレート） ---
fprintf('4. 3D 実空間の回転 (imrotate3) を実行中...\n');

% --- 4.1. original_img (extend_org) の回転 ---
large_org_Re = real(extend_org);
large_org_Im = imag(extend_org);
rotated_org_Re = imrotate3(large_org_Re, theta, rotation_axis, 'linear', 'crop');
rotated_org_Im = imrotate3(large_org_Im, theta, rotation_axis, 'linear', 'crop');
extend_org_rotated = complex(rotated_org_Re, rotated_org_Im);
clear large_org_Re large_org_Im rotated_org_Re rotated_org_Im;

% --- 4.2. Highpass_img (extend_high) の回転 ---
large_high_Re = real(extend_high);
large_high_Im = imag(extend_high);
rotated_high_Re = imrotate3(large_high_Re, theta, rotation_axis, 'linear', 'crop');
rotated_high_Im = imrotate3(large_high_Im, theta, rotation_axis, 'linear', 'crop');
extend_high_rotated = complex(rotated_high_Re, rotated_high_Im);
clear large_high_Re large_high_Im rotated_high_Re rotated_high_Im;

% --- 4.3. Back_img (extend_back) の回転 ---
large_back_Re = real(extend_back);
large_back_Im = imag(extend_back);
rotated_back_Re = imrotate3(large_back_Re, theta, rotation_axis, 'linear', 'crop');
rotated_back_Im = imrotate3(large_back_Im, theta, rotation_axis, 'linear', 'crop');
extend_back_rotated = complex(rotated_back_Re, rotated_back_Im);
clear large_back_Re large_back_Im rotated_back_Re rotated_back_Im;

fprintf('   3D 実空間の回転が完了しました。\n'); 


%% --- 5. 3D k空間データのハイブリッド化 (磁場不均一シミュレーション込み) ---
fprintf('5. 3D k空間のハイブリッド化とReadout歪みシミュレーションを実行中...\n');

% --- 5.1 回転なしデータ (Ground Truth) の作成 ---
fprintf('    回転前のk空間 (hybrid_space) を作成中 (fftn)...\n');
base_and_simulate_space = fft2_3d_slice_by_slice(extend_org); % 回転なし
base_and_direct_space = base_and_simulate_space; % 回転なし

% --- 5.2 回転ありデータ (Distortionなし) の作成 ---
fprintf('    回転後のk空間 (direct_space) を作成中 (fftn)...\n');
direct_space = fft2_3d_slice_by_slice(extend_org_rotated);  % 単純な回転のみ

% --- 5.3 回転ありデータ + 磁場不均一 (Physics Simulation) の作成 ---
fprintf('    物理シミュレーションによるk空間 (simulate_space) を計算中...\n');
fprintf('    (注意: 計算負荷が高いため、時間がかかります)\n');

% 出力用変数の初期化
simulate_space = complex(zeros(size(extend_high_rotated)));

% 時間ベクトルの作成: x方向をReadout方向と仮定
% k空間の中心(DC)が t=0 (エコー中心 = TE) となるように設定
Nx = size(extend_high_rotated, 1); % 512
t_vector = ((0:Nx-1) - floor(Nx/2)) * params.dwell_time; 

% 背景位相マップ (Angle) の抽出
% extend_back_rotated は exp(i * gamma * dB * TE) と仮定
% 回転による補間で振幅が変化している可能性があるため、純粋な位相のみ抽出して使用
Background_Phase_Map_At_TE = angle(extend_back);
Current_TE = params.TE(1);

% --- Readout Simulation Loop ---
% 周波数エンコードの各ステップ(時間)ごとに、磁場の影響を受けた画像を計算しFFTする
for kx_idx = 1:Nx
    if mod(kx_idx, 50) == 0, fprintf('    Line %d / %d 処理中...\n', kx_idx, Nx); end
    
    % 現在のサンプリング時刻 (TEからの相対時間)
    t = t_vector(kx_idx);
    
    % 時間比率 (t / TE)
    % Phase(t) = (Phase_at_TE / TE) * t
    time_ratio = t / Current_TE;
    
    % その時刻における誤差位相項を作成
    % exp(i * theta * ratio)
    Error_Phase_Term = exp(1i * Background_Phase_Map_At_TE * time_ratio);
    
    % 物体(Highpass)に誤差位相を乗算 (Geometric Distortionの発生源)
    % Distortion_Image = Rho(x,y) * exp(-i * gamma * dB * t)
    % ここでは extend_high_rotated が物体、Error_Phase_Term が位相誤差
    Temp_Image_At_t = extend_high_rotated .* Error_Phase_Term;
    
    % 3D FFTを実行し、k空間へ変換
    K_Temp_3D = fftshift(fftn(Temp_Image_At_t));
    
    % k空間の「現在の読み出しライン (kx_idx行目)」のみを抽出して保存
    % これにより、信号読み出し中に位相が蓄積していく効果を再現
    simulate_space(kx_idx, :, :) = K_Temp_3D(kx_idx, :, :);
end
fprintf('    物理シミュレーション完了。\n');


% --- 5.4 ハイブリッド化 (データの置換) ---
x_center = floor(matrix_x / 2) + 1; 
x_start_cut = x_center - floor(cutted_matrix_x / 2); 
y_center_org = floor(magnification / 2) + 1; 
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2); 
y_end_org_cut = y_start_org_cut + cutted_matrix_y - 1; 

% ハイブリッド化を行う行のインデックス
hybrid_row_indices = (x_start_cut + pix_start_row - 1) : (x_start_cut + pix_start_row + width - 2); 
hydrid_y_indices = y_start_org_cut: y_end_org_cut;
 
fprintf('    k空間データをハイブリッド化 (置換) しています...\n');
% 指定された行を、シミュレーション生成したk空間データで上書き
base_and_simulate_space(hybrid_row_indices, hydrid_y_indices, :) = simulate_space(hybrid_row_indices, hydrid_y_indices, :);
base_and_direct_space(hybrid_row_indices, hydrid_y_indices, :) = direct_space(hybrid_row_indices, hydrid_y_indices, :);

clear direct_space; % メモリ節約 (simulate_spaceは後で比較に使うかも？)



%% --- 6. アーティファクト画像の再構成 ---
fprintf('6. アーティファクト画像を再構成中 (ifftn)...\n');
[max_val, max_idx] = max(abs(base_and_simulate_space(:)));
p0_factor = base_and_simulate_space(max_idx) / max_val;
base_and_simulate_space_p0 = base_and_simulate_space / p0_factor;

[max_val, max_idx] = max(abs(base_and_direct_space(:)));
p0_factor = base_and_direct_space(max_idx) / max_val;
base_and_direct_space_p0 = base_and_direct_space / p0_factor;

artifact_img_ext = ifft2_3d_slice_by_slice(base_and_simulate_space_p0);
direct_img_ext = ifft2_3d_slice_by_slice(base_and_direct_space_p0);
clear base_and_simulate_space_p0 base_and_direct_space_p0; 

% 拡張したY次元 (788) を元のY次元 (512) に戻す (中央を切り出す)
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
    imshow(abs(djrect_img_slice), []);
    title(sprintf('Direct Rotation (Slice %d)', slice_idx));
     
    subplot(1, 2, 2);
    imshow(abs(artifact_img_slice), []);
    title(sprintf('With B0 Distortion (Slice %d)', slice_idx));
    drawnow;
end

% ファイル保存
filename_base = sprintf('real_2D_rotate_artifact_th%.1f', theta);
save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(artifact_img));
save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(artifact_img));

fprintf('結果は %s に保存されました。\n', save_path);

% -------------------------------------------------------------------
function save_raw_data(filepath, data)
    fid = fopen(filepath, 'w');
    if fid == -1, error('ファイルが開けませんでした: %s', filepath); end
    fwrite(fid, data, 'double'); 
    fclose(fid);
end

function k_space_3d = fft2_3d_slice_by_slice(image_3d)
% FFT2_3D_SLICE_BY_SLICE 3D画像に対してスライスごとの2D-FFTを実行する関数
%
%   [Input]
%       image_3d : (Nx, Ny, Nz) の3次元画像データ (実数または複素数)
%
%   [Output]
%       k_space_3d : (Nx, Ny, Nz) のk空間データ (複素数)
%                    各スライスに対して fftshift(fft2(...)) が適用されています。

    % 1. 入力サイズの取得
    % (MATLABでは変数名を数字で始めることはできないため、image_3dとしています)
    [rows, cols, num_slices] = size(image_3d);

    % 2. 出力用配列の事前確保 (メモリ効率化のため)
    % 結果は必ず複素数になるため、complexで初期化します
    k_space_3d = complex(zeros(rows, cols, num_slices));

    % 3. スライスごとのループ処理
    for slice_idx = 1:num_slices
        % 各スライスを取り出し、2次元フーリエ変換 -> 中心シフトr
        current_slice = image_3d(:, :, slice_idx);
        k_space_3d(:, :, slice_idx) = fftshift(fft2(current_slice));
    end
end


function image_3d = ifft2_3d_slice_by_slice(k_space_3d)
% IFFT2_3D_SLICE_BY_SLICE 3D k空間データに対してスライスごとの2D逆FFTを実行する関数
%
%   [Input]
%       k_space_3d : (Nx, Ny, Nz) のk空間データ (通常は複素数)
%
%   [Output]
%       image_3d   : (Nx, Ny, Nz) の実空間画像 (複素数)
%                    各スライスに対して ifft2(ifftshift(...)) が適用されています。

    % 1. 入力サイズの取得
    [rows, cols, num_slices] = size(k_space_3d);

    % 2. 出力用配列の事前確保
    % 逆変換結果も複素数になるため complex で初期化
    image_3d = complex(zeros(rows, cols, num_slices));

    % 3. スライスごとのループ処理
    for slice_idx = 1:num_slices
        % 各スライスを取り出す
        current_k_slice = k_space_3d(:, :, slice_idx);
        
        % ifftshift (DC成分を中心から四隅へ戻す) -> ifft2 (逆フーリエ変換)
        image_3d(:, :, slice_idx) = ifft2(ifftshift(current_k_slice));
    end
end