%==================================================================================================
% [V5] RAWデータ (Real_0ch...) を読み込み、
% 3D k空間 (fftn) に対して 3D回転 (imrotate3) を行いシミュレートする
%
% 修正点:
% 1. [最重要] 3D FFT (fftn) を使用し、3D k空間データを生成します。
% 2. Image Processing Toolbox の 'imrotate3' を使用し、3D k空間全体を
%    Z軸周りに回転させます。
% 3. 3D k空間のままハイブリッド化を行います。
% 4. ギブスリンギングを抑えるため、k空間にHamming窓を適用します。
% 5. 3D k空間全体を 'ifftn' で画像空間に再構成します。
% 6. 最後にスライスごとのループで表示・保存します。
%==================================================================================================

fprintf('スクリプトを開始します (V5: 3D-FFT/3D-Rotateアーキテクチャ + Toolbox使用)\n');
clear variables;
close all;

%% --- 1. 初期設定 ---
fprintf('1. パラメータを設定しています...\n');

% パス設定
image_file_1 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 

% 読み込みパスと保存パスを定義
load_base_path = fullfile(image_file_1, image_file_2);
load_mask_path = fullfile(image_file_1, image_file_3);
save_path = fullfile(image_file_1, image_file_4);

if ~exist(save_path, 'dir')
    mkdir(save_path);
    fprintf('保存フォルダを作成しました: %s\n', save_path);
end

% 入力ファイル名 (拡張子なし)
input_Re_name = 'Real_0ch__1_1_1_1_1_0_0_1_23_1_1_1';
input_Im_name = 'Imgn_0ch__1_1_1_1_1_0_0_1_23_1_1_1';

% --- params構造体 ---
params = struct();
params.matrix_size = [512, 512, 23];        % !! 要変更 !!

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
rotation_axis = [0 0 1]; % Z軸まわりに回転

% --- k空間データ領域のインデックス ---
x_center_ext = floor(orig_matrix_x / 2) + 1;
x_start_ext = x_center_ext - floor(cutted_matrix_x / 2); % 145
x_end_ext = x_start_ext + cutted_matrix_x - 1; % 368

y_center_ext = floor(magnification / 2) + 1;
y_start_ext = y_center_ext - floor(cutted_matrix_y / 2); % 416
y_end_ext = y_start_ext + cutted_matrix_y - 1; % 767

% k空間のハイブリッド化を行う「絶対インデックス」
hybrid_row_indices = 112:(112 + 112 - 1); % 112行目から112行分
hybrid_slice = 12;


%% --- 2. RAWデータの読み込み ---
fprintf('2. RAWデータを読み込んでいます...\n');
filename_input_Re = fullfile(load_base_path, [input_Re_name, '.raw']);
filename_input_Im = fullfile(load_base_path, [input_Im_name, '.raw']);

% 実数部
fileID_Re = fopen(filename_input_Re, 'r');
if fileID_Re == -1, error('ファイルが開けませんでした: %s', filename_input_Re); end
data_vector_re = fread(fileID_Re, inf, 'single');
fclose(fileID_Re);
Slice = numel(data_vector_re) / (orig_matrix_x * orig_matrix_y);
if mod(Slice, 1) ~= 0, error('実数部のファイルサイズが不正です。'); end
original_img_Re = reshape(data_vector_re, [orig_matrix_x, orig_matrix_y, Slice]);

% 虚数部
fileID_Im = fopen(filename_input_Im, 'r');
if fileID_Im == -1, error('ファイルが開けませんでした: %s', filename_input_Im); end
data_vector_im = fread(fileID_Im, inf, 'single');
fclose(fileID_Im);
original_img_Im = reshape(data_vector_im, [orig_matrix_x, orig_matrix_y, Slice]);

% 3D複素数画像（画像空間）
orig_img_3d = complex(original_img_Re, original_img_Im);
fprintf('%d x %d x %d の画像を正常に読み込みました。\n', orig_matrix_x, orig_matrix_y, Slice);
clear original_img_Re original_img_Im data_vector_re data_vector_im;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- 4. 3D k空間の回転（モーションのシミュレート） ---
% [修正] 'imrotate3' を使用して 3D k空間全体を回転させます
fprintf('4. 3D original空間の回転 (imrotate3) を実行中...\n');

% imrotate3 は複素数データを直接扱えないため、実数部と虚数部に分けます
original_Re = real(orig_img_3d);
original_Im = imag(orig_img_3d);

% 'linear' は 'bilinear' (2D) や 'trilinear' (3D) に相当
% 'crop' で、はみ出た部分を0にします
fprintf('    実数部を回転中...\n');
rotated_k_Re = imrotate3(original_Re, theta, rotation_axis, 'linear', 'crop');
fprintf('    虚数部を回転中...\n');
rotated_k_Im = imrotate3(original_Im, theta, rotation_axis, 'linear', 'crop');

% 再び複素数k空間データに結合
org_rotated = complex(rotated_k_Re, rotated_k_Im); % これが「回転後」
% の3D original空間だが，サンプリングされていないところを含んでいる．
smp_org_rotaled = complex(zeros([orig_matrix_x, orig_matrix_y, Slice],'double'));
smp_org_rotaled(x_start_ext : x_end_ext, y_start_ext : y_end_ext , :) = org_rotated(x_start_ext : x_end_ext, y_start_ext : y_end_ext , :);
clear original_Re original_Im rotated_k_Re rotated_k_Im org_rotated;
fprintf('    3D original 空間の回転が完了しました。\n');

%% --- 5. 3D k空間データのハイブリッド化 ---
fprintf('5. 3D original 空間のハイブリッド化を実行中...\n');
% k_space_original をコピー
hybrid_org = orig_img_3d;
    
% 指定された行 (hybrid_row_indices) を、回転後のk空間データで上書き
% Y方向 (:,) と Z方向 (Slice) にわたって置き換えます
hybrid_org(hybrid_row_indices, :, hybrid_slice) = smp_org_rotaled(hybrid_row_indices, :, hybrid_slice);

clear smp_org_rotaled; % メモリ節約

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('3. 3D k空間への変換 (fftn) とP0補正を行っています...\n');
k_space_hybrid = fftshift(fftn(hybrid_org));
k_space_org= fftshift(fftn(orig_img_3d));
clear orig_img_3d hybrid_org; % メモリ節約

[max_val, max_idx] = max(abs(k_space_hybrid(:)));
p0_factor = k_space_hybrid(max_idx) / max_val;
k_space_p0 = k_space_hybrid / p0_factor;
[max_2_val, max_2_idx] = max(abs(k_space_org(:)));
p0_factor_org = k_space_org(max_2_idx) / max_2_val;
k_space_org_p0 = k_space_org / p0_factor_org;
clear k_space_hybrid k_space_org;

% 拡大したk_spaceを用意している．
extend_k_space = complex(zeros([orig_matrix_x magnification Slice],'double'));
extend_k_space_org = complex(zeros([orig_matrix_x magnification Slice],'double'));
x_center = floor(orig_matrix_x / 2) + 1;
x_start_cut = x_center - floor(cutted_matrix_x / 2);
x_end_cut = x_start_cut+cutted_matrix_x - 1;

y_center = floor(magnification / 2) + 1;
y_start_cut = y_center - floor(cutted_matrix_y / 2);
y_end_cut = y_start_cut+cutted_matrix_y - 1;

y_center_org = floor(orig_matrix_y / 2) + 1;
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2);
y_end_org_cut = y_start_org_cut+cutted_matrix_y - 1;
extend_k_space(x_start_cut:x_end_cut, y_start_cut:y_end_cut, :) = k_space_p0(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut, :);
extend_k_space_org(x_start_cut:x_end_cut, y_start_cut:y_end_cut, :) = k_space_org_p0(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut, :);
clear k_space_org_p0 k_space_p0;

% --- 6. 最終画像の再構成と保存 ---
% 最後に一度だけ逆フーリエ変換を行い、実空間画像に戻します。
extend_img_shifted = ifftn(ifftshift(extend_k_space));
extend_img_shifted_org = ifftn(ifftshift(extend_k_space_org));
clear extend_k_space extend_k_space_org;
% imshow(real(extend_img_shifted(:,:,1)));

% 論文などでよく見られる画像の向きに合わせる (XとYを転置)
y_start_final = y_center - floor(final_matrix_y / 2);
y_end_final = y_start_final+final_matrix_y - 1;

final_img_shifted = extend_img_shifted(:, y_start_final:y_end_final, :);
final_img_shifted_org = extend_img_shifted_org(:, y_start_final:y_end_final, :);
final_img = permute(final_img_shifted, [2 1 3]);
final_img_org = permute(final_img_shifted_org, [2 1 3]);
clear final_img_shifted final_img_shifted_org;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- 7. スライスごとの表示と保存 ---
fprintf('7. スライスごとの表示と保存を開始します...\n');

for slice_idx = 7:Slice - 7
    
    fprintf('  --- スライス %d / %d を処理中 ---\n', slice_idx, Slice);

    % 現在のスライスを3Dボリュームから抽出
    original_img_slice = final_img_org(:,:,slice_idx);
    artifact_img_slice = final_img(:,:,slice_idx);

    % (マスク適用は省略)

    % 各スライスごとに新しいFigureを作成する
    figure('Name', sprintf('Slice %d Comparison', slice_idx));
    
    % 1行2列のグリッドの1番目（左側）
    subplot(1, 2, 1);
    imshow(abs(original_img_slice), []);
    title(sprintf('Original (Slice %d)', slice_idx));
    
    % 1行2列のグリッドの2番目（右側）
    subplot(1, 2, 2);
    imshow(abs(artifact_img_slice), []);
    title(sprintf('Artifact (Slice %d)', slice_idx));

    % Figureを強制的に今すぐ描画する
    drawnow;

end % --- スライスループ (for slice_idx) の終了 ---


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % ファイルを保存
% save_path = fullfile('.', image_file_1, image_file_2);
% filename_base = sprintf('1st_2DGE_0deg');
% % filename_base = sprintf('1st_%dx%d_to_%dx%d', cutted_matrix_x, cutted_matrix_y, final_matrix_x, final_matrix_y);
% save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(final_img));
% save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(final_img));
% save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(final_img));
% save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(final_img));
% 
% disp('処理が完了しました。');
% 
% % -------------------------------------------------------------------
% % スクリプトの最後にローカル関数を定義します
% % -------------------------------------------------------------------
% function save_raw_data(filepath, data)
%     fid = fopen(filepath, 'w');
%     if fid == -1
%         error('ファイルが開けませんでした: %s', filepath);
%     end
%     fwrite(fid, data, 'double');
%     fclose(fid);
% end


