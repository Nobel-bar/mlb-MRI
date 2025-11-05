clear variables;
close all;

image_file_00 = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!
image_file_000 = "C:\Users\hamaguchi\Downloads\matlab";
image_file_1 = '1_data';
image_file_2 = '2_original_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 
image_file_5 = '5_fitting_output_data'; 

image_file_0 = image_file_000; % local用


% 読み込みパスを定義
load_base_path = fullfile(image_file_0, image_file_2);

% 入力ファイル名 (拡張子なし)
input_Re_name = 'Real_0ch__1_1_1_1_1_0_0_1_23_1_1_1';
input_Im_name = 'Imgn_0ch__1_1_1_1_1_0_0_1_23_1_1_1';

% サイズに関するパラメータ
orig_matrix_x = 512; % 元データのマトリクスサイズ
orig_matrix_y = 768;
cutted_matrix_x = 224; % 実際に収集されたk空間の有効データサイズ
cutted_matrix_y = 352;
final_matrix_x = 512; % 最終的に出力する画像のサイズ
final_matrix_y = 512; 
extention = 2.0/1.3;
magnification = round(orig_matrix_y * extention);

% --- 2. データの読み込み ---
filename_input_Re = fullfile('.', image_file_1, image_file_2, [input_Re_name]);
filename_input_Im = fullfile('.', image_file_1, image_file_2, [input_Im_name]);

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

orig_img = complex(original_img_Re, original_img_Im);
fprintf('%d x %d x %d の画像を正常に読み込みました。\n', orig_matrix_x, orig_matrix_y, Slice);

% --- 3. k空間への変換とP0補正 ---
k_space_orig = fftshift(fftn(orig_img));
[max_val, max_idx] = max(abs(k_space_orig(:)));
[kk, mm, nn] = ind2sub(size(k_space_orig), max_idx);
fprintf('k空間の最大値は座標 (%d, %d, %d) にあります。\n', kk, mm, nn);
p0_factor = k_space_orig(max_idx) / max_val;
k_space_p0 = k_space_orig / p0_factor;


% 拡大したk_spaceを用意している．
extend_k_space = complex(zeros([orig_matrix_x magnification Slice],'double'));
x_center = floor(orig_matrix_x / 2) + 1;
x_start_cut = x_center - floor(cutted_matrix_x / 2);
x_end_cut = x_start_cut+cutted_matrix_x - 1;

y_center_org = floor(orig_matrix_y / 2) + 1;
y_start_org_cut = y_center_org - floor(cutted_matrix_y / 2);
y_end_org_cut = y_start_org_cut+cutted_matrix_y - 1;
extend_k_space(x_start_cut:x_end_cut, y_start_cut:y_end_cut, :) = k_space_p0(x_start_cut:x_end_cut, y_start_org_cut:y_end_org_cut, :);

% --- 6. 最終画像の再構成と保存 ---
% 最後に一度だけ逆フーリエ変換を行い、実空間画像に戻します。
extend_img_shifted = ifftn(ifftshift(extend_k_space));

% imshow(real(extend_img_shifted(:,:,1)));

% 論文などでよく見られる画像の向きに合わせる (XとYを転置)
y_start_final = y_center - floor(final_matrix_y / 2);
y_end_final = y_start_final+final_matrix_y - 1;

final_img_shifted = extend_img_shifted(:, y_start_final:y_end_final, :);
final_img = permute(final_img_shifted, [2 1 3]);

%imshow(abs(final_img(:,:,1)));
%imshow((real(k_space_p0(:,:,1))));

% ファイルを保存
save_path = fullfile('.', image_file_1, image_file_2);
filename_base = sprintf('1st_2DGE_0deg');
% filename_base = sprintf('1st_%dx%d_to_%dx%d', cutted_matrix_x, cutted_matrix_y, final_matrix_x, final_matrix_y);
save_raw_data(fullfile(save_path, [filename_base, '_Re.raw']), real(final_img));
save_raw_data(fullfile(save_path, [filename_base, '_Im.raw']), imag(final_img));
save_raw_data(fullfile(save_path, [filename_base, '_mag.raw']), abs(final_img));
save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), angle(final_img));

disp('処理が完了しました。');

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

