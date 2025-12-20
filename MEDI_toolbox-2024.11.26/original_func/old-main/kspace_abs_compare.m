%==================================================================================================
% kspceのabsを静止と体動で比較する
% 【改変版: 3段階 階段状ROI置換 + Figure 2枚分割表示】
%==================================================================================================

fprintf('スクリプトを開始します (Figure 2枚分割版)\n');
clear variables;
close all;

% --- パス設定 (スクリプトの先頭に追加) ---
% 現在のスクリプトがある場所を取得
% current_dir = fileparts(mfilename('fullpath'));
% 
% % 関数が入っているフォルダのパスを作る
% % '..' は「1つ上の階層」を意味します。
% % 構成が main と original_func が兄弟フォルダ(同じ階層)の場合:
% func_path = fullfile(current_dir, '..', 'original_func'); 
% 
% % パスに追加する
% addpath(genpath(func_path));
% % Windowsのパスを直接指定する例
addpath('C:\Users\hamaguchi\MATLB\Projects\lkb-MRI-sub\MEDI_toolbox-2024.11.26\original_func');

% 確認用表示 (必要なければ消してOK)
fprintf('関数フォルダをパスに追加しました: %s\n', func_path);

% ★★★ 計算負荷設定 ★★★
DS_FACTOR = 1; 
fprintf('ダウンサンプリング係数: %d\n', DS_FACTOR);

% %% --- 1. パラメータ設定 ---
% image_file_2DGE_0deg_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_0deg_H'; % !! 要変更 !!
% image_file_2DGE_1_2_Rotate_H = 'F:\hamaguchi\copy\20241205_RawData_H\Volunteer_Rotate_H\2DGE_1-2_Rotate_H'; % !! 要変更 !!


image_file_2DGE_0deg_H = 'C:\Users\hamaguchi\Downloads\matlab\2DGE_0deg_H'; % !! 要変更 !!
image_file_2DGE_1_2_Rotate_H = 'C:\Users\hamaguchi\Downloads\matlab\2DGE_1-2_Rotate_H'; % !! 要変更 !!
image_file_2DGE_30sec_Rotate_H = 'C:\Users\hamaguchi\Downloads\matlab\2DGE_30sec_Rotate_H';
image_file_0 = '/Users/nori/Downloads/matlab/'; % !! 要変更 !!

image_file_1 = '1_data';
image_file_3 = '3_output_data'; 
image_file_4 = '4_rolate_output_data'; 

mag_filename = '1st_2DGE_0deg_mag.raw';
phase_filename = '1st_2DGE_0deg_phase.raw';

F_mag_filename = '1st_2DGE_1_2_Rotate_mag.raw';
F_phase_filename = '1st_2DGE_1_2_Rotate_phase.raw';

S_mag_filename = '2DGE_30sec_Rotate_H_mag.raw';
S_phase_filename = '2DGE_30sec_Rotate_H_phase.raw';


% --- パラメータ ---
alpha = 3;       % 分割数
beta  = 352;     % 全列数
gamma = 100;     % 開始列
pix_start_row = 116; 
pix_start_col = gamma; 
target_slice = 20; % 比較対象スライス

params = struct();
params.original_matrix_size = [512, 512, 23];
extention = 2.0/1.3;
theta = -18.6;
rotation_axis = [0 0 1];

filename_base = sprintf('Artifact_th%.1f_StepShape_Alpha%d_Beta%d_Gamma%d', theta, alpha, beta, gamma);

load_base_path = fullfile(image_file_2DGE_0deg_H, image_file_1);
load_mask_path = fullfile(image_file_2DGE_0deg_H, image_file_3);
F_load_base_path = fullfile(image_file_2DGE_1_2_Rotate_H, image_file_1);
F_load_mask_path = fullfile(image_file_2DGE_1_2_Rotate_H, image_file_3);
S_load_base_path = fullfile(image_file_2DGE_30sec_Rotate_H, image_file_1);
S_load_mask_path = fullfile(image_file_2DGE_30sec_Rotate_H, image_file_3);
% save_path = fullfile(image_file_0, image_file_4);
% if ~exist(save_path, 'dir'), mkdir(save_path); end

% --- ハイブリッド化パラメータ ---
cutted_matrix_x = round(224 / DS_FACTOR);
cutted_matrix_y = round(beta / DS_FACTOR);
pix_start_row = round(pix_start_row / DS_FACTOR);
pix_start_col = round(pix_start_col / DS_FACTOR);
gamma = round(gamma / DS_FACTOR);


%% --- 2. データ読み込み & ダウンサンプリング ---
fprintf('\n2. データ読み込み中...\n');
dims_orig = params.original_matrix_size;
precision = 'double=>double';


% 関数呼び出し
[iMag_4D, iPhase_4D] = load_mri_data(...
    load_base_path, ...      % base_path
    mag_filename, ...        % mag_filename
    phase_filename, ...      % phase_filename
    dims_orig, ...                   % dims (必須)
    'DS_FACTOR', 1, ...              % オプション
    'MaskPath', load_mask_path, ... % オプション (別フォルダの場合)
    'Precision', precision ...       % オプション
);

original_img = iMag_4D .* exp(1i * iPhase_4D);


[iMag_4D, iPhase_4D] = load_mri_data(...
    F_load_base_path, ...      % base_path
    F_mag_filename, ...        % mag_filename
    F_phase_filename, ...      % phase_filename
    dims_orig, ...                   % dims (必須)
    'DS_FACTOR', 1, ...              % オプション
    'MaskPath', F_load_mask_path, ... % オプション (別フォルダの場合)
    'Precision', precision ...       % オプション
);

F_img = iMag_4D .* exp(1i * iPhase_4D);

[iMag_4D, iPhase_4D] = load_mri_data(...
    S_load_base_path, ...      % base_path
    S_mag_filename, ...        % mag_filename
    S_phase_filename, ...      % phase_filename
    dims_orig, ...                   % dims (必須)
    'DS_FACTOR', 1, ...              % オプション
    'MaskPath', S_load_mask_path, ... % オプション (別フォルダの場合)
    'Precision', precision ...       % オプション
);
S_img = iMag_4D .* exp(1i * iPhase_4D);

params.matrix_size = size(iMag_4D);
matrix_x = params.matrix_size(1);
matrix_y = params.matrix_size(2);
num_slices = params.matrix_size(3);


%% --- 3. 拡張 & 4. 回転準備 ---
fprintf('\n3-4. 実空間拡張とパラメータ準備...\n');
magnification = round(matrix_y * extention);
y_center_final_ext = floor(magnification / 2) + 1;
y_start_final = y_center_final_ext - floor(matrix_y / 2);
y_end_final = y_start_final + matrix_y - 1;

extend_org = complex(zeros(matrix_x, magnification, num_slices));
extend_org(:,y_start_final:y_end_final, :) = original_img;
F_extend_org = complex(zeros(matrix_x, magnification, num_slices));
F_extend_org(:,y_start_final:y_end_final, :) = F_img;
S_extend_org = complex(zeros(matrix_x, magnification, num_slices));
S_extend_org(:,y_start_final:y_end_final, :) = S_img;

% 回転角度
angles = [theta / alpha, (2 * theta) / alpha, (3 * theta) / alpha];
fprintf('回転角度設定: %.2f度, %.2f度, %.2f度\n', angles(1), angles(2), angles(3));


%% --- 7. 正解データとの比較 (画像 & k空間) ---
fprintf('\n7. 正解データとの比較画像を生成中...\n');

% Log変換用無名関数
get_log_mag = @(k) log(abs(k) + 1);

% 回転関数
rot_func_crop = @(img, ang) complex(...
    imrotate(real(img), ang, 'bilinear', 'crop'), ...
    imrotate(imag(img), ang, 'bilinear', 'crop'));

% ----------------------------------------------------
% A. リファレンス(正解)データの準備
% ----------------------------------------------------
% 実画像 (Magnitude)
F_Slice_Img = F_extend_org(:,:,target_slice);
% k空間 (Log Magnitude)
F_Slice_K_Raw = fftshift(fft2(F_extend_org(:,:,target_slice)));
F_Slice_K_Log = get_log_mag(F_Slice_K_Raw);

S_Slice_Img = S_extend_org(:,:,target_slice);
% k空間 (Log Magnitude)
S_Slice_K_Raw = fftshift(fft2(S_extend_org(:,:,target_slice)));
S_Slice_K_Log = get_log_mag(S_Slice_K_Raw);


% ----------------------------------------------------
% B. 比較対象(シミュレーション)データの生成と差分計算
% ----------------------------------------------------
Base_Complex_Slice = extend_org(:,:,target_slice); % Input (0deg)

% --- 1. 回転なし (0度) ---
Img_0_ext = rot_func_crop(Base_Complex_Slice, 0);
% 実画像差分
Diff_Img_0 = imabsdiff(abs(F_Slice_Img), abs(Img_0_ext));
% k空間差分
Img_0_K_Log = get_log_mag(fftshift(fft2(Img_0_ext)));
Diff_K_0 = abs(F_Slice_K_Log - Img_0_K_Log);

% --- 2. theta/3 回転 ---
Img_1_ext = rot_func_crop(Base_Complex_Slice, 0);
% 実画像差分
Diff_Img_1 = imabsdiff(abs(F_Slice_Img), abs(Img_1_ext));
% k空間差分
Img_1_K_Log = get_log_mag(fftshift(fft2(Img_1_ext)));
Diff_K_1 = abs(F_Slice_K_Log - Img_1_K_Log);

% --- 3. 2*theta/3 回転 ---
Img_2_ext = rot_func_crop(Base_Complex_Slice, 0);
% 実画像差分
Diff_Img_2 = imabsdiff(abs(F_Slice_Img), abs(Img_2_ext));
% k空間差分
Img_2_K_Log = get_log_mag(fftshift(fft2(Img_2_ext)));
Diff_K_2 = abs(F_Slice_K_Log - Img_2_K_Log);

% --- 4. theta 回転 ---
Img_3_ext = rot_func_crop(Base_Complex_Slice, 0);
% 実画像差分
Diff_Img_3 = imabsdiff(abs(F_Slice_Img), abs(Img_3_ext));
% k空間差分
Img_3_K_Log = get_log_mag(fftshift(fft2(Img_3_ext)));
Diff_K_3 = abs(F_Slice_K_Log - Img_3_K_Log);

% --- 5 
% 実画像差分
Diff_Img_4 = imabsdiff(abs(F_Slice_Img), abs(S_Slice_Img));
% k空間差分
Diff_K_4 = abs(F_Slice_K_Log - S_Slice_K_Log);


% ----------------------------------------------------
% C. Figure表示 (2枚に分割, 各2x2配置)
% ----------------------------------------------------
% スケーリング決定 (それぞれの最大誤差に合わせる)
max_diff_img = max([max(Diff_Img_0(:)), max(Diff_Img_1(:)), max(Diff_Img_2(:)), max(Diff_Img_4(:))]);
if max_diff_img == 0, max_diff_img = 1; end

max_diff_k = max([max(Diff_K_0(:)), max(Diff_K_1(:)), max(Diff_K_2(:)), max(Diff_K_4(:))]);
if max_diff_k == 0, max_diff_k = 1; end


% --- Figure 1: 実空間 (Image) の差分 (2行2列) ---
figure('Name', 'Comparison 1: Real Image Space', 'WindowState', 'maximized');
colormap(gray(256));

subplot(2, 2, 1); imshow(Diff_Img_0, [0, max_diff_img]);
title({'[Image Diff]', '0 deg vs Correct'}); colorbar;

subplot(2, 2, 2); imshow(Diff_Img_1, [0, max_diff_img]);
title({'[Image Diff]', '\theta/3 deg vs Correct'}); colorbar;

subplot(2, 2, 3); imshow(Diff_Img_2, [0, max_diff_img]);
title({'[Image Diff]', '2\theta/3 deg vs Correct'}); colorbar;

subplot(2, 2, 4); imshow(Diff_Img_4, [0, max_diff_img]);
title({'[Image Diff]', '\theta deg vs Correct'}); colorbar;

sgtitle(sprintf('Real Image Difference (Slice: %d)\nBlack=Match, White=Mismatch', target_slice));


% --- Figure 2: k空間 (K-Space) の差分 (2行2列) ---
figure('Name', 'Comparison 2: K-Space Log-Mag', 'WindowState', 'maximized');
colormap(gray(256));

subplot(2, 2, 1); imshow(Diff_K_0, [0, max_diff_k]);
title({'[K-Space Log Diff]', '0 deg vs Correct'}); colorbar;

subplot(2, 2, 2); imshow(Diff_K_1, [0, max_diff_k]);
title({'[K-Space Log Diff]', '\theta/3 deg vs Correct'}); colorbar;

subplot(2, 2, 3); imshow(Diff_K_2, [0, max_diff_k]);
title({'[K-Space Log Diff]', '2\theta/3 deg vs Correct'}); colorbar;

subplot(2, 2, 4); imshow(Diff_K_4, [0, max_diff_k]);
title({'[K-Space Log Diff]', '\theta deg vs Correct'}); colorbar;

sgtitle(sprintf('K-Space Log Difference (Slice: %d)\nBlack=Match, White=Mismatch', target_slice));


% %% --- ローカル関数 ---
% function save_raw_data(filepath, data)
%     fid = fopen(filepath, 'w');
%     if fid == -1, error('ファイルが開けませんでした'); end
%     fwrite(fid, data, 'double'); fclose(fid);
% end


function [img_4d, phase_4d] = load_mri_data(base_path, mag_name, phase_name, dims, varargin)
% LOAD_MRI_DATA MRIのバイナリデータと関連するMATファイルを読み込む関数
%
% [img_4d, phase_4d] = load_mri_data(base_path, mag_name, phase_name, dims)
% [img_4d, phase_4d] = load_mri_data(..., 'DS_FACTOR', 1, 'MaskPath', '...', 'Precision', 'single')
%
% 入力:
%   base_path  : バイナリファイルがあるフォルダパス
%   mag_name   : Magnitudeファイル名 (.binなど)
%   phase_name : Phaseファイル名 (.binなど)
%   dims       : 元データの次元 [Height, Width, Depth, (Vol)]
%
% オプション (Name-Valueペア):
%   'DS_FACTOR': ダウンサンプリング係数 (デフォルト: 1)
%   'MaskPath' : phase.mat 等があるパス (デフォルト: base_pathと同じ)
%   'Precision': freadの精度 (デフォルト: 'double')

    % --- オプション引数の解析 ---
    p = inputParser;
    addRequired(p, 'base_path', @ischar);
    addRequired(p, 'mag_name', @ischar);
    addRequired(p, 'phase_name', @ischar);
    addRequired(p, 'dims', @isnumeric);
    addParameter(p, 'DS_FACTOR', 1, @isnumeric);
    addParameter(p, 'MaskPath', base_path, @ischar); % 指定なければbase_pathと同じ
    addParameter(p, 'Precision', 'double', @ischar); % デフォルトはsingle(float32)

    parse(p, base_path, mag_name, phase_name, dims, varargin{:});
    
    % 変数への展開
    ds_factor = p.Results.DS_FACTOR;
    mask_path = p.Results.MaskPath;
    precision = p.Results.Precision;
    
    % --- 1. Magnitudeデータの読み込み ---
    mag_full_path = fullfile(base_path, mag_name);
    fid = fopen(mag_full_path, 'rb');
    if fid == -1
        error('Magnitudeファイルが開けません: %s', mag_full_path);
    end
    raw_mag = fread(fid, inf, precision);
    fclose(fid);
    
    % Reshape (オリジナルサイズ)
    mag_orig = reshape(raw_mag, dims);
    
    % --- 2. Phaseデータの読み込み (バイナリ) ---
    % ※元のコードでは読み込んでいますが、最終的なphase_4dにはphase.matのiFreqを使っています。
    %   ここでは一応読み込み処理を残しますが、必要なければ削除可能です。
    phase_full_path = fullfile(base_path, phase_name);
    fid = fopen(phase_full_path, 'rb');
    if fid == -1
        warning('Phaseファイルが開けません: %s', phase_full_path);
        % バイナリがない場合は何もしない（下のtry-catchに委ねる）
    else
        % 読み込むが、元のロジック通りならここは使われない可能性がある
        fread(fid, inf, precision); 
        fclose(fid);  
    end
    
   % --- 3. Phase/Freq情報 (.mat)  の読み込みと結合 ---
    % 元コードのロジック: バイナリではなく mask_path 内の phase.mat (iFreq) を優先使用
    try
        tmp = load(fullfile(mask_path, 'phase.mat'), 'iFreq');
        iFreq = tmp.iFreq;
        % PDF.mat (RDF) も読み込んでいましたが、出力に使われていないため省略します
        % 必要ならここで読み込んでください
    catch
        warning('phase.mat (iFreq) が見つかりません。ゼロで埋めます。');
        iFreq = zeros(dims);
    end
    
    % --- 4. ダウンサンプリング処理と出力 ---
    % Magnitude
    img_4d = mag_orig(1:ds_factor:end, 1:ds_factor:end, :);
    
    % Phase (iFreqを使用)
    phase_4d = iFreq(1:ds_factor:end, 1:ds_factor:end, :);

end