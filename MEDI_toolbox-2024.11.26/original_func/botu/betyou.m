% 1. 画像データの読み込み (例: DICOMまたはNIfTI)
% Magnitude画像（3次元配列）を読み込みます
mag_img = dicomread('あなたの画像ファイル.dcm'); 
% ※または niftiread('mag.nii') など
% ※もし4次元データなら、最初の1時点目だけを取り出してください: mag_img = raw_img(:,:,:,1);

% 2. ピクセルサイズ（Voxel Size）の確認
% DICOM情報などから調べて入力します (例: 1mm x 1mm x 2mm)
voxel_size = [1, 1, 2]; 

% 3. 行列サイズ（Matrix Size）の取得
matrix_size = size(mag_img);

% 4. BETの実行 (ここが重要！)
% マスクを作成します。この BET 関数が bet2.mexw64 を呼び出します。
Mask = BET(mag_img, matrix_size, voxel_size);

% 結果の確認
montage(Mask); % マスク画像を表示して確認
