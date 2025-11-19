%% --- 9. RMSEの計算 ---
fprintf('9. 3Dボリューム全体のRMSEを計算しています...\n');

% 比較する2つの3D複素数ボリューム
vol_1 = original_img;
vol_2 = artifact_img;

% 1. 絶対値（大きさ）の差のRMSEを計算
% (最も一般的な画像の差の評価方法)
fprintf('    画像の大きさ(Magnitude)に基づいたRMSEを計算中...\n');

% 大きさを計算
mag_1 = abs(vol_1);
mag_2 = abs(vol_2);

% 差の2乗
squared_error_mag = (mag_1 - mag_2) .^ 2;

% 2乗誤差の平均
mean_squared_error_mag = mean(squared_error_mag(:));

% 平均の平方根
rmse_mag = sqrt(mean_squared_error_mag);

fprintf('    -> 大きさ(Magnitude)のRMSE: %f\n', rmse_mag);

% 2. (オプション) 複素数そのもののRMSEを計算
% (位相情報も含む誤差評価)
fprintf('    複素数(Complex)に基づいたRMSEを計算中...\n');

% 複素数の差の2乗 (絶対値の2乗)
squared_error_complex = abs(vol_1 - vol_2) .^ 2;

% 2乗誤差の平均
mean_squared_error_complex = mean(squared_error_complex(:));

% 平均の平方根
rmse_complex = sqrt(mean_squared_error_complex);

fprintf('    -> 複素数(Complex)のRMSE: %f\n', rmse_complex);

fprintf('RMSEの計算が完了しました。\n');