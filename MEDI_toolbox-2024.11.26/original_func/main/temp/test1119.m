%==========================================================================
% MRI信号収集シミュレーション (Matlab) - 最終修正版
% 背景磁場不均一性を大幅に増幅し、歪みを可視化
%==========================================================================

% ... (1. パラメータ設定は前回と同じ) ...
params = struct();
params.matrix_size = [512, 512]; 
params.TE = 0.005; % エコー時間 [s]
params.gamma = 2 * pi * 42.57e6; % プロトンの磁気回転比 [rad/(s*T)]
params.FOV = 0.256; % 視野 (Field of View) [m]
N_x = params.matrix_size(1);
N_y = params.matrix_size(2);
dx = params.FOV / N_x; % ボクセルサイズ [m]
% G_read (傾斜磁場強度) の計算値は約 81.8 mT/m

% ... (2. 入力データと不均一磁場マップの準備) ...

% A. スピン密度 (rho) は前回と同じ円形ファントムを使用

% B. 背景磁場不均一性マップ (dB_map) の再設定
dB_map = zeros(N_x, N_y);
x_coords_m = ((1:N_x) - N_x/2) * dx; % x空間座標 [m]

% ★★★ 歪み可視化のため、dB_mapの最大値を約 400 uT に設定 ★★★
% FOVの端で 400 uT のオフセットが発生する線形勾配を設定
dB_max_target = 400e-6; % 400 uT
dB_slope_T_per_m_FINAL = dB_max_target / (params.FOV / 2);

dB_linear_x = dB_slope_T_per_m_FINAL * x_coords_m;
dB_map = repmat(dB_linear_x, N_y, 1); % x軸方向にのみ不均一性が存在するマップ [T]

fprintf('新しい最大 dB: %.1f uT\n', max(max(dB_map)) * 1e6);
fprintf('最大周波数オフセット: %.1f Hz\n', params.gamma * max(max(dB_map)) / (2*pi));

% --- 3. k空間信号の計算 (修正後のロジックを使用) ---

% 静的磁場不均一性による位相誤差マップ
phase_error_map = params.gamma * dB_map * params.TE; 

% 理想的な信号に位相誤差を乗算した、補正前の空間ドメインの信号
weighted_signal = rho .* exp(-1i * phase_error_map);

% k空間信号を計算 (空間ドメインでのFT)
k_space_data = fftshift(fft2(ifftshift(weighted_signal)));

% --- 4. 画像再構成 ---
reconstructed_image = ifftshift(ifft2(fftshift(k_space_data)));
reconstructed_image = abs(reconstructed_image);

% --- 5. 結果の表示 ---
figure;
subplot(1, 3, 1);
imagesc(rho); axis image; colormap gray; title('入力ファントム (T1WI)');

subplot(1, 3, 2);
imagesc(dB_map * 1e6); axis image; colormap default; colorbar;
title('背景磁場不均一性マップ (\mu T)');

subplot(1, 3, 3);
imagesc(reconstructed_image); axis image; colormap gray;
title('再構成画像 (dB不均一性による歪み - 最終修正)');