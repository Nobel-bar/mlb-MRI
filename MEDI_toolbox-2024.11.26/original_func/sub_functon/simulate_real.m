% =========================================================
% ★★★ 行列演算版シミュレーション関数 (位相完全補正版) ★★★
% =========================================================
function k_space_line_signal = simulate_kx_line_fast(image_data, background_phase, kx_row_index, ky_col_indices)
    [Nx, Ny] = size(image_data);
    
    % --- 1. 画像をシフトして原点(1,1)に合わせる ---
    % これを行わないと、中心にある物体の位相計算が不安定になります。
    % 一旦、物体を左上の原点に持っていくことで、安定したDFT計算を行います。
    rho_eff = ifftshift(image_data .* background_phase); 
    
    % --- 2. 空間座標定義 (0 ～ 1) ---
    x_vec = ((0:Nx-1) / Nx).'; % [Nx x 1] (列ベクトル)
    y_vec = ((0:Ny-1) / Ny);   % [1 x Ny] (行ベクトル)
    
    % --- 3. k空間インデックスの変換 (物理周波数へ) ---
    k_x_scalar = (kx_row_index) - (floor(Nx/2) + 1);
    k_y_vec = (ky_col_indices(:)) - (floor(Ny/2) + 1);
    
    % --- 4. 信号計算 (行列演算) ---
    % まず、シフトされた画像に対して素直なDFTを計算します
    
    % X方向(縦)の積分
    exp_kx_vec = exp(-1i * 2 * pi * k_x_scalar * x_vec); 
    term_x = exp_kx_vec.' * rho_eff; 
    
    % Y方向(横)の変換
    exponent_y = -1i * 2 * pi * (y_vec.' * k_y_vec.'); 
    E_ky = exp(exponent_y);
    
    k_line_signal_smooth = term_x * E_ky;
    
    % --- 5. 【重要】位相補正 (-1)^(kx + ky) ---
    % 手順1で画像をシフトした分、k空間上の位相を元に戻す必要があります。
    % シフト量(Nx/2, Ny/2)に対応する位相回転項 (-1)^(kx + ky) を掛けます。
    % これにより、fftshift(fft2(元の画像)) と位相が完全に一致します。
    
    % (-1)^(kx) の項
    phase_shift_x = (-1) .^ k_x_scalar;
    
    % (-1)^(ky) の項 (ベクトル)
    phase_shift_y = (-1) .^ k_y_vec.';
    
    % 補正実行
    k_space_line_signal = k_line_signal_smooth .* phase_shift_x .* phase_shift_y;
    
end