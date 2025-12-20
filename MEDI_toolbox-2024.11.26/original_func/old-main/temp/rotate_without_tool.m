
function vol_out = perform_padded_rotation(vol_in, theta)
    [rows, cols, num_slices] = size(vol_in);
    vol_out = complex(zeros(rows, cols, num_slices));

    pad_rows_total = rows * 3;
    pad_cols_total = cols * 3;
    
    [X_in, Y_in] = meshgrid(1:pad_cols_total, 1:pad_rows_total);
    
    centerX = (pad_cols_total + 1) / 2;
    centerY = (pad_rows_total + 1) / 2;
    
    theta_rad = theta * (pi/180);
    cosT = cos(theta_rad);
    sinT = sin(theta_rad);
    
    X_shifted = X_in - centerX;
    Y_shifted = Y_in - centerY;
    
    X_orig = X_shifted * cosT + Y_shifted * sinT + centerX;
    Y_orig = -X_shifted * sinT + Y_shifted * cosT + centerY;
    
    start_r = floor((pad_rows_total - rows)/2) + 1;
    end_r   = start_r + rows - 1;
    start_c = floor((pad_cols_total - cols)/2) + 1;
    end_c   = start_c + cols - 1;

    for z = 1:num_slices
        padded_img = complex(zeros(pad_rows_total, pad_cols_total));
        padded_img(start_r:end_r, start_c:end_c) = vol_in(:,:,z);
        
        % 実部と虚部を個別に補間
        r_re = interp2(X_in, Y_in, real(padded_img), X_orig, Y_orig, 'linear', 0);
        r_im = interp2(X_in, Y_in, imag(padded_img), X_orig, Y_orig, 'linear', 0);
        
        % --- 【修正】一時変数を使ってから切り出し ---
        rotated_full = complex(r_re, r_im);
        vol_out(:,:,z) = rotated_full(start_r:end_r, start_c:end_c);
    end
end