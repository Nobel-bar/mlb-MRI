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