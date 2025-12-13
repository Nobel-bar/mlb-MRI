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