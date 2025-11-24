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