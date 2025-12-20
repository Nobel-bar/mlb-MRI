% function Mask_ROI_CSF = extract_CSF(R2s, Mask, voxel_size, flag_erode, thresh_R2s)
% 
% if isempty(R2s)
%     Mask_ROI_CSF=[];
%     return
% end
% 
% if nargin < 5
%     thresh_R2s = 5;
% end
% if nargin < 4
%     flag_erode = 1;
% end
% 
% n_region_cen = 3;
% 
% matrix_size = size(Mask);
% 
% 
% 
% % Center region (sphere)
% [X,Y,Z] = ndgrid((1:matrix_size(1))*voxel_size(1),(1:matrix_size(2))*voxel_size(2),(1:matrix_size(3))*voxel_size(3));
% X_cen = sum(X(Mask>0))/sum(Mask(:));
% Y_cen = sum(Y(Mask>0))/sum(Mask(:));
% Z_cen = sum(Z(Mask>0))/sum(Mask(:));
% 
% radius_cen = 30;
% Mask_cen = sqrt(abs(X - X_cen).^2 + abs(Y - Y_cen).^2 + abs(Z - Z_cen).^2) <= radius_cen;
% 
% 
% if flag_erode
%     Mask = SMV(Mask, matrix_size, voxel_size, 10)>0.999;
% end
% 
% 
% 
% Mask_raw_1 = (R2s < thresh_R2s).*Mask_cen;
% CC = bwconncomp(Mask_raw_1,6);
% numPixels = cellfun(@numel,CC.PixelIdxList);
% [numPixels_sorted,idxs] = sort(numPixels,2,'descend');
% ROIs_region_cen = zeros(matrix_size);
% for i = 1:n_region_cen
%     numPixels_sorted(i);
%     idx = idxs(i);
%     ROIs_region_cen(CC.PixelIdxList{idx}) = i;
% end
% 
% Mask_raw_2 = (R2s < thresh_R2s).*Mask;
% CC = bwconncomp(Mask_raw_2,6);
% numPixels = cellfun(@numel,CC.PixelIdxList);
% [numPixels_sorted,idxs] = sort(numPixels,2,'descend');
% ROIs_region = zeros(matrix_size);
% for i = 1:length(idxs)
%     numPixels_sorted(i);
%     idx = idxs(i);
%     ROIs_region(CC.PixelIdxList{idx}) = i;
% end
% 
% % Choose regions which appear at center
% Mask_ROI_CSF = zeros(matrix_size);
% for i = transpose(unique(ROIs_region(ROIs_region_cen>0 & ROIs_region > 0)))
%     Mask_ROI_CSF(ROIs_region == i) = 1;
% end
% 
% Mask_ROI_CSF(Mask == 0) = 0;
% 
% end
% コードのご提示ありがとうございます。 エラーの原因は、コード内の n_region_cen = 3;（中心部のCSF領域を大きい順に3つ探す）という固定設定に対し、実際のデータでは条件に合う領域が 3つ未満（0個、1個、または2個）しか見つからなかったため、存在しない配列インデックス（3番目など）を参照して停止しています。
% 
% 日立のデータ（あるいは今回の被験者）の場合、脳室の形状や位置、あるいはR2*の値によって、候補が少なく検出されているようです。
% 
% 以下の修正版コードを使ってください。 42行目〜47行目の for ループの回数を、実際に見つかった数に合わせて自動調整するように書き換えました。
function Mask_ROI_CSF = extract_CSF(R2s, Mask, voxel_size, flag_erode, thresh_R2s)

if isempty(R2s)
    Mask_ROI_CSF=[];
    return
end

if nargin < 5
    thresh_R2s = 5;
end
if nargin < 4
    flag_erode = 1;
end

% ★変更点1：見つからなかった時の安全策として初期値を定義
n_region_cen = 3; 

matrix_size = size(Mask);

% Center region (sphere)
[X,Y,Z] = ndgrid((1:matrix_size(1))*voxel_size(1),(1:matrix_size(2))*voxel_size(2),(1:matrix_size(3))*voxel_size(3));
X_cen = sum(X(Mask>0))/sum(Mask(:));
Y_cen = sum(Y(Mask>0))/sum(Mask(:));
Z_cen = sum(Z(Mask>0))/sum(Mask(:));

radius_cen = 30;
Mask_cen = sqrt(abs(X - X_cen).^2 + abs(Y - Y_cen).^2 + abs(Z - Z_cen).^2) <= radius_cen;


if flag_erode
    % ここでマスクを削りすぎている可能性があります。
    % もし修正後もCSFが取れない場合は、呼び出し側で flag_erode = 0 にしてみてください。
    Mask = SMV(Mask, matrix_size, voxel_size, 10)>0.999;
end


% --- エラーが発生していた箇所の修正 ---
Mask_raw_1 = (R2s < thresh_R2s).*Mask_cen;
CC = bwconncomp(Mask_raw_1,6);
numPixels = cellfun(@numel,CC.PixelIdxList);
[numPixels_sorted,idxs] = sort(numPixels,2,'descend');

ROIs_region_cen = zeros(matrix_size);

% ★修正箇所: 実際に見つかった数(length(idxs))と、欲しい数(n_region_cen)の小さい方に合わせる
num_found = length(idxs);
if num_found == 0
    warning('中心領域にCSF候補が見つかりませんでした。(thresh_R2sが厳しすぎるか、マスクが狭すぎます)');
else
    % ループ回数を安全な範囲に制限
    loop_count = min(n_region_cen, num_found);
    
    for i = 1:loop_count
        % numPixels_sorted(i); % ← この行は表示用あるいは無意味な参照なので削除してもOKですが、念のためコメントアウト
        idx = idxs(i);
        ROIs_region_cen(CC.PixelIdxList{idx}) = i;
    end
end
% ------------------------------------


Mask_raw_2 = (R2s < thresh_R2s).*Mask;
CC = bwconncomp(Mask_raw_2,6);
numPixels = cellfun(@numel,CC.PixelIdxList);
[numPixels_sorted,idxs] = sort(numPixels,2,'descend');
ROIs_region = zeros(matrix_size);
for i = 1:length(idxs)
    % numPixels_sorted(i); % ここも同様に削除可
    idx = idxs(i);
    ROIs_region(CC.PixelIdxList{idx}) = i;
end

% Choose regions which appear at center
Mask_ROI_CSF = zeros(matrix_size);
for i = transpose(unique(ROIs_region(ROIs_region_cen>0 & ROIs_region > 0)))
    Mask_ROI_CSF(ROIs_region == i) = 1;
end

Mask_ROI_CSF(Mask == 0) = 0;

end