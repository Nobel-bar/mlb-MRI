clear all; clc;

%% 1. 設定：ログファイルからの電流値入力 (手動入力)
% ここにログファイルから読み取った値を入力してください
% 順序の例: [X, Y, Z, Z2, ZX, ZY, X2-Y2, 2XY] (装置の定義に合わせてください)

% ファントム(Shimmed)撮影時の電流値 [単位: 任意(mAやA)]
% ※ No-Shim時はすべて0と仮定します
Currents_Phantom = [100, -50, 20, 10, 0, 0, 0, 0]; 

% 被験者撮影時の電流値
Currents_Subject = [120, -40, 25, 15, 5, -2, 0, 0]; 

% データの読み込み (事前に作成した変数をロードすると仮定)
% b0_phantom_noshim: ファントム(シムなし) [Hz]
% b0_phantom_shim:   ファントム(シムあり) [Hz]
% mask_phantom:      ファントムマスク
% voxel_size:        [2.0, 2.0, 2.0] など

%% 2. 座標系の作成と正規化 (球面調和関数の基底作成)
sz = size(b0_phantom_noshim);
[x, y, z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));

% 座標中心を画像の中心に合わせ、-1~1に正規化 (球半径rで割るのが理想)
center = sz / 2 + 0.5;
x = (x - center(1)) / (min(sz)/2); % 単純化のため最小辺で正規化
y = (y - center(2)) / (min(sz)/2);
z = (z - center(3)) / (min(sz)/2);

% マスク内のデータのみ抽出
ind = find(mask_phantom > 0);
xp = x(ind); yp = y(ind); zp = z(ind);

%% 3. デザイン行列 A の作成 (球面調和関数項)
% 電流値のベクトルの順番に対応させる必要があります。
% ここでは一般的な1次+2次の8項とします。
% 1: X
% 2: Y
% 3: Z
% 4: Z2  (2z^2 - x^2 - y^2) 
% 5: ZX
% 6: ZY
% 7: X2-Y2 (x^2 - y^2)
% 8: 2XY
% (※定数項(B0 offset)は別途扱います)

A = [xp, ...                       % 1. X
     yp, ...                       % 2. Y
     zp, ...                       % 3. Z
     2*zp.^2 - xp.^2 - yp.^2, ...  % 4. Z2
     zp.*xp, ...                   % 5. ZX
     zp.*yp, ...                   % 6. ZY
     xp.^2 - yp.^2, ...            % 7. X2-Y2
     2*xp.*yp];                    % 8. 2XY

% 定数項（オフセット）を追加
A_const = [ones(size(xp)), A]; 

%% 4. 係数の算出 (Fitting)
fprintf('係数を算出しています...\n');

% ファントム(No Shim) の係数算出 -> これが「磁石本来の成分」
% b = A * c  => c = A \ b
vals_noshim = b0_phantom_noshim(ind);
coeffs_noshim = A_const \ vals_noshim; 

% ファントム(Shimmed) の係数算出
vals_shim = b0_phantom_shim(ind);
coeffs_shim = A_const \ vals_shim;

% 最初の項は定数項なので、シム係数としては2番目以降を取り出す
C_noshim_vec = coeffs_noshim(2:end);
C_shim_vec   = coeffs_shim(2:end);

%% 5. シム感度 (Sensitivity) の特定
% 「電流1単位あたり、係数がどれだけ変化するか」
% Sensitivity = (Shimmed係数 - NoShim係数) ./ 電流値
% ※割り算は要素ごと。電流0の項は計算できないので注意(NaN回避)

delta_Coeffs = C_shim_vec - C_noshim_vec;
Sensitivity = zeros(size(delta_Coeffs));

% 電流が流れていたチャンネルのみ計算
valid_idx = abs(Currents_Phantom') > 1e-5; % 0でない場所
Sensitivity(valid_idx) = delta_Coeffs(valid_idx) ./ Currents_Phantom(valid_idx)';

% 結果確認
disp('算出されたシム感度 (Hz/Unit):');
disp(Sensitivity);

%% 6. 被験者用EFの再構築
fprintf('被験者の電流値に基づいてEFを再構築中...\n');

% 被験者のシム設定による係数の予測
% 予測係数 = NoShim係数 + (感度 * 被験者電流)
C_subject_vec = C_noshim_vec + Sensitivity .* Currents_Subject';

% 空間全体での再構築
% 全画素の座標ベクトルを作成
x_all = x(:); y_all = y(:); z_all = z(:);

A_all = [x_all, ...
         y_all, ...
         z_all, ...
         2*z_all.^2 - x_all.^2 - y_all.^2, ...
         z_all.*x_all, ...
         z_all.*y_all, ...
         x_all.^2 - y_all.^2, ...
         2*x_all.*y_all];

% 定数項(Magnet offset) + シム項
% coeffs_noshim(1) は定数項(全体的な周波数オフセット)
EF_vec = coeffs_noshim(1) + A_all * C_subject_vec; 

% 3次元画像に戻す
EF_estimated = reshape(EF_vec, sz);

%% 7. 最終的なSF (Susceptibility Field)
% b0_subject: 被験者の観測B0マップ
SF_subject = b0_subject - EF_estimated;

% 表示
figure;
subplot(1,3,1); imagesc(b0_subject(:,:,round(sz(3)/2))); axis image; colorbar; title('Measured Subject B0');
subplot(1,3,2); imagesc(EF_estimated(:,:,round(sz(3)/2))); axis image; colorbar; title('Estimated EF (From Logs)');
subplot(1,3,3); imagesc(SF_subject(:,:,round(sz(3)/2))); axis image; colorbar; title('Result: SF');