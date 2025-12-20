
function corrected_kspace = correct_global_phase(kspace_3d)
% CORRECT_GLOBAL_PHASE 3D k空間データのGlobal Phase (p0) 補正を行う関数
%
%   [corrected_kspace, p0_factor] = correct_global_phase(kspace_3d)
%
%   入力:
%       kspace_3d : 3次元の複素数配列 (k空間データ)
%
%   出力:
%       corrected_kspace : 補正後のデータ
%       p0_factor        : 適用された補正係数 (複素数)
%
%   処理内容:
%       データの絶対値が最大となる点を探し、その点の位相成分(p0_factor)で
%       全体を除算することで、最大点の位相を0にします。

    % データの整合性チェック
    if isempty(kspace_3d)
        error('入力データが空です。');
    end

    % 1. 絶対値の最大値とそのインデックスを取得
    %    (:) を使うことで多次元配列を1次元ベクトルとして扱います
    [max_val, max_idx] = max(abs(kspace_3d(:)));

    % 2. ゼロ除算の防止
    if max_val == 0
        warning('データの最大値が0です。補正は行われません。');
        p0_factor = 1;
        corrected_kspace = kspace_3d;
        return;
    end

    % 3. 補正係数 (p0_factor) の計算
    %    最大点の複素数値を最大振幅で割ることで、位相項 (exp(i*theta)) を抽出
    val_at_max = kspace_3d(max_idx);
    p0_factor = val_at_max / max_val;

    % 4. 全体に補正を適用
    corrected_kspace = kspace_3d / p0_factor;

    % (デバッグ用出力: 必要なければコメントアウトしてください)
    fprintf('   p0補正完了: Max Val=%.2e, Phase=%.2f rad\n', ...
        max_val, angle(p0_factor));
end