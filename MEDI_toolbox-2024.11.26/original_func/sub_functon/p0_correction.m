function corrected_kspace = correct_global_phase(kspace_3d)
% グローバル位相の補正: 最大強度の点の位相を0度にする
    [max_val, max_idx] = max(abs(kspace_3d(:)));
    if max_val == 0
        corrected_kspace = kspace_3d;
        return;
    end
    p0_factor = kspace_3d(max_idx) / max_val;
    corrected_kspace = kspace_3d / p0_factor;
end