% --- 既存のパラメータ設定はそのまま使用 ---
% ... (省略) ...

%% 5. 🧬 k空間ハイブリッド化と歪みシミュレーション (拡張版)
% --- k空間のベース作成、歪みシミュレーション、DC補正、ハイブリッド化 ---
fprintf('\n5. k空間のハイブリッド化とReadout歪みシミュレーションを実行中...\n');

% --- 5.0 モーションとB0変化に関する追加パラメータ ---
motion_param = struct();
motion_param.B0_shift_per_TR_Hz = 50;  % 1回のTRで最大何Hz B0が変化するか (例: 50Hz)
motion_param.B0_drift_type = 'random'; % 'random' or 'sinusoidal'
motion_param.B0_drift_frequency_TR = 10; % 'sinusoidal'の場合、何TRで1サイクルか
motion_param.B0_spatial_gradient_change_Hz_per_m = 0; % B0の変化に空間勾配も付加 (Hz/m)

% --- k空間のサイズ決定 ---
Ny_total = magnification; % 拡張後のYサイズ (k空間のライン数に相当)
Nx_total = matrix_x;      % k空間の列数

% --- 空のk空間を初期化 ---
k_space_simulated_motion = complex(zeros(Nx_total, Ny_total, num_slices));

% --- 各 k_y ラインを個別にシミュレーションするループ ---
fprintf('    各 k_y ラインの収集をシミュレート中...\n');

% k_yグラディエント強度をシミュレートするためのスケール
% k_y は -Ny/2 から Ny/2-1 までのインデックスを持つと考える
k_y_indices_sim = ((0:Ny_total-1) - floor(Ny_total/2)); 

% FOV_y を計算
FOV_y = params.FOV * (magnification / matrix_y); % 拡張後のFOV

% k空間の最大値 (k_y,max = 1 / (2*dy)) を想定
% 簡略化のため、k_y の物理単位は省略し、相対的な強度として扱う
% 厳密には k_y = (gamma / (2*pi)) * G_PE_amp * t_PE

% 各スライス、各 k_y ラインに対してシミュレーション
for slice_idx = 1:num_slices
    fprintf('        スライス %d/%d を処理中...\n', slice_idx, num_slices);
    
    current_original_slice = extend_org_rotated(:,:,slice_idx); % 回転済みのスライス画像
    current_highpass_slice = extend_high_rotated(:,:,slice_idx); % 背景磁場影響下のスライス (位相誤差源)

    for ky_line_idx = 1:Ny_total % 各 k_y ライン（PEステップ）ごとに
        
        % --- 5.1. 体動によるB0変化の決定 (このTRサイクル中) ---
        % k_yライン収集時間スケールでのB0変化をモデル化
        current_b0_shift_Hz = 0; % このTRでの中心B0オフセット

        if strcmp(motion_param.B0_drift_type, 'random')
            current_b0_shift_Hz = motion_param.B0_shift_per_TR_Hz * (2*rand - 1); % -B0_shift から +B0_shift
        elseif strcmp(motion_param.B0_drift_type, 'sinusoidal')
            current_b0_shift_Hz = motion_param.B0_shift_per_TR_Hz * sin(2*pi * (ky_line_idx / motion_param.B0_drift_frequency_TR));
        end
        
        % --- 5.2. 位相エンコード (PE) のシミュレーション ---
        % 現在のk_yラインに対応するPEグラディエントの「意図された」強度
        % 実際にはG_PEが印加されてからリードアウトが始まるまでに位相が蓄積される
        intended_pe_phase_factor = exp(1i * k_y_indices_sim(ky_line_idx) * (2*pi / Ny_total) * matrix_y); % PE軸上の位置エンコード

        % 背景磁場変化による位相誤差の蓄積 (PEグラディエント印加中と仮定)
        % これは、PEパルスが印加される短い時間 t_PE の間に B0 が変動したと近似
        t_PE = params.dwell_time * 10; % 例: PEパルスがFEの10倍の時間かかると仮定
        b0_induced_pe_error = exp(1i * params.gamma * current_b0_shift_Hz * t_PE);

        % 全体のPE位相に B0_induced_pe_error が乗る (PE軸方向のゴーストやぼけ)
        % Image_after_PE = current_highpass_slice .* intended_pe_phase_factor .* b0_induced_pe_error; % 位相が位置とB0で複雑に絡む

        % もしくは簡略化: PEグラディエントの「設定値」に、B0変化によるオフセットが加わり、結果的に異なるk_yを収集してしまう
        % ここでは、PEパルス後のデータ全体に B0 変化による位相エラーが乗ると仮定
        
        % --- 5.3. 周波数エンコード (FE) のシミュレーション ---
        % Readout期間中のB0変化を、周波数オフセットとしてモデル化
        
        % B0変化による周波数オフセット (Hz)
        delta_f_b0 = current_b0_shift_Hz; 
        
        % B0変化が空間勾配も持つ場合 (オプション)
        % 例: ΔB/Δx の空間勾配が追加された場合
        b0_spatial_gradient_hz_per_pixel = motion_param.B0_spatial_gradient_change_Hz_per_m * (params.FOV / matrix_x);
        x_coords_m = ((0:Nx_total-1) - floor(Nx_total/2)) * (params.FOV / matrix_x);
        
        b0_spatial_offset = x_coords_m * b0_spatial_gradient_hz_per_pixel; % x軸に沿った追加の周波数オフセット
        
        % k空間ラインごとのデータ生成
        % 空間画像に直接B0変化の影響を加える
        
        % B0変化による位相オフセット (時間依存)
        % Readout期間中にB0が変化すると、各点(x,y)で異なる位相誤差が蓄積
        % -> これはFE方向の歪みとPE方向のゴーストの両方に寄与
        
        % シミュレーション簡略化のため、まず空間ドメインでB0変化による周波数シフトを画像に適用
        % （ここでは、B0変化が一定の周波数オフセットをもたらすと近似）
        
        % FE軸歪みのシミュレーション:
        % B0変化による周波数オフセットがリードアウト中に加わった場合、
        % その周波数オフセットは画像のFE方向のシフトとして現れる。
        % これは空間ドメインでの画像を FE軸方向にシフトさせるFFTのプロパティを使う
        
        % FE軸シフトは空間ドメインでの位相勾配に相当するが、
        % ここでは周波数エンコード後のk空間に直接影響を与える形で表現する
        
        % まず、B0変化による位相誤差を時間領域で考える
        % -> ここは非常に複雑になるため、より簡略化した方法を取る
        
        % 【簡略化されたアプローチ】
        % 1. PE軸への影響: k空間ライン全体にかかるB0変化による位相オフセット (ゴースト)
        % 2. FE軸への影響: リードアウト中のB0変化による周波数オフセット (歪み)
        
        % --- 空間画像にB0変化の影響を適用 ---
        % 背景磁場変化による位相誤差は、空間画像にそのまま乗算される
        % Δφ_B0 = γ * ΔB * TE_eff
        % ただし、ΔBはTRごとに変動するため、ここでは各k_yラインで異なるものとする
        
        % B0_induced_phase_map: 空間的に変化するB0オフセットが、リードアウト時間(TE_eff)中に蓄積する位相
        % Readout期間中にB0が変化すると、各x方向のデータ収集時に周波数オフセットと位相オフセットが生じる
        % これを、ここでは「空間画像に適用される位相誤差マップ」としてモデル化
        
        % まず、ベースとなる背景磁場位相マップ (Highpass_img にはRDFがあるが、ここではBack_imgを使う)
        Background_Phase_Map_Rotated_slice = angle(extend_back_rotated(:,:,slice_idx));
        
        % この背景磁場マップに、さらにTRごとのB0シフトを加える
        % このB0シフトが「TR中にどこで」「どのくらいの時間」かかっているかをモデル化する必要がある
        % ここでは、単純にこのTR中の B0_shift_Hz が背景磁場の位相に加わると仮定
        
        % FE軸上の位置 x におけるB0変化による周波数オフセット
        % このオフセットが、FEエンコード後の信号に影響を与える
        
        % 実際のB0変化による空間的な位相誤差マップ (FE軸歪みの原因)
        % これは、空間的なB0不均一性と、このTR中のグローバルなB0オフセットを組み合わせる
        
        % リードアウト期間中のB0変化による位相誤差をモデル化
        % Readout期間 (T_read) 中に B0_shift_Hz が継続したと近似
        T_read = Nx_total * params.dwell_time; % リードアウト時間
        
        % B0変化による位相誤差マップ (空間に依存)
        % TR_B0_change_map(x,y) = current_b0_shift_Hz + B0_spatial_gradient_change_Hz_per_m * x
        % この周波数オフセットが T_read の間に位相として蓄積される
        
        % 各ピクセルにおける周波数オフセット (Hz)
        pixel_b0_freq_offset = current_b0_shift_Hz + b0_spatial_offset; % (1xNx_total)
        
        % 各ピクセルに蓄積される位相誤差
        % (y方向は均一と仮定して repmat)
        phase_error_map_per_pixel = repmat(exp(1i * 2*pi * pixel_b0_freq_offset * T_read), [Ny_total, 1])'; % Nx x Ny
        
        % --- 空間画像にモーションの影響を乗せる ---
        % Highpass_img (背景磁場以外の情報) + B0変化による歪み
        Image_with_Motion_Distortion = current_highpass_slice .* phase_error_map_per_pixel;
        
        % --- k空間へ変換 (このk_yラインのみを抽出) ---
        k_space_temp = fft2(fftshift(Image_with_Motion_Distortion));
        
        % k_space_simulated_motion(:, ky_line_idx, slice_idx) = k_space_temp(:, ky_line_idx); % <- これは違う
        
        % 目的: この k_y_line_idx に対応する k_space_temp の「行」を
        % k_space_simulated_motion の ky_line_idx 行にコピーする
        
        % k_space_temp は画像全体からFFTされたk空間。
        % この中から、現在の「意図された」k_yラインインデックスに対応する行を抽出する
        % 通常、fft2(image) の結果は kx_idx, ky_idx の順
        
        % k_space_simulated_motion の y軸インデックスは ky_line_idx
        % k_space_temp は Nx_total x Ny_total
        
        % ここが最も重要なポイント：
        % B0変化によって、実際に収集される k_y ラインは「意図された」k_y_line_idx とは異なるかもしれない。
        % しかし、装置は「意図された」k_y_line_idx にデータを入れる。
        % そのため、Image_with_Motion_Distortion を FFT して得られた k_space_temp の
        % 特定の k_y ラインを、シミュレーション結果の k_space_simulated_motion の
        % 現在の ky_line_idx に代入するのが正しい。
        
        % FFT結果のk空間のy軸（第2次元）の中から、
        % 現在処理しているPEステップに対応するデータ（またはそれに近いもの）を抽出する
        % 最も単純には、Image_with_Motion_Distortion全体のFFT結果を現在のk_yラインにコピーする
        % (これは非常に簡略化されたアプローチ)
        
        % 別の考え方: B0変化はPE軸とFE軸両方の位相誤差を生み出す
        % FE軸の位相誤差 -> 歪み
        % PE軸の位相誤差 -> ゴースト
        
        % ここでは、Readout期間中にB0が変化した場合のk空間データ収集をモデル化する
        % 1. 画像にB0変化による空間位相誤差を適用
        % 2. 空間画像をFFTし、そのk空間データを収集したk_yラインに代入

        % 簡略化されたk空間ライン生成（B0変化が空間ドメインの位相に影響を与えるモデル）
        k_space_line_from_image = fft(Image_with_Motion_Distortion(ky_line_idx, :)); % 現在のy方向の空間ラインをFFT
        k_space_simulated_motion(:, ky_line_idx, slice_idx) = k_space_line_from_image;
        
    end % end of ky_line_idx loop
end % end of slice_idx loop


% --- 5.4 ハイブリッド化 (今回はB0変化シミュレーションに集中するため省略 or 調整) ---
% ... (既存のハイブリッド化ロジックは、必要に応じてこのk_space_simulated_motionに適用) ...
% 例: base_and_simulate_space = k_space_simulated_motion;

base_and_simulate_space = k_space_simulated_motion; % モーションシミュレーション結果を使用
base_and_direct_space = fft2_3d_slice_by_slice(extend_org_rotated); % 直接回転後のk空間

% DC成分補正 (必要であればここで適用)
% ... (省略) ...

clear k_space_simulated_motion; % メモリ解放

fprintf('    物理シミュレーション完了。\n');

% ... (以降のステップ6と7はそのまま) ...

% ===================================================================
% 🛠️ Local Functions (補助関数)
% ===================================================================
% ... (既存のLocal Functionsはそのまま) ...