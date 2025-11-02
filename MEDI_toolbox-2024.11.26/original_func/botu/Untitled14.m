
% [★修正: `Slice` -> `num_slices`]
for slice_idx = 1:23 % (デバッグのためスライス12のみ)
    
    fprintf('   --- スライス %d / %d を処理中 ---\n', slice_idx, num_slices);

    % 現在のスライスを3Dボリュームから抽出
    % [★修正: `Mag_4D` -> `iMag_4D`]
    original_img_slice = iMag_4D(:,:,slice_idx, echo_idx); % 強度画像を表示
    artifact_img_slice = artifact_img(:,:,slice_idx);
    direct_artifact_img_slice = direct_artifact_img(:,:,slice_idx);
    

    % (マスク適用は省略)

    % 各スライスごとに新しいFigureを作成する
    figure('Name', sprintf('Slice %d Comparison', slice_idx), 'WindowState', 'maximized');
     
    % 1行2列のグリッドの1番目（左側）
    subplot(1, 3, 1);
    imshow(original_img_slice, []);
    title(sprintf('Original (Slice %d)', slice_idx));
     
    % 1行2列のグリッドの2番目（中央）
    subplot(1, 3, 2);
    imshow(abs(artifact_img_slice), []);
    title(sprintf('Artifact (Slice %d)', slice_idx));
    
    % 1行2列のグリッドの3番目（右側）
    subplot(1, 3, 3);
    imshow(abs(direct_artifact_img_slice), []);
    title(sprintf('Direct　Artifact (Slice %d)', slice_idx));
    % Figureを強制的に今すぐ描画する
    drawnow;
     
    %% --- 8. (修正) 結果のファイル保存 (ループ内に移動) ---
     
    % 最終画像を permute (転置) する
    % (注: 以前の 2D スクリプトでは転置していた)
    % final_artifact_img_slice = permute(artifact_img_slice, [2 1]);
    final_artifact_img_slice = artifact_img_slice; % 3D処理では転置不要と仮定

end % --- スライスループ (for slice_idx) の終了 ---
