%% --- 5. ヒストグラムでの比較 ---
fprintf('ヒストグラムを比較します...\n');

figure('Name', 'Histogram Comparison');
sgtitle('Field Map Histogram Comparison', 'FontWeight', 'bold');

% マスク内の値のみを抽出
iFreq_values = iFreq(Mask == 1);
RDF_values = RDF(Mask == 1);

% 1. 左側: PDF適用前 (iFreq)
subplot(1, 2, 1);
histogram(iFreq_values, 2000); % 200本のビンで表示
title('Before PDF (iFreq)');
xlabel('Field Map Value (a.u.)');
ylabel('Voxel Count');
grid on;

% 2. 右側: PDF適用後 (RDF)
subplot(1, 2, 2);
histogram(RDF_values, 2000);
title('After PDF (RDF)');
xlabel('Local Field Value (a.u.)');
ylabel('Voxel Count');
grid on;

fprintf('ヒストグラムの比較表示が完了しました。\n');