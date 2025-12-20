% ファイル名の指定
filename = '"F:\hamaguchi\20251215\dual_echo\24\1_original_data\1.2.392.200036.9123.100.12.12.49941.90251215100115033358907077"'; 

% 1. メタデータ（患者情報や撮像条件）を読む
info = dicominfo(filename);
disp(info.PatientID); % IDなどを確認できます