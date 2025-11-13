filename_base = sprintf(['1f'], theta);
save_raw_data(fullfile(save_path, [filename_base, '_phase.raw']), abs(ll));

% ===================================================================
% ローカル関数定義
% ===================================================================

function save_raw_data(filepath, data)
    fid = fopen(filepath, 'w');
    if fid == -1, error('ファイルが開けませんでした: %s', filepath); end
    fwrite(fid, data, 'double');
    fclose(fid);
end
