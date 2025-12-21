function [shimStruct, dicomInfo] = get_dicom_info(targetPath)
    shimStruct = struct(); dicomInfo = [];
    if ~isfolder(targetPath), error(['No folder: ' targetPath]); end
    files = dir(fullfile(targetPath, '*'));
    for i = 1:length(files)
        fname = files(i).name;
        if startsWith(fname,'.') || files(i).isdir, continue; end
        try
            info = dicominfo(fullfile(files(i).folder, fname));
            if isfield(info, 'Private_0029_1022')
                raw = info.Private_0029_1022;
                if isa(raw,'uint8')||isa(raw,'int8'), raw=char(raw'); else, raw=string(raw); end
                pts = strsplit(raw, ',');
                for k=1:length(pts)
                    it=strtrim(pts{k});
                    if contains(it,'='), kv=strsplit(it,'='); 
                        v=str2double(kv{2}); if ~isnan(v), shimStruct.(strtrim(kv{1}))=v; end
                    end
                end
                dicomInfo = info; break;
            end
        catch, continue; end
    end
end