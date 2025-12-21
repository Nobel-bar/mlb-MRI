function [m, f] = load_mat_data(fm, fp)
    d=load(fm); if isfield(d,'iMag'),m=d.iMag; elseif isfield(d,'Mask'),m=d.Mask; else, m=[]; end
    d=load(fp); if isfield(d,'iFreq'),f=d.iFreq; elseif isfield(d,'phase'),f=d.phase; else, f=[]; end
    m=double(m); f=double(f);
end