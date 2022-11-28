function [low_threshold] = getThreshold(matrix_profile,lowp)
    corr = matrix_profile(:);
    sz = length(corr);
    for i = 1:sz
        if isnan(corr(i))
            corr(i) = 0;
        end
    end
    low_threshold = prctile(corr,lowp);
    low_threshold = round(low_threshold,3);
    clear corr;
end