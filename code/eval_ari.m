function [aris, matches] = eval_ari(js, labels)
    aris = 0;
    matches = cell(2, length(js.Ls));
    for i = 1:length(js.Ls)
        L = js.Ls(i);
        match = zeros(length(js.ids), L);
        for j= 1:length(js.ids)
            idx = js.ids(j);
            match(j, :) = labels(idx, js.subseqs(j, i):js.subseqs(j, i)+L-1);
        end
        matches{1, i} = match;
        base = mode(match);
        match(match~=base) = -1;
        match(match==base) = 1;
        match(match==-1) = 0;
        matches{2, i} = match;
        aris = aris + sum(match, 'all')/length(match(:));
    end
    aris = aris/length(js.Ls);

end