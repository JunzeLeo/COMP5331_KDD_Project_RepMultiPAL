function [avgc] = eval_c(data, js, c, draw)
%EVAL_C 此处显示有关此函数的摘要
%   此处显示详细说明
% 
total_c = 0;
for i = 1:length(js.Ls)
    len = js.Ls(i);
    subs_idx = js.subseqs(:, i);
    subs = zeros(length(js.ids), len);
    if draw
        close all;
    end
    for j = 1:length(js.ids)
        sub_idx = subs_idx(j);
        subs(j, :) = data(js.ids(j), sub_idx:sub_idx+len-1);
        if draw
            hold on;
            plot(subs(j, :));
        end
    end
    if draw
        pause;
    end
    corrs = 0;
    cnt = 0;
    for j = 1:size(subs, 1)-1
        for k = j+1:size(subs, 1)
            corrs = corrs + c(subs(j, :), subs(k, :));
            cnt = cnt + 1;
        end
    end
    total_c = total_c + corrs/cnt;
end
avgc = total_c/length(js.Ls);
end

