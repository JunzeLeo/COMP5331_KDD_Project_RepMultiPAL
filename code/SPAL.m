function [single_joins] = SPAL(N,m,L,mp,mp_I,C,Mt, p,single_joins)
    n_join = size(single_joins,2);
    [ids, taken_subseqs] = get_taken_subseqs(N, L, single_joins);
    RN = size(ids,1);
    if RN < 2
        return;
    end
    total_subseq = m-L+1;
    matching = zeros(RN,total_subseq);
    match_limit = 0.1;
    for i = 1:RN
        for j = 1:total_subseq
            for k = 1:n_join
                if taken_subseqs(ids(i),2*k-1) == 0
                    continue;
                end
                mxst = max(j,taken_subseqs(ids(i),2*k-1));
                mned = min(j+L-1,taken_subseqs(ids(i),2*k));
                count = mned-mxst+1;
                if count > 0
                    matching(i,j) = max(matching(i,j), count/L);
                end
            end
        end
    end
    candidates = cell(RN*total_subseq,6);
    candidate_size = 0;
    for i = 1:RN
        for j = 1:total_subseq
            if matching(i,j) < match_limit
                mx_corr = C;
                mx_idx = -1;
                for k = 1:RN %% original 1:RN
                    idx = mp_I(ids(i),ids(k),j);
                    corr = mp(ids(i),ids(k),j);
                    if i ~= k && idx > 0 && corr >= mx_corr && ...
                            matching(k, idx) < match_limit
                        mx_corr = corr;
                        mx_idx = k;
                    end
                end
                if mx_idx > 0
                    candidate_size = candidate_size+1;
                    if p == 0
                        candidates(candidate_size,:) = {[i;mx_idx],j, mp_I(ids(i),ids(mx_idx),j),mx_corr,mx_idx,[1, mp(ids(i),ids(mx_idx),j);mp(ids(i),ids(mx_idx),j), 1]};
                    else
                        candidates(candidate_size,:) = {[i;mx_idx],[j, mp_I(ids(i),ids(mx_idx),j)], mp_I(ids(i),ids(mx_idx),j),mx_corr,mx_idx,[1, mp(ids(i),ids(mx_idx),j);mp(ids(i),ids(mx_idx),j), 1]};
                    end
                end
            end
        end
    end
    if candidate_size == 0
        clear candidates;
        return;
    end
    candidates = sortrows(candidates(1:candidate_size,:), 4, 'descend');
    candidate_size = min(candidate_size, Mt);
    candidates = candidates(1:candidate_size,:);
    if RN > 2
        iteration = 2;
        while iteration <= RN
            is_bad = zeros(candidate_size,1);
            all_bad = 1;
            for i = 1:candidate_size
                cand = candidates(i,:);
                sz = size(cand{1},1);
                newset = zeros(sz+1,1);
                newset(1:sz) = cand{1};
                is_used = zeros(RN,1);
                for j = 1:sz
                    is_used(newset(j)) = 1;
                end
                k = cand{5};
                corr_mat = cand{6};
                idx = cand{3};
                mx_id = -1;
                mx_corr = C;
                mx_subid = -1;
                for j = 1:RN
                    if is_used(j) == 1
                        continue;
                    end
                    idx1 = mp_I(ids(k),ids(j),idx);
                    corr = mp(ids(k),ids(j),idx);
                    if idx1 > 0 && corr >= mx_corr && matching(j,idx1) < match_limit
                        mx_corr = corr;
                        mx_id = j;
                        mx_subid = idx1;
                    end
                end
                if mx_id > 0
                    newset(sz+1) = mx_id;
                    if sz+1<=p
                        newcorr_mat = ones(sz+1, sz+1);
                        newcorr_mat(1:sz, 1:sz) = corr_mat;
                        for s = 1:sz
                            newcorr_mat(s, sz+1) = mp(cand{1}(s), mx_id, cand{2}(s));
                            newcorr_mat(sz+1, s) = mp(mx_id, cand{1}(s), mx_subid);
                        end
                        [~, new_repid] = max(sum(newcorr_mat));
                        newk = newset(new_repid);
                        cand{2} = [cand{2}, mx_subid];
                        newk_idx = cand{2}(new_repid);
                        candidates(i,:)= {newset,cand{2},newk_idx,min(cand{4},mx_corr), newk, newcorr_mat};
                    elseif p~= 0
                        candidates(i,:)= {newset,[cand{2}, mx_subid],cand{3},min(cand{4},mx_corr), cand{5}, cand{6}};
                    else 
                        candidates(i,:)= {newset,[cand{2}, mx_subid],mx_subid,min(cand{4},mx_corr), mx_id, cand{6}};
                    end
                    all_bad = 0;
                else
                    is_bad(i) = 1;
                end
            end
            if all_bad == 0
                cur = 0;
                for i = 1:candidate_size
                    if is_bad(i) == 0
                        cur = cur+1;
                        candidates(cur,:) = candidates(i,:);
                    end
                end
                candidate_size = cur;
                candidates = candidates(1:candidate_size,:);
                iteration = iteration+1;
                if RN == iteration
                    break;
                end
            else
                break;
            end
        end
    end
    mx_idx = -1;
    mx_corr = 0;
    for i = 1:candidate_size
        cand = candidates(i,:);
        if cand{4} > mx_corr
            mx_idx = i;
            mx_corr = cand{4};
        end
    end
    if mx_idx == -1
        clear candidates;
        return;
    end
    n_join = n_join+1;
    single_joins(n_join).sublen = L;
    single_joins(n_join).C = C;
    cand = candidates(mx_idx,:);
    sz = size(cand{1},1);
    tmp_ids = cand{1};
    idx = cand{2}(1);
    act_ids = zeros(sz,1);
    subseqs = zeros(sz,1);
    act_ids(1) = ids(tmp_ids(1));
    subseqs(1) = idx;
    for i = 2:sz
        act_ids(i) = ids(tmp_ids(i));
        idx = mp_I(act_ids(i-1),act_ids(i),idx);
        subseqs(i) = idx;
    end
    single_joins(n_join).ids = act_ids;
    single_joins(n_join).subseqs = subseqs;
    clear candidates;
end

function [ids,taken_subseqs] = get_taken_subseqs(N, cur_sublen, single_joins)
    n_join = size(single_joins,2);
    taken_subseqs = zeros(N,n_join*2);
    used = zeros(N,1);
    for i = 1:n_join
        sublen = single_joins(i).sublen;
        ids = single_joins(i).ids;
        subseqs = single_joins(i).subseqs;
        n_ids = size(ids,1);
        for j = 1:n_ids
            taken_subseqs(ids(j),2*i-1) = subseqs(j);
            taken_subseqs(ids(j),2*i) = subseqs(j)+sublen-1;
            if sublen == cur_sublen
                used(ids(j)) = 1;
            end
        end
    end
    RN = N-sum(used);
    ids = zeros(RN,1);
    j = 0;
    for i = 1:N
        if used(i) == 0
            j = j+1;
            ids(j) = i;
        end
    end
end
