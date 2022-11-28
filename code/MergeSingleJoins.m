function [multi_joins, n_join] = MergeSingleJoins(N, single_joins,K)
    multi_joins = struct('Ls',{},'Cs',{},'ids',{},'subseqs',{});
    n_join = length(single_joins);
    for k = K:-1:1
        all_combinations = nchoosek(1:n_join,k);
        len = size(all_combinations,1);
        mx_count = 0;
        mx_sublens = zeros(k,1);
        mx_C = zeros(k,1);
        mx_ids = zeros(N,1);
        mx_subseqs = zeros(N,k);
        mn_overlap_score = Inf;
        for i = 1:len
            join_ids = all_combinations(i,:);
            visited = zeros(N,1);
            all_subseqs = zeros(N,k);
            tmp_sublens = zeros(k,1);
            tmp_C = zeros(k,1);
            for j = 1:k
                aid = join_ids(j);
                ids = single_joins(aid).ids;
                tmp_subseqs = single_joins(aid).subseqs;
                tmp_sublens(j) = single_joins(aid).sublen;
                tmp_C(j) = single_joins(aid).C;
                n_ids = size(ids,1);
                for l = 1:n_ids
                    visited(ids(l)) = visited(ids(l))+1;
                    all_subseqs(ids(l),j) = tmp_subseqs(l);
                end
            end
            cnt = 0;
            tmp_ids = zeros(N,1);
            tmp_subseqs = zeros(N,k);
            for j = 1:N
                if visited(j) == k
                    cnt = cnt+1;
                    tmp_ids(cnt) = j;
                    tmp_subseqs(cnt,:) = all_subseqs(j,:);
                end
            end
            tmp_ids = tmp_ids(1:cnt);
            tmp_subseqs = tmp_subseqs(1:cnt,:);
            if cnt >= 2
                if cnt > mx_count
                    mx_count = cnt;
                    mx_sublens = tmp_sublens;
                    mx_C = tmp_C;
                    mx_ids = tmp_ids;
                    mx_subseqs = tmp_subseqs;
                    mn_overlap_score = computeOverlapScore(mx_count,k,mx_sublens,mx_subseqs);
                elseif cnt == mx_count
                    overlap_score = computeOverlapScore(cnt,k,tmp_sublens,tmp_subseqs);
                    if compareScores(overlap_score, mn_overlap_score, tmp_sublens, mx_sublens) == 1
                        mx_sublens = tmp_sublens;
                        mx_C = tmp_C;
                        mx_ids = tmp_ids;
                        mx_subseqs = tmp_subseqs;
                        mn_overlap_score = overlap_score;
                    end
                end
            end
        end
        if mx_count > 0
            multi_joins(1).Ls = mx_sublens;
            multi_joins(1).Cs = mx_C;
            multi_joins(1).ids = mx_ids;
            multi_joins(1).subseqs = mx_subseqs;
            n_join = k;
            break;
        end
    end
end

function [ok] = compareScores(overlap_score, mn_overlap_score, sublens, mx_sublens)
    ok = 0;
    if overlap_score < mn_overlap_score
        ok = 1;
    elseif overlap_score == mn_overlap_score
        dist = sum(abs(diff(sublens)));
        dist1 = sum(abs(diff(mx_sublens)));
        if dist > dist1
            ok = 1;
        elseif dist == dist1 && max(sublens) > max(mx_sublens)
            ok = 1;
        end
    end
end

function [overlap_score] = computeOverlapScore(n_id, n_subseq, sublens, subseqs)
    overlap_score = 0;
    cnt = 0;
    for i = 1:n_id
        for j = 1:n_subseq
            for k = j+1:n_subseq
                st = subseqs(i,j);
                ed = st+sublens(j)-1;
                st1 = subseqs(i,k);
                ed1 = st1+sublens(k)-1;
                stmx = max(st,st1);
                mned = min(ed,ed1);
                if stmx <= mned
                    overlap_score = overlap_score+(2*(mned-stmx+1))/(sublens(j)+sublens(k));
                end
                cnt = cnt+1;
            end
        end
    end
    if cnt > 0
        overlap_score = overlap_score/cnt;
    end
end