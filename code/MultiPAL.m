%{

    MultiPAL computes k-different multi-way joins.
    Inputs:
		data = a 2-D nxm matrix where each row represents a time series.
        N = total number of time series in the dataset.
        m = length of each time series.
        k = the number of different multi-way joins the user wants.
        L = [] or a set of subsequence lengths if the user has prior knowledge
            of expected patterns lengths.
        Mt = memory threshold.
        scamp_location = if L is given and the matrix profiles for each length
            in L are already computed, set "", otherwise, give the location of
            the SCAMP executable file,which is required to compute the matrix profiles.
    Output:
    	multi_joins = a structure containing at most k-different multi-way joins,
            which has the following format:
        multi_joins = struct('Ls',{},'Cs',{},'ids',{},'subseqs',{});
        where,
            Ls = a set subsequence lengths containing at most k lengths,
                meaning |Ls| different multi-way joins are found.
            Cs = minimum correlation coefficients for each length in Ls.
            ids = the set of joined time series ids
            subseqs = a 2-D |ids|x|Ls| matrix containining the subsequences
                in the join. The i-th column of subseqs represents i-th multi-way join.
%}
function [multi_joins] = MultiPAL(data, N, m, k, L, Mt, p, use_dtw)
    tic;
    size_of_L = length(L);
    if size_of_L ~= 0
        for i = 1:size_of_L
            if ~isfile(strcat('matrix_profile_',int2str(L(i)),'.mat')) && ~use_dtw
                % Computing the matrix profiles for all pairs of time series
                [~, ~] = computeAllPairsMatrixProfiles(data, N, m, L(i));
            end
        end
        [multi_joins,n_join] = findKMultiJoins(N, m, k, L, Mt, p,use_dtw);
    else
        size_of_L = k;
        low_b = 5;
        up_b = 20;
        iteration = 1;
        multi_joins = struct('Ls',{},'Cs',{},'ids',{},'subseqs',{});
        n_join = 0;
        while iteration <= 5
            L = zeros(size_of_L,1);
            interval = (up_b-low_b)/size_of_L;
            for i = 1:size_of_L
                L(i) = ceil(((low_b+0.5+(interval*(i-1)))/100)*m);
            end
            for i = 1:size_of_L/2
                tmp = L(i);
                L(i) = L(size_of_L-i+1);
                L(size_of_L-i+1) = tmp;
            end
            for i = 1:size_of_L
                if ~isfile(strcat('matrix_profile_',int2str(L(i)),'.mat'))&& ~use_dtw
                    % Computing the matrix profiles for all pairs of time series
                    [~, ~] = computeAllPairsMatrixProfiles(scamp_location, data, N, m, L(i));
                end
            end
            [tmp_multi_joins,tmp_n_join] = findKMultiJoins(N,m,k,L,Mt,p, use_dtw);
            if tmp_n_join == k 
                multi_joins = tmp_multi_joins;
                n_join = k;
                break;
            else
                if tmp_n_join > n_join
                    multi_joins = tmp_multi_joins;
                    n_join = tmp_n_join; 
                end
            end
            iteration = iteration+1;
            low_b = low_b-1;
            up_b = up_b+1;
            size_of_L = size_of_L+2;
        end
    end
    toc;
%     fprintf('Runtime = %f ms\n', (end_time-start_time)/1000);
    save(strcat(int2str(n_join),'-joins.mat'),'multi_joins');
end

function [multi_joins,n_join] = findKMultiJoins(N,m,k,L,Mt,p,use_dtw)
    single_joins = struct('sublen',{},'C',{},'ids',{},'subseqs',{});
    n_join = 0;
    size_of_L = length(L);
     for i = 1:size_of_L
        sublen = L(i);
        if ~use_dtw
            load(strcat('matrix_profile_',int2str(sublen),'.mat'),'mp','mp_I');
        else
            load(strcat('matrix_profile_dtw_',int2str(sublen),'.mat'),'mp','mp_I');
        end
        [C] = getThreshold(mp,85);
        while 1
            single_joins = SPAL(N,m,sublen,mp,mp_I,C,Mt,p,single_joins);
            if n_join == length(single_joins)
                break;
            end
            n_join = n_join+1;
        end
        clear matrix_profile;
        clear mp_index;
    end
    [multi_joins,n_join] = MergeSingleJoins(N, single_joins, k);
end