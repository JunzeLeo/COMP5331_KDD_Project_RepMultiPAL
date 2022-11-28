function [mp, mp_I] = computeAllPairsMatrixProfiles(data, N, m, L)
    data_location = 'data_for_mp/';
    if isfolder(data_location)
        rmdir(data_location,'s');
    end
    mkdir(data_location);
    for i = 1:N
        fileID = fopen(strcat(data_location,int2str(i),'.txt'), 'w');
        for j = 1:m
            fprintf(fileID, '%f\n', data(i,j));
        end
        fclose(fileID);
    end
    total_subseqs = m-L+1;
    mp = zeros(N, N, total_subseqs);
    mp_I = zeros(N, N, total_subseqs);
    for i = 1:N
        for j = 1:N
            if i ~= j
                file1 = strcat(data_location,int2str(i),'.txt');
                file2 = strcat(data_location,int2str(j),'.txt');
                [dist,indx] = computeMatrixProfile(L, file1, file2);
                dist = 1-((dist.^2)./(2*L)); % converting distances to pearsons correlations
                if length(dist) < total_subseqs
                    for k = length(dist)+1:total_subseqs
                        dist(k) = NaN;
                    end
                end
                mp(i,j,:) = dist;
                mp_I(i,j,:) = indx;
            end
        end
    end
    save(strcat('matrix_profile_',int2str(L),'.mat'), 'mp', 'mp_I');
    rmdir(data_location,'s');
end