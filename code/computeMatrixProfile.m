%{
    Compute matrix profile between two time series.
    Inputs:
    scampLoc = SCAMP executable file with path location.
    sublen = subsequence length
    file1,file2 = input file names with  path location.
    withGPU = 1 or 0 depending on the availability of GPU
    cpuThreads = number of cpu threads used for computing if no gpu is available.
    Outputs:
    mp = matrix profile
    mp_I = matrix profile index values
%}
function [dist, indices] = computeMatrixProfile(sublen, file1, file2)
    system(join(["python scamp_abjoin.py ", string(file1), string(file2), int2str(sublen), ...
                " mp_dist.txt mp_index.txt"]));
    dist = importdata('mp_dist.txt');
    indices = importdata('mp_index.txt');
    indices = indices+1;
    delete mp_dist.txt;
    delete mp_index.txt;
end