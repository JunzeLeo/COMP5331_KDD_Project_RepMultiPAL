%% Run Demo
clear;
addpath('./correlation_functions/')
% A demo on how to run RepMultiPAL with given subsequence lengths
% RepMultiPAL inputs
load('../benchmark_data/birds.mat');
p = 0;      % number of prefix time series considered when pick represent
N = size(data, 1);      % number of time series
m = size(data, 2);    % length of each time series
K = 3;      % number of desired multi-way joins
Mt = 1000000;  % Memory threshold
L = [35 30 20];   % given subsequence lengths
% RepMultiPAL inputs done

delete matrix_profile_*.mat % if you want to re-use matrix profiles, delete this line
[multi_joins] = RepMultiPAL(data,N,m,K,L,Mt,p);

% displaying the joins
displayJoin(multi_joins, data, 1);

%% Evaluation Demo
if exist('labels', 'var')
    ari = eval_ari(multi_joins, labels); % compute ari 
end

% the function for calculating correlation, should be one of [IPD, DCO, kDDTW]
correlation_function = @IPD; 

% v is the final output for correlation
v = eval_c(data, multi_joins, correlation_function, false); 