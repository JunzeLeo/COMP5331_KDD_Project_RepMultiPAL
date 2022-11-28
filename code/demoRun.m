clear;
addpath('./correlation_functions/')
% Demo 1
% A demo on how to run MultiPAL with given subsequence lengths
% MultiPAL inputs
load('../evaluation_data/birds.mat');
p = 0;      % number of prefix time series
N = size(data, 1);      % number of time series
m = size(data, 2);    % length of each time series
K = 3;      % number of desired multi-way joins
Mt = 1000000;  % Memory threshold
L = [35 30 20];   % given subsequence lengths
use_dtw = false;
% MultiPAL inputs done

% delete matrix_profile_*.mat
[multi_joins] = MultiPAL(data,N,m,K,L,Mt,p,use_dtw);

% displaying the joins
displayJoin(multi_joins, data, 1);
% [ari, m] = eval_ari(multi_joins, labels);
v = eval_c(data, multi_joins, @IPD, false);
%     pause;
%     close all;

%{
% Demo 2
% A demo on how to run MultiPAL with no prior knowledge of subsequence lengths
% MultiPAL inputs
load('sinusoidal.mat');
scamp_location = '/usr/local/SCAMP/build/SCAMP'; % use to compute matrix profiles
N = 4;  % Number of time series
m = 291;  % length of each time series
k = 3;    % number of desired multi-way joins
Mt = 1000000;  % Memory threshold
% MultiPAL inputs done

[multi_joins] = MultiPAL(data,N,m,k,[],Mt,scamp_location);

% displaying the joins
if length(multi_joins) > 0
    displayJoin(multi_joins, data, 1);
    pause;
    close all;
else
    fprintf('No Join found\n');
end
%}
% data = zeros(10, 733);
% for i=1:100:10*100
%     data(1+(i-1)/100,:) = ussimplifiedcsv((i-1)*733+1:i*733);
% end
