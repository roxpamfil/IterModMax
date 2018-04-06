%% Add paths to GenLouvain code and helper functions
addpath(genpath('../../GenLouvain2.1/'));
addpath('../');
addpath('../HelperFunctions/');

%% Load data and show drip plot of planted partition
load('temporal_example.mat');
[N, T] = size(S_true);
drip_plot(S_true);

%% Iterative algorithm
gamma0 = 1; omega0 = 1; 
[gamma, omega, ~, S] = it_mod_max_temporal(A, gamma0, omega0);
drip_plot(S);

%% Iterative algorithm with layer-dependent parameters
gamma0 = 1; omega0 = 1; beta0 = 1;
[gamma, omega, beta, S] = ...
  it_mod_max_temporal(A, gamma0, omega0, beta0, [], 0);
drip_plot(S);

%% Naive modularity maximisation
gamma0 = 1; omega0 = 1;
PP = @(S) postprocess_ordinal_multilayer(S, T);
[B, ~] = multiord(A, gamma0, omega0);
[S_naive, ~, ~] = iterated_genlouvain(B, [], [], 1, 'moverandw', [], PP);
%[S_naive, ~] = genlouvain(B, [], 0, 1, 'moverandw');
S_naive = reshape(S_naive, N, T);
drip_plot(S_naive);

