%% Add paths (change as needed)
addpath(genpath('../../GenLouvain2.1/'));
addpath('../');
addpath('../HelperFunctions/');

%% Load bipartite multilevel network
%  Top layer has 3 type-1 nodes and 3 type-2 nodes
%  Bottom layer has 3 type-1 nodes and 6 type-2 nodes
load('multilevel_bipartite_example.mat');

%% Run iterative modularity maximisation
%  There are two fixed points, one with p = 1 (for larger omega0), and one 
%    with p = 0.67 (for smaller omega0).
gamma0 = 0.5;
omega0 = 0.1;
K_max = 10;
[gamma, omega, beta, S, Q, converged] = ...
  it_mod_max_multilevel(A, pi_map, gamma0, omega0, [], K_max, ...
  1, 'it', 'moverandw', 0, 'bipartite');
