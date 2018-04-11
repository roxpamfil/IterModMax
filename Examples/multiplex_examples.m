%% Add paths (change as needed)
addpath('../');
addpath('../HelperFunctions/')
addpath(genpath('../../GenLouvain2.1/'));

%% Load and process Lazega Law Firm data
%  (downloaded from http://deim.urv.cat/~manlio.dedomenico/data.php)
filename = 'Lazega-Law-Firm_multiplex.edges';
fileID = fopen(filename, 'r');
M = fscanf(fileID, '%d %d %d %d', [4 Inf])';
T = max(M(:, 1));  % number of layers
N = max(max(M(:, [2 3])));  % number of nodes
A = cell([1 T]);
for t = 1:T
  L = M(M(:, 1) == t, [2 3]);
  At = full(sparse(L(:, 1), L(:, 2), 1, N, N));
  A{t} = At;
end
fclose(fileID);

%% Run multiplex optimisation
gamma0 = 1.5; omega0 = 1;
K_max_alg = 20;
[gamma, omega, ~, S, Q, converged] = it_mod_max_multiplex(A, ...
  gamma0, omega0, [], K_max_alg, 1, '', '', 0, 'directed');
drip_plot(S);
