% Iterative modularity maximisation algorithm for multilevel networks
% which updates estimates of gamma and omega until convergence.
%
% Inputs:
%   - A: cell array of length T (the number of layers) such that A{t}
%     is the Nt x Nt intralayer adjacency matrix for layer t
%   - Pi: cell array of length T - 1, where Pi{t} is a vector of length
%     N_t (the number of nodes in layer t) whose ith entry gives the node index
%     that node i maps to in layer t + 1 
%   - gamma0: starting value for the resolution gamma
%   - omega0: starting value for the coupling omega
%   - beta0: starting value for the layer weights beta (not needed when
%     uniform = 1, see below)
%
% Optional inputs:
%   - K_max: when the number of communities detected by GenLouvain is above
%     this value, the algorithm decreases gamma by 20% (keeping omega
%     fixed); this helps in situations in which gamma is too large, leading 
%     to the detection of many small communities, which results in an even
%     larger gamma, and so on; default value K_max = 30
%     smaller communities being detected
%   - uniform: set to 1 if the parameters are uniform across layers, and to
%     0 otherwise; default value uniform = 1 
%   - algtype: set to one of 'st' (standard GenLouvain), 'it' (iterated
%     GenLouvain); 'itPP' (iterated GenLouvain with post-processing) not
%     currently implemented for multilevel networks; default value algtype = 'it'
%   - movetype: set to one of 'move', 'moverand', or 'moverandw' (see
%     GenLouvain options); default value movetype = 'moverandw'
%   - drawflag: set to 1 to show a drip plot of the multilayer community
%     structure at each step of the iteration; default value drawflag = 0
%
% Outputs:
%   - gamma, omega, beta: parameter estimates at convergence, or parameter
%     values yielding the highest-modularity partition if iteration did not
%     converge (after discarding first burn_in iterations, see below)
%   - S: multilayer community structure at convergence, or
%     highest-modularity community structure if iteration did not converge
%     (after discarding first burn_in iterations, see below)
%   - Q: modularity value at convergence, or highest modularity encountered
%     if iteration did not converge (after discarding first burn_in 
%     iterations, see below)
%   - converged: binary flag, equal to 1 if iteration converged to
%     prescribed tolerance
%
% Other notes:
%   - by default, convergence tolerance is 1e-2 for gamma and 5e-2 for
%     omega; the larger tolerance for the latter is justified by the fact 
%     that modularity maximisation is less sensitive to the choice of omega
%   - max_iter is the maximum number of iterations, set to 8 by default
%   - when the algorithm fails to converge in max_iter iterations, it does
%     two things: it discards the first few iterations (set by the
%     burn_in parameter and equal to 10 by default); and then it returns 
%     the multilayer partition S with the highest value of the normalised 
%     modularity Q from the remaining (max_iter - burn_in) iterations
% 
% Dependencies:
%   - genlouvain.m and iterated_genlouvain.m, which should be in the Matlab 
%     path
%   - estimate_SBM_parameters_multilevel.m, optimal_gamma.m, optimal_omega.m,
%     optimal_beta.m, which should be in the same directory as this file
%
function [gamma, omega, beta, S, Q, converged] = ...
  it_mod_max_multilevel(A, Pi, gamma0, omega0, beta0, K_max, ...
  uniform, algtype, movetype, drawflag, networktype)

  % Set parameter defaults
  if nargin < 5 || isempty(beta0)
    beta0 = 1;
  end
  if nargin < 6 || isempty(K_max)
    K_max = 30;
  end
  if nargin < 7 || isempty(uniform)
    uniform = 1;  
  end
  if nargin < 8 || isempty(algtype)
    algtype = 'it';
  end
  if nargin < 9 || isempty(movetype)
    movetype = 'moverandw';
  end
  if nargin < 10 || isempty(drawflag)
    drawflag = 0;
  end
  if nargin < 11 || isempty(networktype)
    networktype = 'undirected';
  end

  % Process inputs
  T = length(A);  % number of layers
  N = zeros(1, T);  % number of nodes in each layer
  if strcmp(networktype, 'bipartite')
    for t=1:T
      N(t) = sum(size(A{t}));
    end
  else
    for t=1:T
      N(t) = size(A{t}, 1);
    end
  end
  PP = @(S) postprocess_ordinal_multilayer(S, T);  % post-processing function

  gamma = gamma0; omega = omega0; beta = beta0;
  
  % Convergence settings, change if desired
  tol = 1e-2;     % convergence tolerance
  max_iter = 20;  % maximum number of iterations 
  
  % If iteration does not converge, discard first burn_in steps and take
  % the highest-modularity run from among the rest
  burn_in = 10;   
  
  % Initialise counters
  converged = 0;
  iter = 0;
  
  % Same parameters for all layers
  if uniform  
    fprintf('Initialisation: gamma = %.2f, omega = %.2f\n', ...
        gamma0, omega0);
      
    % Initialise arrays for storing partial results
    gamma_all = zeros(1, max_iter);
    omega_all = zeros(1, max_iter);
    Q_all = zeros(1, max_iter);
    S_all = cell(1, max_iter);
      
    while ~converged && iter < max_iter
      iter = iter + 1;
      
      % Run multilayer modularity with current values of omega and gamma
      
      % Construct modularity matrix B 
      if strcmp(networktype, 'undirected')
        [B, twom] = multilevel(A, Pi, gamma, omega);  
      elseif strcmp(networktype, 'directed')
        error('Not implemented!');
        %[B, twom] = multiorddir_f(A, gamma, omega);
      elseif strcmp(networktype, 'bipartite')
        [B, twom] = multilevelbipartite(A, Pi, gamma, omega);
      end
      
      % Run GenLouvain algorithm
      if strcmp(algtype, 'st')
        [S_vec, Q] = genlouvain(B, [], 0, 1, movetype);
      elseif strcmp(algtype, 'it')
        [S_vec, Q, ~] = iterated_genlouvain(B, [], [], 1, movetype);
      elseif strcmp(algtype, 'itPP')
        error('Post-processing not implemented for multilevel networks.');
        [S_vec, Q, ~] = iterated_genlouvain(B, [], [], 1, movetype, [], PP);
      end
      
      % Reshape S as a cell array and normalise Q 
      S = cell(1, T);
      for t=1:T
        S{t} = S_vec(sum(N(1:t-1))+1:sum(N(1:t)));
      end
      Q = Q / twom;
      
      % Save results from current run
      gamma_all(iter) = gamma;
      omega_all(iter) = omega;
      Q_all(iter) = Q;
      S_all{iter} = S;
      
      % Visualisation of multilayer partition
      if drawflag
        error('Visualisation not currently implemented for multilevel networks.');
      end
      
      % Use output to estimate th_in, th_out, p, K
      [th_in, th_out, p, K] = ...
        estimate_SBM_parameters_multilevel(A, S, Pi, 1, networktype);

      % Display results of current run
      fprintf('Iteration %d: gamma = %.2f, omega = %.2f, p = %.2f, K = %d\n', ...
        iter, gamma, omega, p, K);
      
      % Re-estimate gamma and omega
      if K == 1
        % Increase gamma, decrease omega 
        gamma_new = gamma * 1.2;
        omega_new = omega * 0.9;
      elseif K > K_max
        % Reduce gamma, keep omega fixed
        gamma_new = gamma * 0.8;
        omega_new = omega;
      else
        gamma_new = optimal_gamma(th_in, th_out, gamma);
        omega_new = optimal_omega(p, K, th_in, th_out, omega);
      end

      % Check for convergence (note higher tolerance for omega!)
      if max(abs(gamma_new / gamma - 1), ...
          abs(omega_new / omega - 1) / 5) < tol
        converged = 1;
      else
        gamma = gamma_new; omega = omega_new;
      end      
    end
    
    % Save highest-modularity results (ignoring first burn_in runs) if 
    % iteration did not converge
    if ~converged
      [~, idx] = max(Q_all(burn_in+1:end));
      idx = idx + burn_in;
      S = S_all{idx};
      Q = Q_all(idx);
      gamma = gamma_all(idx);
      omega = omega_all(idx);
    end
      
    % Estimated gamma and omega from final modularity-maximisation run
    fprintf('Final values: gamma = %.2f, omega = %.2f\n', ...
      gamma_new, omega_new);
    
    fprintf('\nOptimal values: gamma = %.2f, omega = %.2f\n\n', ...
      gamma, omega);
  else  % layer-dependent parameters 
    error('Layer-dependent case not currently implemented for multilevel networks');
  end
end
