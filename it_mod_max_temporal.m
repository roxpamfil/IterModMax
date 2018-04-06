% Iterative modularity maximisation algorithm for temporal networks
% which updates estimates of gamma and omega until convergence.
%
% Inputs:
%   - A: cell array of length T (the number of layers) such that A{t}
%     is the adjacency matrix for layer t
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
%     GenLouvain), or 'itPP' (iterated GenLouvain with post-processing);
%     default value algtype = 'itPP'
%   - movetype: set to one of 'move', 'moverand', or 'moverandw' (see
%     GenLouvain options); default value movetype = 'moverandw'
%   - drawflag: set to 1 to show a drip plot of the multilayer community
%     structure at each step of the iteration; default value drawflag = 0
%   - networktype: one of 'undirected' (default), 'directed', or 'bipartite'
%     (which implicitly assumes undirected edges)
%
% Outputs:
%   - gamma, omega, beta: parameter estimates at convergence, or parameter
%     values yielding the highest-modularity partition if iteration did not
%     converge
%   - S: multilayer community structure at convergence, or
%     highest-modularity community structure if iteration did not converge
%   - Q: modularity value at convergence, or highest modularity encountered
%     if iteration did not converge
%   - converged: binary flag, equal to 1 if iteration converged to
%     prescribed tolerance
%
% Other notes:
%   - by default, convergence tolerance is 1e-2 for gamma and 5e-2 for
%     omega; the larger tolerance for the latter is justified by the fact 
%     that modularity maximisation is less sensitive to the choice of omega
%   - max_iter is the maximum number of iterations, set to 20 by default
%   - when the algorithm fails to converge in max_iter iterations, it 
%     returns the multilayer partition S with the highest value of the 
%     normalised modularity Q
% 
% Dependencies:
%   - genlouvain.m and iterated_genlouvain.m, which need to be in a 
%     directory "../GenLouvain2.1" or otherwise in the Matlab path
%   - estimate_SBM_parameters.m, optimal_gamma.m, optimal_omega.m,
%     optimal_beta.m, which need to be in the same directory as this file
%
function [gamma, omega, beta, S, Q, converged] = ...
  it_mod_max_temporal(A, gamma0, omega0, beta0, K_max, ...
  uniform, algtype, movetype, drawflag, networktype)

  % Add path to GenLouvain code
  addpath(genpath('../GenLouvain2.1/'));
  addpath('./HelperFunctions/');

  % Set parameter defaults
  if nargin < 4 || isempty(beta0)
    beta0 = 1;
  end
  if nargin < 5 || isempty(K_max)
    K_max = 30;
  end
  if nargin < 6 || isempty(uniform)
    uniform = 1;  
  end
  if nargin < 7 || isempty(algtype)
    algtype = 'itPP';
  end
  if nargin < 8 || isempty(movetype)
    movetype = 'moverandw';
  end
  if nargin < 9 || isempty(drawflag)
    drawflag = 0;
  end
  if nargin < 10 || isempty(networktype)
    networktype = 'undirected';
  end

  % Process inputs
  T = length(A);  % number of layers
  N = length(A{1});  % number of nodes in each layer
  if strcmp(networktype, 'bipartite')
    N = sum(size(A{1}));
  end
  PP = @(S) postprocess_ordinal_multilayer(S, T);  % post-processing function

  gamma = gamma0; omega = omega0; beta = beta0;
  
  % Convergence settings, change if desired
  tol = 1e-2;     % convergence tolerance
  max_iter = 20;  % maximum number of iterations 
  converged = 0;
  iter = 0;
  
  fprintf('Initialisation: gamma = %.2f, omega = %.2f\n', ...
        gamma0, omega0);
  
  % Same parameters for all layers
  if uniform  
    Q_best = 0;
    while ~converged && iter < max_iter
      iter = iter + 1;
      
      % Run multilayer modularity with current values of omega and gamma
      
      % Construct modularity matrix B
      if strcmp(networktype, 'undirected')
        [B, twom] = multiord(A, gamma, omega);  
      elseif strcmp(networktype, 'directed')
        [B, twom] = multiorddir_f(A, gamma, omega);
      elseif strcmp(networktype, 'bipartite')
        [B, twom] = multiordbipartite(A, gamma, omega);
      end
      
      % Run GenLouvain algorithm
      if strcmp(algtype, 'st')
        [S, Q] = genlouvain(B, [], 0, 1, movetype);
      elseif strcmp(algtype, 'it')
        [S, Q, ~] = iterated_genlouvain(B, [], [], 1, movetype);
      elseif strcmp(algtype, 'itPP')
        [S, Q, ~] = iterated_genlouvain(B, [], [], 1, movetype, [], PP);
      end
      
      % Reshape S and normalise Q
      S = reshape(S, N, T);
      Q = Q / twom;
      
      % Save run if best so far
      if Q > Q_best
        Q_best = Q; 
        S_best = S;
        gamma_best = gamma;
        omega_best = omega;
      end
      
      % Visualisation of multilayer partition
      if drawflag
        drip_plot(S); 
        title(strcat('$\gamma$=', num2str(gamma, '%.2f'), ...
          ', $\omega=$', num2str(omega, '%.2f'), ...
          ', $Q=$', num2str(Q, '%.2f')), ...
          'Interpreter', 'LaTeX');
        drawnow
      end
      
      % Use output to estimate th_in, th_out, p, K
      [th_in, th_out, p, K] = estimate_SBM_parameters(A, S);

      % Re-estimate gamma and omega
      if K == 1
        % Increase gamma, decrease omega 
        gamma_new = gamma * 1.2;
        omega_new = omega * 0.9;
      elseif K > K_max
        % Reduce gamma, keep omega fixed
        gamma_new = gamma * 0.8;
        omega_new = omega;
        %omega_new = omega * 1.1;
      else
        gamma_new = optimal_gamma(th_in, th_out, gamma);
        omega_new = optimal_omega(p, K, th_in, th_out, omega);
      end

      % Display outputs
      fprintf('Iteration %d: gamma = %.2f, omega = %.2f, p = %.2f, K = %d\n', ...
        iter, gamma_new, omega_new, p, K);

      % Check for convergence (note higher tolerance for omega!)
      if max(abs(gamma_new / gamma - 1), ...
          abs(omega_new / omega - 1) / 5) < tol
        converged = 1;
      end
      
      gamma = gamma_new; omega = omega_new;
      
    end
    
    % Save highest-modularity results if iteration did not converge
    if ~converged
      S = S_best;
      Q = Q_best;
      gamma = gamma_best;
      omega = omega_best;
    end
      
    fprintf('\n')
  else  % layer-dependent parameters 
    % Create vectors of appropriate length if scalars are given
    if length(gamma) == 1
      gamma = repmat(gamma, T, 1);
    end
    if length(omega) == 1
      omega = repmat(omega, T - 1, 1);
    end
    if length(beta) == 1
      beta = repmat(beta, T, 1);
    end
    
    Q_best = 0;
    
    while ~converged && iter < max_iter
      iter = iter + 1;

      % Run multilayer modularity with current values of omega and gamma
      [B, twom] = multiordgen(A, gamma, omega, beta); 
      if strcmp(algtype, 'st')
        [S, Q] = genlouvain(B, [], 0, 1, movetype);
      elseif strcmp(algtype, 'it')
        [S, Q, ~] = iterated_genlouvain(B, [], [], 1, movetype);
      elseif strcmp(algtype, 'itPP')
        [S, Q, ~] = iterated_genlouvain(B, [], [], 1, movetype, [], PP);
      end
      
      % Reshape S and normalise Q
      S = reshape(S, N, T);
      Q = Q / twom;
      
      % Save run if best so far
      if Q > Q_best
        Q_best = Q; 
        S_best = S;
        gamma_best = gamma;
        omega_best = omega;
        beta_best = beta;
      end
      
      % Show visualisation of multilayer partition
      if drawflag
        drip_plot(S); 
        title(strcat('$\langle\gamma\rangle$=', num2str(mean(gamma), '%.2f'), ...
          ', $\langle\omega\rangle=$', num2str(mean(omega), '%.2f'), ...
          ', $Q=$', num2str(Q, '%.2f')), ...
          'Interpreter', 'LaTeX');
        drawnow
      end
      
      % Use output to estimate th_in, th_out, p, K
      [th_in, th_out, p, K] = estimate_SBM_parameters(A, S, 0);
      
      % Re-estimate gamma, omega, beta
      gamma_new = optimal_gamma(th_in, th_out, gamma, 0);
      omega_new = optimal_omega(p, K, th_in, th_out, omega, 0);
      beta_new = optimal_beta(th_in, th_out);
      
      % Display average parameter values and maximum change
      fprintf('Iteration %d average values: gamma = %.2f, omega = %.2f, p = %.2f, K = %.1f\n', ...
        iter, mean(gamma_new), mean(omega_new), mean(p), mean(K));
      fprintf('..max change: delta(gamma) = %.2f, delta(omega) = %.2f, delta(beta) = %.2f\n', ...
        max(abs(gamma_new - gamma)), max(abs(omega_new - omega)), max(abs(beta_new - beta)));
      %fprintf('Gamma:'); gamma_new
      %fprintf('Omega:'); omega_new
      %fprintf('Beta:'); beta_new
      
      % Check for convergence
      if max(max(max(abs(gamma_new ./ gamma - 1)), ...
          max(abs(omega_new ./ omega - 1)) / 5), ...
          max(abs(beta_new ./ beta - 1)) / 1) < tol
        converged = 1;
      end

      gamma = gamma_new; omega = omega_new; beta = beta_new;
    end
    
    % Save highest-modularity results if iteration did not converge
    if ~converged
      S = S_best;
      Q = Q_best;
      gamma = gamma_best;
      omega = omega_best;
    end
      
    fprintf('\n');
  end
end
