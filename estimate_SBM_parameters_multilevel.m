% Function for estimating edge probabilities theta_in and theta_out,
% copying probability p, and number of communities K given a multilevel
% network A, community assignments matrix S, and hierarchical structure pi_map.
%
% Inputs:
%   - A is a cell array of length T (the number of layers) such that A{t}
%     is the Nt x Nt adjacency matrix for layer t
%   - S is a cell array of length T, where S{t} is a vector of length Nt
%     giving the community labels in layer t
%   - pi_map is a cell array of length T - 1, where pi_map{t} is a vector of length
%     N_t whose ith entry gives the node index that node i maps to in layer t + 1 
%
% Outputs:
%   - th_in and th_out, which correspond to the intra-community and
%     inter-community edge probabilities (note that the model assumes
%     degree correction); note that these parameters are global estimates,
%     for all the layers
%   - p, which is the probability that a node copies its community label
%     from the previous layer in a temporal network; also estimated
%     globally
%   - K, which is the number of communities
%
% This function calls a number of helper functions below.
%
function [th_in, th_out, p, K] = ...
    estimate_SBM_parameters_multilevel(A, S, pi_map, uniform, networktype)
  
  if nargin < 4 || isempty(uniform)
    uniform = 1;
  end
  if nargin < 5 || isempty(networktype)
    networktype = 'undirected';
  end
  
  if uniform
    [th_in, th_out] = estimate_th_in_th_out(A, S, networktype);
    K = estimate_K(S);
    p = estimate_p(S, K, pi_map);
  else
    T = length(A);  % number of layers
    th_in = zeros(T, 1);
    th_out = zeros(T, 1);
    p = zeros(T - 1, 1);
    K = zeros(T, 1);
    
    for t=1:T
      [th_in_t, th_out_t] = estimate_th_in_th_out({A{t}}, S(:, t), networktype);
      th_in(t) = th_in_t;
      th_out(t) = th_out_t;
      K(t) = estimate_K(S(:, t));
    end
    
    for t=1:T-1
      p(t) = estimate_p(S(:, t:t+1), K(t+1));
    end
  end
end

% TODO: input twom and degrees so as not to re-calculate each time
function [th_in, th_out] = estimate_th_in_th_out(A, S, networktype)
  T = length(A);
  
  if ~strcmp(networktype, 'bipartite')
    % Unipartite network
    
    % theta_in and theta_out assuming degree-corrected model
    num_in = 0; denom_in = 0;
    num_out = 0; denom_out = 0;
    for t = 1:T
      St = S{t};
      twom = nnz(A{t});
      deg = sum(A{t}, 1);
      deg_r = zeros(max(St), 1);
      twom_in = 0;

      % sum of degrees in each group, and number of within-community edges
      for r = 1:max(St)
        idx = St == r;
        deg_r(r) = sum(deg(idx));
        twom_in = twom_in + nnz(A{t}(idx, idx));
      end

      num_in = num_in + twom_in;
      num_out = num_out + twom - twom_in;
      term = sum(deg_r .^ 2) / twom;
      denom_in = denom_in + term;
      denom_out = denom_out + twom - term;
    end
    th_in = num_in / denom_in;
    th_out = num_out / denom_out;
  else
    % Bipartite network
    
    % theta_in and theta_out assuming degree-corrected model
    num_in = 0; denom_in = 0;
    num_out = 0; denom_out = 0;
    for t = 1:T
      [M, N] = size(A{t});
      St = S{t};
      mm = nnz(A{t});
      deg_rows = sum(A{t}, 2);
      deg_cols = sum(A{t}, 1);
      deg_r1 = zeros(max(St), 1);
      deg_r2 = zeros(max(St), 1);
      m_in = 0;

      % sum of degrees in each group, and number of within-community edges
      for r = 1:max(St)
        idx = St == r;
        idx1 = idx(1:M);
        idx2 = idx(M+1:M+N);
        deg_r1(r) = sum(deg_rows(idx1));
        deg_r2(r) = sum(deg_cols(idx2));
        m_in = m_in + nnz(A{t}(idx1, idx2));
      end

      num_in = num_in + m_in;
      num_out = num_out + mm - m_in;
      term = sum(deg_r1 .* deg_r2) / mm;
      denom_in = denom_in + term;
      denom_out = denom_out + mm - term;
    end
    th_in = num_in / denom_in;
    th_out = num_out / denom_out;
  end
end

% S is a cell array
function [p] = estimate_p(S, K, pi_map)
  if K == 1
    p = 1;
  else  
    pers = multilevel_persistence(S, pi_map);
    p = max((pers - 1 / K) / (1 - 1 / K), 0);
  end
end

function [K] = estimate_K(S)
  %K = length(unique(S(:)));
  %K = max(S(:));
  K = max(cellfun(@(x) max(x), S));
end
