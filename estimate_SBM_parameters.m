% Function for estimating edge probabilities theta_in and theta_out,
% copying probability p, and number of communities K given a multilayer
% network A and community assignments matrix S.
%
% Inputs:
%   - A is a cell array of length T (the number of layers) such that A{t}
%     is the adjacency matrix for layer t
%   - S is an N x T matrix of community labels, where N is the number of
%     nodes
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
function [th_in, th_out, p, K] = estimate_SBM_parameters(A, S, uniform)
  if nargin < 3
    uniform = 1;
  end
  
  if uniform
    [th_in, th_out] = estimate_th_in_th_out(A, S);
    K = estimate_K(S);
    p = estimate_p(S, K);
  else
    T = length(A);  % number of layers
    th_in = zeros(T, 1);
    th_out = zeros(T, 1);
    p = zeros(T - 1, 1);
    K = zeros(T, 1);
    
    for t=1:T
      [th_in_t, th_out_t] = estimate_th_in_th_out({A{t}}, S(:, t));
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
function [th_in, th_out] = estimate_th_in_th_out(A, S)
  T = length(A);
  
  % theta_in and theta_out assuming degree-corrected model
  num_in = 0; denom_in = 0;
  num_out = 0; denom_out = 0;
  for t = 1:T
    twom = nnz(A{t});
    deg = sum(A{t}, 1);
    deg_r = zeros(max(S(:, t)), 1);
    twom_in = 0;
    
    % sum of degrees in each group, and number of within-community edges
    for r = 1:max(S(:,t))
      idx = S(:, t) == r;
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
end

function [p] = estimate_p(S, K)
  if K == 1
    p = 1;
  else  
    pers = ordinal_persistence(S);
    p = max((K * pers - 1) / (K - 1), 0);
  end
end

function [K] = estimate_K(S)
  %K = length(unique(S(:)));
  K = max(S(:));
end
