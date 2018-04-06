% Function for calculating the optimal interlayer copling omega given
% copying probability p, number of communities K, and SBM probabilities 
% th_in and th_out (degree-corrected model).
%
% If p = 1, the theoretically-optimal omega is infinite. This code sets it
% equal to omega_max, which should be a reasonably large number (currently
% set to 1000).
%
% If omega is layer-dependent, set optional parameter 'uniform' to 0 (set
% to 1 by default).
%
function [omega] = optimal_omega(p, K, th_in, th_out, omega_old, uniform)
  omega_max = 1000;  % omega set to this value if p = 1
  
  if nargin < 5
    omega_old = ones(length(p), 1);
  end
  
  if nargin < 6 || uniform
    % Compute a single omega
    omega = log(1 + K * p / (1 - p)) / (2 * (log(th_in) - log(th_out)));
    if th_out == 0
      omega = log(1 + K * p / (1 - p)) / (2 * (log(th_in)));
    end
    omega = max(min(omega, omega_max), 0);
    
    if isnan(omega)
      % Perturbate by (-10%, +10%)
      omega = omega_old * (1 + 0.2 * rand() - 0.1);
    end
  else
    % Compute omega for each pair of adjacent layers
    T = length(th_in);  % number of layers
    omega = zeros(T - 1, 1);
    norm_factor = 2 * mean(log(th_in) - log(th_out));
    
    for t=1:T-1
      omega(t) = log(1 + K(t + 1) * p(t) / (1 - p(t))) / norm_factor;
      omega(t) = min(omega(t), omega_max);
      
      if isnan(omega)
        % Perturbate old omega by (-10%, +10%)
        omega(t) = omega_old(t) * (1 + 0.2 * rand() - 0.1);
      end
    end
  end
end
