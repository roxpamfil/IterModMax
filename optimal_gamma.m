% Function for calculating the optimal resolution gamma given SBM 
% probabilities th_in and th_out (degree-corrected model).
%
% Optional argument gamma_old is useful when there is only one community
% (meaning that th_out is not defined) or there is some other degeneracy.
% Then the algorithm returns a new gamma which is equal to the old one plus
% some perturbation in the range (-10%, +10%).
% 
function [gamma] = optimal_gamma(th_in, th_out, gamma_old, uniform)  
  if nargin < 3
    gamma_old = 1;
  end
      
  if nargin < 4 || uniform
    gamma = (th_in - th_out) / (log(th_in) - log(th_out));

    if isnan(th_out) || th_out == 0  || isnan(gamma)
      % Perturbate by (-10%, +10%)
      gamma = gamma_old * (1 + 0.2 * rand() - 0.1);
    end
  else
    % Calculate gamma for each layer
    T = length(th_in);
    gamma = zeros(T, 1);
    for t=1:T
      % Call this function recursively
      gamma(t) = optimal_gamma(th_in(t), th_out(t), gamma_old(t), 1);
    end
  end
end