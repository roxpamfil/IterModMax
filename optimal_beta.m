% Function for calculating the optimal layer weights beta in the case when
% gamma and omega are layer-dependent.
function [beta] = optimal_beta(th_in, th_out)
  T = length(th_in);  % number of layers
  beta = zeros(T, 1);
  
  for t=1:T
    beta(t) = log(th_in(t)) - log(th_out(t)); 
  end
  beta = beta / mean(beta);
  
  % Replace NaN values with 1
  beta(isnan(beta)) = 1;
end