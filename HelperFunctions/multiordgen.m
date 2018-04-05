% Modification of multiord in GenLouvain2.1 that weighs intralayer
% modularity matrices using beta parameters, and allows for different omega
% parameters between different layers
function [B, twom] = multiordgen(A, gamma, omega, beta)
  if nargin < 2
    gamma = 1;
  end

  if nargin < 3
    omega = 1;
  end

  if nargin < 4
    beta = 1;
  end
  
  N = length(A{1});
  T = length(A);

  if length(gamma) == 1
    gamma = repmat(gamma, T, 1);
  end
  
  if length(omega) == 1
    omega = repmat(omega, T - 1, 1);
  elseif size(omega, 1) == 1
    omega = omega';  
  end
  
  if length(beta) == 1
    beta = repmat(beta, T, 1);
  end

  B = spalloc(N * T, N * T, N * N * T + 2 * N * (T - 1));
  twom = 0;
  for s = 1:T
    kout = sum(A{s}, 1);
    kin = sum(A{s}, 2);
    mm = sum(kout);
    twom = twom + mm;
    indx = [1:N] + (s - 1) * N;
    B(indx,indx) = beta(s) * ((A{s} + A{s}') / 2 - ...
      gamma(s) / 2 .* ((kin * kout + kout' * kin') / mm));
  end
  omega_vec = repmat(omega', N, 1);
  B = B + spdiags(horzcat([ones(N, 1); omega_vec(:)], [omega_vec(:); ones(N, 1)]), ...
    [N -N], N * T, N * T);
  twom = twom + 2 * N * sum(omega);
end
