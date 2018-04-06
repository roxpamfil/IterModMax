function [pers] = multilevel_persistence(S, pi_map)
  T = length(S);  % number of layers
  pers = 0;
  for t=2:T
    Stmin1 = S{t - 1};
    St = S{t};
    Nt = length(St);
    pers = pers + sum(St == Stmin1(pi_map{t - 1})) / Nt;   
  end
  
  % Normalise persistence
  pers = pers / (T - 1);
end