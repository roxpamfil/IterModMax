function [Q] = multilayer_modularity(S, B, twom)
  S = S(:);
  S_mat = sparse(1:length(S), S, 1);
  Q = trace((S_mat' * B) * S_mat / twom);
end