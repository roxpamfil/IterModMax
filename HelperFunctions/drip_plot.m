% Visualisation of multilayer community structure given by S.
function [] = drip_plot(S)
  figure()
  imagesc(S);
  c = colorcube;
  %c = c(13:13+max(S(:)), :);
  %c = c([8 64], :);
  c = c(10:end, :);
  colormap(c);
  
  set(gca, 'TickLabelInterpreter', 'LaTeX', 'FontSize', 20);
  xlabel('Layers', 'Interpreter', 'LaTeX');
  ylabel('Nodes', 'Interpreter', 'LaTeX');
  xticks(1:size(S, 2));
end