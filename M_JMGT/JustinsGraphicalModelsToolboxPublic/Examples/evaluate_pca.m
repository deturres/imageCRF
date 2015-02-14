function [hor_efeats_ij ver_efeats_ij] = evaluate_pca(data)

% inputs:
%   data           - a cell array of data (should be size of the sidewalkdetector_feature ly x lx)
% output:
%   hor_efeat_ij   - a cell array of edge features in the horizontal link
%   ver_efeats_ij  - a cell array of edge features in the vertical link

N = size(data);
for n = 1:N
    d = data{n};
    % build up a sphere with 8 neighbours around each cell
    [coeff, score, latent] = pca(d_neighbours);
end

if isempty(efeat)
    npairs = size(model.pairs,1);
    efeat  = ones(npairs,1);
end

x = feat;
y = zeros(model.nnodes,1);
z = efeat;
        
end[coeff, score, latent] = pca();
