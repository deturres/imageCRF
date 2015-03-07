function [hor_efeats_ij ver_efeats_ij] = evaluate_pca(data)

% build up a 8-neighbours matrix around each point, not considering the first/last row/coloumn,
% i.e using knn algo: kdtreeNS = KDTreeSearcher(x); [n,d]=knnsearch(kdtreeNS,x);

% inputs:
%   data           - a matrix of data (should be size of the sidewalkdetector_feature ly x lx)
% output:
%   hor_efeat_ij   - the edge feature in the horizontal link
%   ver_efeats_ij  - the edge feature in the vertical link

% Considering the first/last row/coloumn kipping the neighbours for the first/last row/coloumn, so to have always
% 8 neighbours for each cell
[r c] = size(data);
hor_efeats_ij(1:r,1) = 0.2; hor_efeats_ij(1,1:c) = 0.2;
hor_efeats_ij(1:r,c) = 0.2; hor_efeats_ij(r,1:c) = 0.2;
ver_efeats_ij(1:r,1) = 0.1; ver_efeats_ij(1,1:c) = 0.1;
ver_efeats_ij(1:r,c) = 0.1; ver_efeats_ij(r,1:c) = 0.1;

%%
for i=2:r-1
    for j=2:c-1
        % computing neighbours for cell d(i,j) (center of the matrix)        
        % d_neighbours(2,2) = d(i,j);
        fprintf('.');
        index = sub2ind(size(data), i, j);
        [Iadj , Radj, Nfound] = neighbourND(index, size(data)); 
        if Nfound < 8
            break
        end
        
        [I_adj J_adj] = ind2sub(size(data),Iadj);
        index_matrix = [I_adj' J_adj'];
        index_matrix = [index_matrix(1:4,:); i j ; index_matrix(5:end,:)];
        
        fprintf(',\n');

        d_adj = ones(1,9);
        for rr=1:size(index_matrix,1)
            d_adj(rr)=data(index_matrix(rr,1),index_matrix(rr,2));
        end
 %%       
        d_neighbours = ones(3,3);
        d_neighbours = reshape(d_adj,3,3);
        [coeff,~,latent] = pca(d_neighbours);
        if numel(latent)==0
           ver = 0.1;
           hor = 0.1;
        elseif numel(latent)<2 && latent(1)>0.1
           ver = 0.1;
           hor = latent(1);
        elseif numel(latent)<2 && latent(1)<0.1
           ver = latent(1);
           hor = 0.1;
        elseif latent(1) < latent(2)
           ver = latent(1);
           hor = latent(2);
        else
           ver = latent(2);
           hor = latent(1);
        end

        hor_efeats_ij(i,j) = hor;
        ver_efeats_ij(i,j) = ver;
    end
    
    fprintf('end row\n');

end


end
