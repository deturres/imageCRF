function test_crf_linear_independent_noties

ly    = 3;
lx    = 3;
nvals = 2;
nfeat = 3;
rho   = .5;
%lossname = 'ul';
loss_spec = 'trunc_ul_trw_5';

model = gridmodel(ly,lx,nvals);

y = ceil(nvals*rand(model.nnodes,1));
x = randn(model.nnodes,nfeat);

F = randn(nvals,nfeat,model.nnodes);
G = randn(nvals,nvals,model.ncliques);

[L b_ij b_i dF dG] = crf_linear_independent_noties(model,F,G,x,y,rho,loss_spec);

e = 1e-7;
for i=1:size(F,1)
    for j=1:size(F,2)
        for k=1:size(F,3)
            F2 = F;
            F2(i,j,k) = F(i,j,k)+e;
            L2 = crf_linear_independent_noties(model,F2,G,x,y,rho,loss_spec);
            dF2(i,j,k) = (1/e)*(L2-L);
        end
    end
end

'difference of dF and numerical:'
norm(dF(:)-dF2(:))

for i=1:size(G,1)
    for j=1:size(G,2)
        for k=1:size(G,3)
            G2 = G;
            G2(i,j,k) = G(i,j,k)+e;
            L2 = crf_linear_independent_noties(model,F,G2,x,y,rho,loss_spec);
            dG2(i,j,k) = (1/e)*(L2-L);
        end
    end
end

'difference of dG and numerical:'
norm(dG(:)-dG2(:))
