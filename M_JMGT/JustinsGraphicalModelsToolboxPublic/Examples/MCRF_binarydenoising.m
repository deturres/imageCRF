function MCRF_binarydenoising

% parameters of the problem
N     = 7; % number of training images
siz   = 50; % size of training images
rho   = .5; % TRW edge appearance probability
nvals = 2; % this problem is binary

% make a graph for this CRF. (A simple pairwise grid)
model = gridmodel(siz,siz,nvals);

% make a bunch of data. Basically, we make noisy images, then smooth them to make the true (discrete) output values, and then add noise to make the input.
for n=1:N
    
    x{n} = round(imfilter(rand(siz),fspecial('gaussian',50,7),'same','symmetric')); % true label x
    % extremely difficult noise pattern -- from perturbation paper
    t = rand(size(x{n}));
    noiselevel = 1.25; % in perturbation paper 1.25
    y{n} = x{n}.*(1-t.^noiselevel) + (1-x{n}).*t.^noiselevel; % noisy input y
    
end

% make features and labels. The features consist of simply the input image y itslef and a constant of one.
for n=1:N
    feats{n}  = [y{n}(:) 1+0*x{n}(:)];
    labels{n} = x{n}+1;
end

% no edge features here (smootheness of the result)
efeats = []; % none

    % visualization function.
    % This takes a cell array of predicted beliefs as input, and shows them to the screen during training. 
    function viz(b_i)
        % here, b_i is a cell array of size nvals X nvars
        for n=1:N
            subplot(3,N,n    ); imshow(reshape(b_i{n}(2,:),siz,siz));
            title('predicted belief');
            subplot(3,N,n+  N); imshow(reshape(feats{n}(:,1),siz,siz));
            title('input noisy image');
            subplot(3,N,n+2*N); imshow(reshape(labels{n}-1,siz,siz));
            title('label');
            
        end
        xlabel('top: marginals  middle: input  bottom: labels')
        drawnow
    end

%-----------------training----------------%

% We pick a string to specify the loss and inference method. In this case, we choose truncated fitting with the clique logistic loss based on TRW with five iterations.
% Other options include 'pert_ul_trw_1e5' (perturbation, univariate logistic loss, TRW, threshold of 1e-5),
% 'em_mnf_1e5' (Surrogate Expectation-Maximization based on mean-field with a threshold of 1e-5 (simplifies to surrogate likelihood with no hidden variables),
% or 'trunc_em_trwpll_10' (Truncated surrogate EM based on multithreaded TRW with 10 iterations).
loss_spec = 'trunc_cl_trw_5';

% some parameters for the training optimization
crf_type  = 'linear_linear';
options.derivative_check = 'off';
options.viz         = @viz; % function for visualization
options.rho         = rho;
options.print_times = 1;
options.nvals       = nvals;

% we actually optimize
% This prints a visualization while running, using the viz function above.
p = train_crf(feats,efeats,labels,model,loss_spec,crf_type,options);

% The result is a structure array p. It contains two matrices. The first, F, determines the univariate potentials. 
% Specifically, the vector of log-potentials for node i is given by multiplying F with the features for node i. Similarly, G determines the log-potentials for the edge interactions. 
% Since there are no features, though, G just multiplies a constant of 1, meaning that the 4 entries of G are themselves the log-potentials for the four possible values of (x_i,x_j).



%-----------------testing----------------%

% Now that we've trained the image, let's make a new test image, and get example marginals for it.
% make a test image
x = round(imfilter(rand(siz),fspecial('gaussian',50,7),'same','symmetric'));
t = rand(size(x));
y = x.*(1-t.^noiselevel) + (1-x).*t.^noiselevel; 
feats  = [y(:) 1+0*x(:)];
labels = x+1;

[b_i b_ij] = eval_crf(p,feats,efeats,model,loss_spec,crf_type,rho);

b_i = reshape(b_i',[siz siz nvals]);

[~,label_pred] = max(b_i,[],3);
error = mean(label_pred(:)~=labels(:))

end