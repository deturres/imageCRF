function MCRF_binarydenoising(path_name)

%% load the fake data
traindir = [ path_name '/train/'];
train_names = dir([traindir '*.jpg']);

% parameters of the problem
N     = length(train_names);  % size of training/test images
rho   = .5; % (1 = loopy belief propagation) (.5 = tree-reweighted belief propagation)
nvals = 2; % this problem is binary

fprintf('loading fake data and computing input, labels and feature maps...\n');

% in case of random generated data
% N     = 4;  % size of training images random generated
% siz   = 50; % size of training images random generated
% % load a fake dataset or randomly generate it, Basically, making noisy images, then smoothing them to make the true (discrete) output values, and then adding noise to make the input.
x = cell(1,N);
feats = cell(1,N);
efeats = cell(1,N);
for n=1:N
    % x{n} = round(imfilter(rand(siz),fspecial('gaussian',50,7),'same','symmetric')); % true label x
    
    % load your own data as true label x, add noise to create the input y
    I = double(imread(([traindir train_names(n).name])));
    img = rgb2gray(I);
    x{n}  = round(img); % true label x
    figure('Name','Loading input...','NumberTitle','off'); imshow(x{n});
    
    % extremely difficult noise pattern -- from perturbation paper
    t = rand(size(x{n})); 
    noiselevel = 1.25; % in perturbation paper 1.25
    y{n} = x{n}.*(1-t.^noiselevel) + (1-x{n}).*t.^noiselevel; % noisy input y
    figure('Name','Loading label...','NumberTitle','off'); imshow(y{n});
end




%% make features and labels.

% The features consist of simply the input image y itslef and a constant of one.
% The labels representation consists on value from of 1-nvals with 0 for unlabeled
for n=1:N
    feats{n}  = [y{n}(:) 1+0*x{n}(:)];
    labels{n} = x{n}+1;
    
%     % finding the first n max values(value<0.2)to be set to 1 (less weight to be of class 0)
%     [r,c] = find(feats{n}<0.35);
%     for m=1:size(r)        
%         feats{n}(r(m),c(m)) = 1.0;
%     end
%     figure('Name','Loading feature without extreme positive value...','NumberTitle','off'); 
%     imshow(reshape(feats{n}(:,1),size(y{1},1),size(y{1},2)));
end

%% creating the model and the edge features

% make a graph for this CRF. (A simple pairwise grid)
% siz = length(x{1});  % use this in case of squared images
% model = gridmodel(siz,siz,nvals);
[ly lx] = size(y{1});
model = gridmodel(ly,lx,nvals);

% edge features here (affect the smootheness of the result)

fprintf('computing edge features...\n')
edge_params = {{'const'},{'threshangle'},{'pairtypes'}};

% modify the computation of pca to use it as edge feature:
efeats = [];
% efeats = cell(1,N);
% parfor n=1:N
%     efeats{n} = edgeify_im(y{n},edge_params,models{n}.pairs,models{n}.pairtype);
% end

% for n=1:N    
%     data = y{n};
%     size(data)
%     [hor_efeats_ij ver_efeats_ij] = evaluate_pca(data);
%     efeats{n} = [hor_efeats_ij(:) ver_efeats_ij(:) 1+0*x{n}(:)];
% end

    % visualization function.
    % This takes a cell array of predicted beliefs as input, and shows them to the screen during training.         
    
    function viz(b_i)
        % here, b_i is a cell array of size nvals X nvars (univariate marginals)
        for n=1:N
            subplot(N,3,n    ); imshow(reshape(b_i{n}(2,:),ly,lx));
            title('predicted belief');
            subplot(N,3,n+  N); imshow(reshape(feats{n}(:,1),ly,lx));
            title('input noisy image');
            subplot(N,3,n+2*N); imshow(reshape(labels{n}-1,ly,lx));
            title('label');
            
        end
        xlabel('top: marginals  middle: input  bottom: labels')
        drawnow
    end

%% %-----------------training----------------%

fprintf('training the model (this is slow!)...\n')
% We pick a string to specify the loss and inference method. In this case, we choose truncated fitting with the clique logistic loss based on TRW with five iterations.
% Other options include 'pert_ul_trw_1e5' (perturbation, univariate logistic loss, TRW, threshold of 1e-5),
% 'em_mnf_1e5' (Surrogate Expectation-Maximization based on mean-field with a threshold of 1e-5 (simplifies to surrogate likelihood with no hidden variables),
% or 'trunc_em_trwpll_10' (Truncated surrogate EM based on multithreaded TRW with 10 iterations).
loss_spec = 'trunc_cl_trw_5'; % parameter learning and inference

% some parameters for the training optimization
crf_type  = 'linear_linear'; 
options.derivative_check = 'off';
options.viz         = @viz; % function for visualization
options.rho         = rho;
options.print_times = 1;
options.nvals       = nvals;

% figure('Name','Training...','NumberTitle','off');

% we actually optimize
% This prints a visualization while running, using the viz function above.
p = train_crf(feats,efeats,labels,model,loss_spec,crf_type,options);

% The result is a structure array p. It contains two matrices. The first, F, determines the univariate potentials. 
% Specifically, the vector of log-potentials for node i is given by multiplying F with the features for node i. 
% Similarly, G determines the log-potentials for the edge interactions. 
% If there are no edge features, though, G just multiplies a constant of 1, meaning that the 4 entries of G are themselves the log-potentials for the four possible values of (x_i,x_j).


%% %-----------------testing----------------%

% Now that we've trained the image, let's make a new test image, and get example marginals for it.

% Make a test image; use siz as dimension if you are generating dataset or using fake pregenerated
% x = round(imfilter(rand(ly,lx),fspecial('gaussian',50,7),'same','symmetric')); % label 
% t = rand(size(x));
% noiselevel = 1.25;
% y = x.*(1-t.^noiselevel) + (1-x).*t.^noiselevel; % input
% feats  = [y(:) 1+0*x(:)];
% labels = x+1;

% if using the very same image as testing image (todo: we can also load from the
% test folder)
yt=y{1};
xt=x{1};

[ly lx] = size(yt);
labelst = xt+1;
% [hor_efeats_ijt ver_efeats_ijt] = evaluate_pca(yt);
featst  = [yt(:) hor_feats_ij(:) ver_feats_ij(:) 1+0*xt(:)];

[b_i b_ij] = eval_crf(p,featst,efeats,model,loss_spec,crf_type,rho);

b_i_reshape = reshape(b_i',[ly lx nvals]);

[~,label_pred] = max(b_i_reshape,[],3);
error = mean(label_pred(:)~=labelst(:))

% (case 0(black) = yes curb!)
% computing the predicted labels considering value bigger than threshold t=0.75 as
% label 1, otherwise as label 0(yes curb!)
t = 0.9; siz = [size(b_i_reshape,1), size(b_i_reshape,2)];
[i,j] = ind2sub(siz,find(b_i_reshape(:,:,2)<t));
for n=1:size(i)        
    b_i_reshape(i(n),j(n),2) = 0.0;
end

% (case 1(white) = yes curb!) using origin image
% computing the predicted labels considering value smaller than threshold t as
% label 0, otherwise as label 1
% t = 0.80;
% siz = [size(b_i_reshape,1), size(b_i_reshape,2)];
% [i,j] = ind2sub(siz,find(b_i_reshape(:,:,1)>t));
% for n=1:size(i)        
%     b_i_reshape(i(n),j(n),2) = 1.0;
% end

% choose between max value (taking the corresponding index-->class) and
% mean value directly corresponding to the probability predicted label
[~,label_pred] = max(b_i_reshape,[],3);
% [label_pred] = mean(b_i_reshape,3);
% error = mean(label_pred(:)~=labelst(:))

% visualizing final predicted marginal and label
figure('Name','Testing..marginal ','NumberTitle','off');
subplot(2,N,1    ); imshow(reshape(b_i(2,:),ly,lx));
title('predicted marginal belief(class1))');
subplot(2,N,1+  N); imshow(reshape(b_i(1,:),ly,lx));
title('predicted marginal belief(class0))');

figure('Name','Testing..predicted belief ','NumberTitle','off');
imshow(reshape(b_i(2,:),ly,lx));

figure('Name','Testing..labels','NumberTitle','off');
subplot(N,3,1); imshow(reshape(featst(:,1),ly,lx));
title('input');
subplot(N,3,1+  N); imshow(reshape(labelst(:)-1,ly,lx));
title('true label');
hold on; axes();
subplot(N,3,1+  2*N); imshow(reshape(label_pred(:)-1,ly,lx)); % -1 to the label just in case they were computed with max()
title('predicted label');
end