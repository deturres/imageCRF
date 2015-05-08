function MCRF_multiEuropa(path_name)

%% load the data and computing labels and features map
imdir = [ path_name '/train/portion/new/']; % Valid for entire log new dataset(inside trains) or portion(inside trains/portion)
im_names = dir([imdir '*.png']); 
labdir = [ path_name '/labels/portion/new/']; % same name for entire_log_new or portion
lab_names = dir([labdir '*multi4AREA_GT.png']);

% parameters of the problem
N     = 4; %length(im_names);  % size of training images
rho   = .5; % TRW edge appearance probability
nvals = 4; % curb/sidewalkWall/buildings/background (0 in case unlabeled)
rez    = .6; % how much resolution percentage to use
cmap = [1 1 1; 1 0 0; 0 0 1; 0 1 0]; % to represent the label

fprintf('loading data and computing feature maps...\n');
% load true label x and input image from europa2_sidewalktetector
ims = cell(1,N);
labelsRGB = cell(1,N);
labels0 = cell(1,N);
labels = cell(1,N);
feats = cell(1,N);
efeats = cell(1,N);
for n=4:N
    
    % load input images
    I = double(imread(([imdir im_names(n).name])))/255;    
    img = rgb2gray(I);
    ims{n}  = img; % input images x
    figure('Name','Loading input...','NumberTitle','off'); imshow(ims{n});
    
    % compute the distance transform image
%     [D_euclidean,D_quasiEuclidean] = distance_map(I);
  
    % load labels
    L = double(imread(([labdir lab_names(n).name])))/255;
%     limg = rgb2gray(L);
    labelsRGB{n}  = L; % true label GT
    figure('Name','Loading label RGB...','NumberTitle','off'); imshow(labelsRGB{n});
    
end

%%
% compute the distance transform image
    [D_euclidean,D_quasiEuclidean] = distance_map(ims{4});

%% 
% The labels representation consists on values from  1 to nvals, with 0 for unlabeled
for n=1:N
    fprintf('new image\n');
    % reduce resolution for speed (mostly in case we use the different GRIDMAPS images)
    % DO NOT REDUCE when using small portion of the entire log if it's 0.5
    ims{n}    = imresize(ims{n},rez,'bilinear');

    % compute the label as a ly*lx matrix, whose values are classes depending on the 3-rgb channels original labels0 images
    [ly lx lz] = size(labelsRGB{n});
    l = zeros(ly,lx);
    for i=1:ly
        for j=1:lx
            l_r = labelsRGB{n}(i,j,1);
            l_g = labelsRGB{n}(i,j,2);
            l_b = labelsRGB{n}(i,j,3);
            if(l_r>0.8 & l_g<0.8 & l_b<0.8)
                l(i,j) = 2; % red means street (if not labeling area: curb side of the sidewalk)
            elseif(l_b>0.8 & l_r<0.8 & l_g<0.8)
                l(i,j) = 3; % blue means sidewalk (if not labeling area: wall side of the sidewalk)
            elseif(l_g>0.8 & l_r<0.8 & l_b<0.8)
                l(i,j) = 4; % green means building and background next to the sidewalk
            else
                l(i,j) = 1; % white means background(no points from the cloud)
            end
        end
    end
%     figure('Name','Loading label...');
%     colormap(cmap); miximshow(reshape(l,ly,lx),nvals);
    labels0{n} = l;
    % DO NOT REDUCE when using small portion of the entire log if it's 0.5
    labels{n} = imresize(labels0{n},rez,'nearest');
%     labels{n} = labels0{n};
    fprintf('label computed\n');
end
    %%
    % The features consist of simply the input image y itslef, the first two 
    % principal component for each cell in the center of  8-neighbours adjacent
    % cells, and a constant of one.
for n=1:N
    
    [hor_efeats_ij ver_efeats_ij] = evaluate_pca(ims{n});
    feats{n}  = [ims{n}(:) hor_efeats_ij(:) ver_efeats_ij(:) 1+0*labels{n}(:)]; %hor_efeats_ij(:) ver_efeats_ij(:)
    fprintf('features computed\n');    
%     % finding the first n max values(value<0.2)to be set to 1 (less weight to be of class 0)
%     [r,c] = find(feats{n}(1,1)<0.35);
%     for m=1:size(r)        
%         feats{n}(r(m),c(m)) = 1.0;
%     end
%     figure('Name','Loading feature without extreme positive value...','NumberTitle','off'); 
%     imshow(reshape(feats{n}(:,1),size(y{1},1),size(y{1},2)));
end
    fprintf('end dataset\n');

%% creating the modelyes
% the images come in slightly different sizes, so we need to make many models
% use a "hashing" strategy to not rebuild.  Start with empty giant array

% very big model_hash in case of entire_log_new
%   model_hash = repmat({[]},1000,1150);

%very small model_hash in case of small portion of entire_log_new
model_hash = repmat({[]},500,500);

fprintf('building models...\n')
for n=1:N
    [ly lx lz] = size(ims{n});
    if isempty(model_hash{ly,lx});
        model_hash{ly,lx} = gridmodel(ly,lx,nvals);
    end
end
models = cell(N,1);
for n=1:N
    [ly lx lz] = size(ims{n});
    models{n} = model_hash{ly,lx};
end

%% edge features (affect the smootheness of the result) plus splitting train and test set

fprintf('computing edge features...\n')
% instead of 'diffthresh', we will have the edge feature coming from the angles wrt the road direction
edge_params = {{'const'},{'diffthresh'},{'pairtypes'}}; 

% add here the evaliuation of the historical angle between cells
efeats = [];
% efeats = cell(1,N);
% parfor n=1:N
%     efeats{n} = edgeify_im(y{n},edge_params,models{n}.pairs,models{n}.pairtype);
% end

fprintf('splitting data into a training and a test set...\n')
% split everything into a training and test set

% with entire logs 1 and 2: one(log2) train, one(log1) test: 
% k = 2;
% [who_train who_test] = kfold_sets(N,2,k)

% with gridmaps for each log 1 and/or 2, the second parameter is how many training images you want to take into account 
k = 1;
K = 4;
[who_train who_test] = kfold_sets(N,K,k)


ims_train     = ims(who_train);
feats_train   = feats(who_train);
efeats_train  = []; % efeats(who_train);
labels_train  = labels(who_train);
labels0_train = labels0(who_train);
models_train  = models(who_train);
% imshow(ims_train{1})

ims_test     = ims(who_test);
feats_test   = feats(who_test);
efeats_test  = []; % efeats(who_test);
labels_test  = labels(who_test);
labels0_test = labels0(who_test);
models_test  = models(who_test);
% imshow(ims_test{:})

    % visualization function.
    % This takes a cell array of predicted beliefs as input, and shows them to the screen during training.         
    
%     function viz(b_i)
%         % here, b_i is a cell array of size nvals X nvars (univariate marginals)
%         for n=1:N
%             subplot(N,3,n    ); imshow(reshape(b_i{n}(2,:),ly,lx));
%             title('predicted belief');
%             subplot(N,3,n+  N); imshow(reshape(feats{n}(:,1),ly,lx));
%             title('input noisy image');
%             subplot(N,3,n+2*N); imshow(reshape(labels{n}-1,ly,lx));
%             title('label');
%             
%         end
%         xlabel('top: marginals  middle: input  bottom: labels')
%         drawnow
%     end

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
% options.viz         = @viz; % function for visualization
options.rho         = rho;
options.print_times = 1;
options.nvals       = nvals;

% figure('Name','Training...','NumberTitle','off');

% we actually optimize
% This prints a visualization while running , using the viz function above.
p = train_crf(feats_train,efeats_train,labels_train,models_train,loss_spec,crf_type,options);

% if using the entire set to train
% p = train_crf(feats,efeats,labels,models,loss_spec,crf_type,options);

% The result is a structure array p. It contains two matrices. The first, F, determines the univariate potentials. 
% Specifically, the vector of log-potentials for node i is given by multiplying F with the features for node i. 
% Similarly, G determines the log-potentials for the edge interactions. 
% If there are no edge features, though, G just multiplies a constant of 1, meaning that the 4 entries of G are themselves the log-potentials for the four possible values of (x_i,x_j).


%% %-----------------testing----------------%

fprintf('testing the model...\n')
fprintf('get the marginals for test images...\n');
for n=1:length(feats_test)
    [ly lx] = size(labels_test{n});
    [b_i b_ij] = eval_crf(p,feats_test{n},efeats_test,models_test{n},loss_spec,crf_type,rho);
    % choose between max value (taking the corresponding index-->class) and
    % mean value directly corresponding to the probability predicted label
    % [label_pred] = mean(b_i_reshape,3);
    [~,label_pred] = max(b_i,[],1);
    b_i_reshape = reshape(b_i',[ly lx nvals]);
    label_pred = reshape(label_pred,ly,lx);
    error_downsample = mean(label_pred(:)~=labels_test{n}(:))

    % Accuracy: pixelwise error
    label0 = labels0_test{n};
    % upsample predicted images to full resolution
    label_pred  = imresize(label_pred,size(label0),'nearest');
    E(n) = sum(label_pred(label0(:)>0)~=label0(label0(:)>0));
    T(n) = sum(label0(:)>0);
    error_upsample = mean(label_pred(:)~=label0(:))
    fprintf('error_upsample on test data,pred~GT: %f \n', error_upsample)
    fprintf('total pixelwise error on test data: %f \n', sum(E)/sum(T))
    
    figure('Name','Testing..input image','NumberTitle','off');
    imshow(reshape(feats_test{n}(:,1),ly,lx));
    
    M = length(feats_test);
    % visualizing final predicted marginal and label
    figure('Name','Testing..predicted marginalbelief','NumberTitle','off');
    subplot(2,2,1); imshow(reshape(b_i(1,:),ly,lx));
    title('(class 1-background)');
    subplot(2,2,2); imshow(reshape(b_i(2,:),ly,lx));
    title('(class 2-street)'); % curb side
    subplot(2,2,3); imshow(reshape(b_i(3,:),ly,lx));
    title('(class 3-sidewalk)'); % wall side
    subplot(2,2,4); imshow(reshape(b_i(4,:),ly,lx));
    title('(class 4-buildings)');
    
    figure('Name','Testing..labels true and predicted','NumberTitle','off');
    colormap(cmap); 
    label_gt = label0(:);
    [ly lx] = size(labels0_test{n});
    subplot(M,2,1   ); miximshow(reshape(label_gt,ly,lx),nvals);
    title('true label');
    subplot(M,2,1+ M); miximshow(reshape(label_pred(:),ly,lx),nvals);
    title('predicted label');
    % subplot(M,3,1+ 2*M);  imshow(reshape(label_pred(:),ly,lx)); % -1 to the label just in case they were computed with max()
end
end