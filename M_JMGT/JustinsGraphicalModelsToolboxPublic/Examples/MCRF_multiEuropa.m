function MCRF_binaryEuropa(path_name)

%% load the data and computing labels and features map
imdir = [ path_name '/train/'];
im_names = dir([imdir '*.png']); % '*.png'. Just in case of entire log: *0.5.png, in case version in black use 0.5_origin_...png
labdir = [ path_name '/labels/'];
lab_names = dir([labdir '*multi3_GT.png']);

% parameters of the problem
N     = length(im_names);  % size of training images
rho   = .5; % TRW edge appearance probability
nvals = 3; % curb/sidewalkWall/background 0 is unlabeled
rez    = .6; % how much to reduce resolution
cmap = [1 1 1; 1 0 0; 0 0 1; 0 1 0]; % to represent the label

fprintf('loading data and computing feature maps...\n');
% load true label x and input image from europa2_sidewalktetector
ims = cell(1,N);
labels0 = cell(1,N);
labels = cell(1,N);
feats = cell(1,N);
efeats = cell(1,N);
for n=1:N
    
    % load input images
    I = double(imread(([imdir im_names(n).name])))/255;    
    img = rgb2gray(I);
    ims{n}  = img; % input images x
    figure('Name','Loading input...','NumberTitle','off'); imshow(ims{n});
    % load labels
    L = double(imread(([labdir lab_names(n).name])))/255;
%     limg = rgb2gray(L);
    labels0{n}  = L; % true label GT
    figure('Name','Loading label...','NumberTitle','off'); imshow(labels0{n});
    
end


%% 
% The labels representation consists on values from  1 to nvals, with 0 for unlabeled
for n=1:N
    fprintf('new image\n');
    % reduce resolution for speed, in case we use the different GRIDMAPS images
    ims{n}    = imresize(ims{n}   ,rez,'bilinear');
    % compute the label as a ly*lx matrix, whose values are classes depending on the 3-rgb channels original labels0 images
    [ly lx lz] = size(labels0{n});
    l = zeros(ly,lx);
    for i=1:ly
        for j=1:lx
            l_r = labels0{n}(i,j,1);
            l_g = labels0{n}(i,j,2);
            l_b = labels0{n}(i,j,3);
            if(l_r>0.8 & l_g<0.8 & l_b<0.8)
                l(i,j) = 1; % red means curb side of the sidewalk
            elseif(l_b>0.8 & l_r<0.8 & l_g<0.8)
                l(i,j) = 2; % blue means wall side of the sidewalk
            elseif(l_g>0.8 & l_r<0.8 & l_b<0.8)
                l(i,j) = 3; % green means building and background next to the sidewalk
            else
                l(i,j) = 0; % background unlabeled
            end
        end
    end
    figure(n);
    colormap(cmap); miximshow(reshape(l+1,ly,lx),nvals+1); % adjust the labeleing to visualize the unlabeled portion in white91 for matlab)
    labels0{n} = l;
    fprintf('label computed\n');
end
    %%
    % The features consist of simply the input image y itslef, the first two 
    % principal component for each cell in the center of  8-neighbours adjacent
    % cells, and a constant of one.
for n=1:N
    labels{n} = imresize(labels0{n},rez,'nearest');
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

model_hash = repmat({[]},1200,1200);
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
edge_params = {{'const'},{'diffthresh'},{'pairtypes'}};

% add here the evaliuation of the historical angle between cells
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

fprintf('splitting data into a training and a test set...\n')
% split everything into a training and test set

% with entire logs 1 and 2,if one(log2) train, one(log1) test: 
k = 2;
[who_train who_test] = kfold_sets(N,2,k)

% with gridmaps for each log 1 and 2, the second parameter is how many training images you want to take into account 

% k = 2;
% [who_train who_test] = kfold_sets(N,N,k)


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
% This prints a visualization while running, using the viz function above.
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

    b_i_reshape = reshape(b_i',[ly lx nvals]);

    [~,label_pred] = max(b_i_reshape,[],3);
    error_downsample = mean(label_pred(:)~=labels_test{n}(:))

    % Accuracy: pixelwise error
    label0 = labels0_test{n};
    size(label0)
    size(label_pred)
    % upsample predicted images to full resolution
    label_pred  = imresize(label_pred,size(label0),'nearest');
    E(n) = sum(label_pred(label0(:)>0)~=label0(label0(:)>0));
    T(n) = sum(label0(:)>0);
    error = mean(label_pred(:)~=label0(:))
    fprintf('error on test data,pred~GT: %f \n', error)
    fprintf('total pixelwise error on test data: %f \n', sum(E)/sum(T))

    % (case 0(black) = yes curb!)
    % computing the predicted labels considering value bigger than threshold t=0.75 as
    % label 1, otherwise as label 0(yes curb!)
    t = 0.85; siz = [size(b_i_reshape,1), size(b_i_reshape,2)];
    [i,j] = ind2sub(siz,find(b_i_reshape(:,:,2)<t));
    for h=1:size(i)        
        b_i_reshape(i(h),j(h),2) = 0.0;
    end

    % choose between max value (taking the corresponding index-->class) and
    % mean value directly corresponding to the probability predicted label
    [~,label_pred] = max(b_i_reshape,[],3);
    % [label_pred] = mean(b_i_reshape,3);
    
    M = length(feats_test);
    % visualizing final predicted marginal and label
    figure('Name','Testing..marginal','NumberTitle','off');
    subplot(M,3,1); imshow(reshape(b_i(1,:),ly,lx));
    title('predicted marginal belief(class 1-background))');
    subplot(M,3,1+   M); imshow(reshape(b_i(2,:),ly,lx));
    title('predicted marginal belief(class 2-curb))');
    subplot(M,3,1+ 2*M); imshow(reshape(b_i(3,:),ly,lx));
    title('predicted marginal belief(class 3-sidewalkWall))');

    figure('Name','Testing..labels','NumberTitle','off');
%     subplot(M,3,1); 
    imshow(reshape(feats_test{n}(:,1),ly,lx));
    title('input');
    
    figure('Name','Testing..labels true and predicted','NumberTitle','off');
    colormap(cmap); 
    subplot(M,2,1   ); miximshow(reshape(labels_test{n}(:),ly,lx),nvals);
    title('true label');
    hold on; axes();
    subplot(M,2,1+ M); miximshow(reshape(label_pred(:),ly,lx),nvals);
    title('predicted label');
    % subplot(M,3,1+ 2*M);  imshow(reshape(label_pred(:),ly,lx)); % -1 to the label just in case they were computed with max()
end
end