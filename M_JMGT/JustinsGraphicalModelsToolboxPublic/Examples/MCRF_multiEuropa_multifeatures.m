function MCRF_multiEuropa_multifeatures(path_name)

%% just temporarly done to save the correct angles features
A = double(imread(('./Dataset/europaData/entire_log_new/res10/01_mapImage0.1_AfromCode.png')))/255;
imgangle = rgb2gray(A);
% figure('Name','Loading angleWRTroadsInVicinity...','NumberTitle','off'); imshow(imgangle);

% load as data
angdir = [ path_name ];
ang_names = dir([angdir '*0.1_A.png.dat']); %% the dat file is built up row_wise
ang = importdata([angdir ang_names.name]);
ang_img = imcomplement(reshape(ang, size(imgangle,2),size(imgangle,1)));
ang_img = ang_img';
figure('Name','Loading angleWRTroadsInVicinity load from dat file...','NumberTitle','off'); imshow(ang_img);
imwrite(ang_img,'./Dataset/europaData/entire_log_new/res10/train/01_mapImage0.1_AREAL.png', 'png');


%% load the data and computing labels and features map
imdir = [ path_name 'train/portion/new/']; % Valid for entire log new dataset(inside trains) or portion(inside trains/portion)
im_names = dir([imdir '*0.1.png']);
imstepdir = [ path_name 'train/portion/new/']; % Valid for entire log new dataset(inside trains) or portion(inside trains/portion)
imstep_names = dir([imstepdir '*0.1_step.png']);
imweightdir = [ path_name 'train/portion/new/']; % Loading the weights images correspondent to the original features
imweight_names = dir([imweightdir '*0.1_W.png']);
imangledir = [ path_name 'train/portion/new/']; % Loading the weights images correspondent to the original features
imangle_names = dir([imangledir '*0.1_A.png']);
labdir = [ path_name 'labels/portion/new/']; % same name for entire_log_new or portion
lab_names = dir([labdir '*multi4AREA_GT.png']);

% parameters of the problem
N     = length(im_names);  % size of training images
rho   = .5; % TRW edge appearance probability
nvals = 4; % curb/sidewalkWall/buildings/background (0 in case unlabeled)
rez    = .5; % how much resolution percentage to use
cmap = [1 1 1; 1 0 0; 0 0 1; 0 1 0]; % to represent the label

fprintf('loading data and computing feature maps...\n');
% load true label and input images from europa2_sidewalktetector
ims = cell(1,N);
ims2 = cell(1,N);
ims3 = cell(1,N);
labelsRGB = cell(1,N);
labels0 = cell(1,N);
labels = cell(1,N);
feats = cell(1,N);
wimgs = cell(1,N);
efeats = cell(1,N);
for n=1:N
    
    % load input images (three main features)
    I = double(imread(([imdir im_names(n).name])))/255;    
    img = rgb2gray(I);
    ims{n}  = img; % input images, first feature (heightGradientChange)
%     figure('Name','Loading input feature 1 [heightGradientChange]...','NumberTitle','off'); imshow(ims{n});
    Istep = double(imread(([imstepdir imstep_names(n).name])))/255;    
    imgstep = rgb2gray(Istep);
    ims2{n}  = imgstep; % input images, second feature (stepHeightInVicinity)
%     figure('Name','Loading input feature 2 [stepHeightInVicinity]...','NumberTitle','off'); imshow(ims2{n});
    % load angles
    A = double(imread(([imangledir imangle_names(n).name])))/255;
%     imgangle = rgb2gray(A);
    ims3{n}  = A; % input images, third feature (angle wrt roads in vicinity)
%     figure('Name','Loading input feature 3 [angleWRTroadsInVicinity]...','NumberTitle','off'); imshow(ims3{n});

    % load weights
    W = double(imread(([imweightdir imweight_names(n).name])))/255;
    wimg = rgb2gray(W);
    wimgs{n}  = wimg; % weights images to be used as filter or unary features
%     figure('Name','Loading weights...','NumberTitle','off'); imshow(wimgs{n});

    % load labels
    L = double(imread(([labdir lab_names(n).name])))/255;
    labelsRGB{n}  = L; % true label GT
%     figure('Name','Loading label RGB...','NumberTitle','off'); imshow(labelsRGB{n});
    
end

%% 
% The labels representation consists on values from  1 to nvals, with 0 for unlabeled
for n=1:N
    fprintf('new image\n');
    % reduce resolution for speed (mostly in case we use the different GRIDMAPS images)
    % DO NOT REDUCE when using small portion of the entire log if it's 0.5
    ims{n}  = imresize(ims{n},rez,'bilinear');
    ims2{n} = imresize(ims2{n},rez,'bilinear');
    ims3{n} = imresize(ims3{n},rez,'bilinear');

    wimgs{n} = imresize(wimgs{n},rez,'bilinear');

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
    
    % compute the distance transform image
    [D_euclidean,D_euclidean_compl,D_euclidean_signed] = distance_map(ims{n});
    [hor_efeats_ij ver_efeats_ij] = evaluate_pca(ims{n});
    % normalizing the distance map value
    %D_euclidean = D_euclidean(:)/norm(D_euclidean(:));
    D_euclidean_compl = D_euclidean_compl(:)/norm(D_euclidean_compl(:));
%     D_euclidean_signed = D_euclidean_signed(:)/norm(D_euclidean_signed(:));

    feats{n}  = [ims{n}(:) ims2{n}(:) ims3{n}(:) D_euclidean_signed(:) hor_efeats_ij(:) ver_efeats_ij(:) 1+0*labels{n}(:)]; %D_euclidean_compl(:)
    fprintf('features computed\n');    
end
    fprintf('end dataset\n');

%% creating the models
% the images come in slightly different sizes, so we need to make many models
% use a "hashing" strategy to not rebuild.  Start with empty giant array

% very big model_hash in case of entire_log_new
%   model_hash = repmat({[]},1000,1150);

% smaller model_hash in case of small portion of entire_log_new
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
K = N;
[who_train who_test] = kfold_sets(N,K,k)

% decide which features to use for the training
using_feats = cell(1,N);
for n=1:N
    using_feats{n} = feats{n}(:,1:5);
end

ims_train     = ims(who_train);
feats_train   = using_feats(who_train);
efeats_train  = []; % efeats(who_train);
labels_train  = labels(who_train);
labels0_train = labels0(who_train);
models_train  = models(who_train);
% imshow(ims_train{1})

ims_test     = ims(who_test);
feats_test   = using_feats(who_test);
efeats_test  = []; % efeats(who_test);
labels_test  = labels(who_test);
labels0_test = labels0(who_test);
models_test  = models(who_test);
% imshow(ims_test{1})

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

%% plot features
fprintf('Plotting features computed\n');    
for n=1:N
    
    [ly lx] = size(labels{n});
    feat1 = reshape(using_feats{n}(:,1),ly,lx);
    feat2 = reshape(using_feats{n}(:,2),ly,lx);
    feat3 = reshape(using_feats{n}(:,3),ly,lx);
    dist_map = reshape(using_feats{n}(:,4),ly,lx);
%     dist_map = dist_map(:)/norm(dist_map(:));
    pca_hor = reshape(using_feats{n}(:,5),ly,lx);
%     pca_ver = reshape(using_feats{n}(:,6),ly,lx);
    figure('Name', 'Features used'), 
    subplot(2,3,1), subimage(mat2gray(feat1)), title('feat1 gradient')
    subplot(2,3,2), subimage(mat2gray(feat2)), title('feat2 steps')
    subplot(2,3,3), subimage(mat2gray(feat3)), title('feat3 angles')

    subplot(2,3,4), subimage(mat2gray(dist_map)), title('Distance map'), % hold on, imcontour(dist_map);
    subplot(2,3,5), subimage(mat2gray(pca_hor)), title('First pca component')
%     subplot(2,3,6), subimage(mat2gray(pca_ver)), title('Second pca component')
end


%% %-----------------training----------------%

fprintf('training the model (this is slow!)...\n')
% We pick a string to specify the loss and inference method. In this case, we choose truncated fitting with the clique logistic loss based on multithreaded TRW with five/ten iterations.
% Other options include 'pert_ul_trw_1e5' (perturbation, univariate logistic loss, TRW, threshold of 1e-5),
% 'em_mnf_1e5' (Surrogate Expectation-Maximization based on mean-field with a threshold of 1e-5 (simplifies to surrogate likelihood with no hidden variables),
% or 'trunc_em_trwpll_10' (Truncated surrogate EM based on multithreaded TRW with 10 iterations).
loss_spec = 'trunc_cl_trwpll_10'; % 'trunc_cl_trw_5'  parameter learning and inference

% some parameters for the training optimization
crf_type  = 'linear_linear'; 
options.derivative_check = 'off';
% options.viz         = @viz; % function for visualization
options.rho         = rho;
options.print_times = 1; % since this is so slow, print stuff to screen
options.nvals       = nvals;

options.gradual     = 1; % use gradual fitting
options.maxiter     = 3000;
options.rho         = rho;
options.reg         = 1e-4;
options.opt_display = 0;

% figure('Name','Training...','NumberTitle','off');

% we actually optimize
% This prints a visualization while running , using the viz function above, if enabled.
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
    
    
    
%  to be included inside a smoothing_result function
    
    % max value (taking the corresponding index-->class)
    [~,label_pred] = max(b_i,[],1);
    ratio_confidence = 0.38;
    for b = 1:size(label_pred,2)
        if(label_pred(b)~=1 && label_pred(b)~=4)
            b_i_ratio = b_i(3,b)/b_i(2,b); % ratio between sidewalk belief and street belief
            if(b_i_ratio>ratio_confidence)
                label_pred(b) = 3;
            end
        end
    end    
    
    % loading weights to give more confidence to label=1(background if the
    % weights features has value = 1 and the ;learning didnt recognized was
    % buildings or background (no information from the sidewalk detector)

    weights = wimgs(who_test);
    w = weights{n}(:);
    for i = 1:size(w)
        if(w(i)>0.99 && label_pred(i)~=4)
          label_pred(i) = 1;
        end
    end
    
% untill here 



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
    subplot(1,2,1); imshow(reshape(feats_test{n}(:,1),ly,lx));
    subplot(1,2,2); imshow(reshape(feats_test{n}(:,2),ly,lx));

    M = length(feats_test);
    % visualizing final predicted marginal and label
    figure('Name','Testing..predicted marginalbelief','NumberTitle','off');
    subplot(2,2,1); imshow(reshape(b_i(1,:),ly,lx));
    title('(class 1-background)');
    subplot(2,2,2); imshow(reshape(b_i(2,:),ly,lx));
    title('(class 2-street)'); % street
    subplot(2,2,3); imshow(reshape(b_i(3,:),ly,lx));
    title('(class 3-sidewalk)'); % sidewalk
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
end
end