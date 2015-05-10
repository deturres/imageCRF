function [D_euclidean,D_euclidean_compl] = distance_map(I)

% comment the function
% 

% inputs:
%   I           - an image in gray scale (should be size of the sidewalkdetector_feature ly x lx)
% output:
%   D_euclidean        - the distance map (unary feature) 
%   D_euclidean_compl  - the complemtary distance map (unary feature)

% define the level according to which compute the binaty tranformation
level = 0.935; % graythresh(I)
% compute the binary image transformation
bwimage = im2bw(I, level);
figure, subplot(1,3,1), imshow(I),
subplot(1,3,2), imshow(bwimage),
subplot(1,3,3), imshow(~bwimage),

% define the  method to be used for the distance map
method = 'euclidean';
% compute the distance transform and the skeleton og the binary image
D_euclidean = bwdist(bwimage,method);
D_euclidean_compl = bwdist(~bwimage,method);
% D_skel = bwmorph(bwimage,'skel', 1);
% D_skel_compl = bwmorph(~bwimage,'skel', 1);

figure, subplot(1,2,1), subimage(mat2gray(D_euclidean)), title('Euclidean distance trasform')
% hold on, imcontour(D_euclidean)
subplot(1,2,2), subimage(mat2gray(D_euclidean_compl)), title('Euclidean distance trasform complementary')
% hold on, imcontour(D_euclidean_compl)
% subplot(2,2,3), subimage(mat2gray(D_skel)), title('Skeleton')
% subplot(2,2,4), subimage(mat2gray(D_skel_compl)), title('Skeleton complementary')

% % compute the gray-weighted distance transform of a greyyscale image
% mask = I<0.2;
% size(I)
% size(mask)
% D_euclidean_gray = graydist(I,mask);
% figure, subimage(mat2gray(D_euclidean_gray)), title('Euclidean gray-weighted distance trasform')
% % hold on, imcontour(D_euclidean_grey)

end