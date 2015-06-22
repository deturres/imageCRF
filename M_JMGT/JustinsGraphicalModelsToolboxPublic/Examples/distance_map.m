function [D_euclidean,D_euclidean_compl,D_euclidean_signed] = distance_map(I)

% comment the function
% inputs:
%   I           - an image in gray scale (should be size of the sidewalkdetector_feature ly x lx)
% output:
%   D_euclidean        - the distance map (unary feature) 
%   D_euclidean_compl  - the complemtary distance map (unary feature)

% define the level according to which compute the binaty tranformation
level = 0.9; % graythresh(I)
% compute the binary image transformation
bwimage = im2bw(I, level);
figure('Name','Original image and corresponding binary image','NumberTitle','off'), 
subplot(1,2,1), imshow(I),
subplot(1,2,2), imshow(bwimage),
% subplot(1,3,3), imshow(~bwimage),

% compute the distance transform and the skeleton og the binary image
D_euclidean = double(bwdist(bwimage,'euclidean'));
D_euclidean_compl = double(bwdist(~bwimage,'euclidean'));
% D_skel = bwmorph(bwimage,'skel', 1);
% D_skel_compl = bwmorph(~bwimage,'skel', 1);


% create a map of points inside the closed surface
M = imfill(bwimage, 'holes');
figure('Name','imfill output','NumberTitle','off'),  imshow(M);
% generate the signed distance transform
D_euclidean_signed = D_euclidean_compl;
D_euclidean_signed(M) = -D_euclidean_compl(M);
% D_euclidean_signed = -D_euclidean_signed;

% plot preview of the tranforms
% figure, subplot(1,2,1), subimage(mat2gray(D_euclidean_signed)), title('Euclidean distance trasform signed')
% hold on, imcontour(D_euclidean_signed)
% subplot(1,2,2), subimage(mat2gray(D_euclidean_compl)), title('Euclidean distance trasform complementary')
% hold on, imcontour(D_euclidean_compl)

% subplot(2,2,3), subimage(mat2gray(D_skel)), title('Skeleton')
% subplot(2,2,4), subimage(mat2gray(D_skel_compl)), title('Skeleton complementary')

%% try another transform, it gives same result
% compute the gray-weighted distance transform of a greyyscale image
% mask = I<0.2;
% size(I)
% size(mask)
% D_euclidean_gray = graydist(I,mask);
% figure, subimage(mat2gray(D_euclidean_gray)), title('Euclidean gray-weighted distance trasform')
% % hold on, imcontour(D_euclidean_grey)

end