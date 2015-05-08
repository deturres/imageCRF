function [D_euclidean,D_quasiEuclidean] = distance_map(I)

method = 'euclidean';
level = graythresh(I);
bwimage = im2bw(I, level);

imshow(I), figure, imshow(bwimage)

D_euclidean = bwdist(~bwimage,method);
method = 'quasi-euclidean';
D_quasiEuclidean = bwdist(~bwimage,method);

figure, subimage(mat2gray(D_euclidean)), title('Euclidean')
% hold on, imcontour(D_euclidean)

% figure, subplot(1,2,1), subimage(mat2gray(D_euclidean)), title('Euclidean')
% hold on, imcontour(D_euclidean)
% subplot(1,2,2), subimage(mat2gray(D_quasiEuclidean)), title('Quasi-Euclidean')
% hold on, imcontour(D_quasiEuclidean)

end