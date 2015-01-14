Please also look at README.txt and licence files. If in doubt contact the authors.
This software is modified from the MRF 1.9/2.0 package available at :
http://vision.middlebury.edu/MRF/code/MRF2.0.zip

The Graph cuts code from the original MRF package is modified slightly to work
for the image segmentation of 2D and 3D CRFs and to run on energies that can be 
negative.
Some changes have been made so Graph cuts runs on 64bit mac os x 10.7 besides 32bit Fedora release 18 / 64bit Fedora release 19

Added code for Max Product runs on general graphs by Chetan Bhole
Added code for Variational Message Passing (meanfield) , Syn etc versions by Jerod Weinmann
Added code for Variational Message Passing for general graph by Chetan Bhole
Added code for dual decomposition using trees and unimodular graphs on binary label problems by Chetan Bhole

If you use this software, read README.txt for citation as well. Depending on which algorithm you 
use, you should cite the appropriate paper.

[1] QPBO : please look at QPBO-v1.3.src/QPBO.h for citation details.

[2] Variational Message Passing (MeanField), SynMeanField, ASyncMeanField, SumProd, SyncSumProd, ASyncSumProd, DenseSumProd, DenseSyncSumProd,  SparseMeanField, SparseAsyncMeanField, SparseAsyncSumProd for 2D grid graphs : 
J. Weinman, L. Tran, C. Pal
Efficiently Learning Random Fields for Stereo Vision with Sparse Message Passing. In proc. European Conf. on Computer Vision (ECCV), Springer-Verlag LNCS, vol. 1, pp. 617-630, 2008 

[3] Variational Message Passing for general graph (non 2D grid) : 
C. Bhole, C. Pal, D. Rim, A. Wismüller
3D segmentation of abdominal CT imagery with graphical models, conditional random fields and learning 
Machine Vision and Applications, pp. 1-25, 2013. 
and 
J. Weinman, L. Tran, C. Pal
Efficiently Learning Random Fields for Stereo Vision with Sparse Message Passing. In proc. European Conf. on Computer Vision (ECCV), Springer-Verlag LNCS, vol. 1, pp. 617-630, 2008 

[4] Dual decomposition using trees and unimodular graphs, Flipper graphs on binary label problems and MaxProdBPTree : 
C. Bhole, J. Domke, D. Gildea
Approximate inference using unimodular graphs in dual decomposition,
In Optimization in Machine Learning, NIPS Workshop, 2013

[5] Max Product (MaxProdBP) runs on general graphs (non 2D grid) :
C. Bhole, C. Pal, D. Rim, A. Wismüller
3D segmentation of abdominal CT imagery with graphical models, conditional random fields and learning 
Machine Vision and Applications, pp. 1-25, 2013. 
and 
M. F. Tappen and W. T. Freeman.
Comparison of Graph Cuts with Belief Propagation for Stereo, using Identical MRF Parameters.
In Proceedings of the Ninth IEEE International Conference on Computer Vision (ICCV),
pages 900-907, 2003

