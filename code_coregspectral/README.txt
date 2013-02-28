Paper: Coregularized multiview spectral clustering, NIPS, 2011.
=================================================================

Code
========
run_synth.m:  example code showing the usage of the two multiview approaches proposed in the paper

baseline_spectral.m:  single view spectral clustering using normalized Laplacian (symmetric)
spectral_mindis.m:  Two view spectral clustering of De Sa et al
spectral_mindis_onkernel.m:  Two view spectral clustering of De Sa et al (on similarity matrices)
spectral_pairwise.m:  Pairwise coregularized spectral clustering
spectral_pairwise_multiview.m:  Pairwise coregularized spectral clustering (more than two views)
spectral_pairwise_multiview_onkernel.m:  Pairwise coregularized spectral clustering (more than two views) on similarity matrices
spectral_centroid.m:  Centroid based coregularized spectral clustering
spectral_centroid_multiview.m:  Centroid based coregularized spectral clustering (more than two views)
spectral_centroid_multiview_onkernel.m: Centroid based coregularized spectral clustering (more than two views) on similarity matrices


Note
=======
1. All the codes use Gaussian kernel for construction of similarity matrix.

2. Pairwise spectral approach reports the clsutering performance on the first view in each iteration. If you have prior knowledge about which view is more informative about the clusters, use that as the first view.


E.g.

X{1} = X1; X{2} = X2;
sigma(1) = sigma1; sigma(2) = sigma2;
lambda = 0.5; numiter = 5; numClust = 2; num_views = 2;
[F P R nmi avgent AR] = spectral_pairwise_multview(X,num_views,numClust,sigma,lambda,truth,numiter);

The performance measures are reported on X{1} in each iteration in the above function call.

If you want the performance to be reported on X{2} instead, please do as follows:
X{1} = X2; X{2} = X1;
sigma(1) = sigma2; sigma(2) = sigma1;
lambda = 0.5; numiter = 5; numClust = 2; num_views = 2;
[F P R nmi avgent AR] = spectral_pairwise_multview(X,num_views,numClust,sigma,lambda,truth,numiter);


3. Performance can be sensitive to the Lambda value chosen.

 