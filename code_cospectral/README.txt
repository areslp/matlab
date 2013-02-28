Paper: A Co-training approach for multiview spectral clustering, ICML, 2011.
=============================================================================

Code
========
run_synth.m:  example code showing the usage of the the multiview spectral clustering approach proposed in the paper

baseline_spectral_onkernel.m:  single view spectral clustering using normalized Laplacian (symmetric)
spectral_cotraining.m: co-trained spectral clustering proposed in the paper

Note
=======
1. All spectral clustering codes provided here use Gaussian kernel for construction of similarity matrix.
 
2. Number of eigenvectors to project on (decided by projev argument in the MATLAB code), can be varied increasing from 1. We have observed that a value in the range of 1 to 1.5 gives reasonable results. 