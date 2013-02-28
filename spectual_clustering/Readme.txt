% ===========================
% Bachelor Thesis
% Author : Ingo Bürk
% Year   : 2011/2012
% Contact: admin@airblader.de
% ===========================

% ===========================================================
% GENERAL INFORMATION
% ===========================================================

The Matlab files provided can create sample data, different
kinds of similarity graphs and perform fast and efficient
spectral clustering algorithms.
If you need help, please read this file first and try typing
'help [Filename]' to get information. If there are still
questions left, I'll be happy to help you upon contacting
me via email.

% ===========================================================
% TECHNICAL INFORMATION
% ===========================================================

If you load and use your own data (adjacency matrix), please
keep in mind that using a sparse matrix will reduce memory
drastically.
All methods in these files work with sparse matrices and
therefore assume a sparse structure.

If your data is not sparse per se, consider using a
similarity graph that will give it a sparse structure.

I recommend using either an epsilon or mutual k-Nearest
neighbors similarity graph, combined with the spectral
clustering algorithm according to Shi and Malik (2000).


% ===========================================================
% REFERENCES
% ===========================================================

This thesis and the resulting work is based on

- Ulrike von Luxburg, "A Tutorial on Spectral Clustering", 
  Statistics and Computing 17 (4), 2007