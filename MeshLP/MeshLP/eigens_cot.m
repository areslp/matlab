nev = 300;
prefix='bighand';

filename=sprintf('./%s.off', prefix);

[W A]=cotlp_matrix(filename);

Am = sparse([1:length(A)], [1:length(A)], A);

[evecs evals] = eigs(W, Am, nev, -1e-5);
evals = diag(evals);

filename=sprintf('./%s_cot', prefix);
save(filename, 'W', 'A', 'evecs', 'evals');
