#ifndef __COMP_MESHLPMATRIX_H__
#define __COMP_MESHLPMATRIX_H__

#include "tmesh.h"
void generate_sym_meshlp_matrix(TMesh& mesh, double h, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV, vector<double>& AAV);
void generate_sym_meshlp_matrix_geod(TMesh& mesh, double h, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV, vector<double>& AAV);

void generate_meshlp_matrix(TMesh& mesh, double h, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV);
void generate_meshlp_matrix_geod(TMesh& mesh, double h, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV);

void generate_meshlp_matrix_adp(TMesh& mesh, double hs, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV);
void generate_meshlp_matrix_adp_geod(TMesh& mesh, double hs, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV);


void generate_Xu_Meyer_laplace_matrix(TMesh& mesh,  vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV, vector<double>& DDV);


//void compute_one2part_Euclidean_vdist(unsigned int vid_start, TMesh& mesh, vector<pair<unsigned int, double> >& vgdists, double maxdist);

#endif //__COMP_MESHLPMATRIX_H__
