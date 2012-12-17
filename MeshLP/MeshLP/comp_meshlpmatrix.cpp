#include <mex.h>
#include "comp_meshlpmatrix.h"
#include "geodesics/geodesic_algorithm_exact.h"

void compute_one2part_Euclidean_vdist(unsigned int vid_start, TMesh& mesh, vector<pair<unsigned int, double> >& vgdists, double maxdist);
void compute_one2part_Geodesic_vdist(unsigned int vid_start, geodesic::Mesh& geod_mesh, geodesic::GeodesicAlgorithmExact& algorithm, vector<pair<unsigned int, double> >& vgdists, double maxdist);


void generate_sym_meshlp_matrix(TMesh& mesh, double h, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV, vector<double>& AAV)
{
	unsigned int nv = mesh.v_count();
	unsigned int nf = mesh.f_count();

	//compute area for each vertexs
	AAV.clear();
	AAV.resize(nv, 0);
		
	unsigned int vid0, vid1, vid2;
	double farea;
	for(unsigned int f = 0; f < nf; f ++){
		vid0 = mesh.facet(f).vert(0);
		vid1 = mesh.facet(f).vert(1);
		vid2 = mesh.facet(f).vert(2);
		VECTOR3 vv = cross( mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord(), 
								  mesh.vertex(vid2).coord() - mesh.vertex(vid0).coord() );
		farea = norm(vv) / 2.0;

		AAV[vid0] += farea / 3;
		AAV[vid1] += farea / 3;
		AAV[vid2] += farea / 3;
	}
	
	//compute the laplacian matrix

	vector<pair<unsigned int, double> >vgdists;
	vector<double> totalweight;
	totalweight.resize(nv, 0);

	double hh = h * h;
	for(unsigned int i = 0; i < nv; i ++){
		//cout<<"i: "<<i<<"\t \r";
		mexPrintf("i: %d\r", i);

		vgdists.clear();		
		compute_one2part_Euclidean_vdist(i, mesh, vgdists, h * rho);
		
		for(unsigned int j = 0; j < vgdists.size(); j ++){
			unsigned int vid = vgdists[j].first;
			if( vid <= i ){
				continue;
			}

			double weight = exp(-vgdists[j].second * vgdists[j].second / hh) * ( 4.0 / (M_PI * hh * hh) );
			//cout<<"vid: "<<vid<<" dist: "<<vgdists[j].second<<" weight: "<<weight<<" vareas: "<<vareas[vid]<<endl;

			weight *=  AAV[vid] * AAV[i];

			IIV.push_back(i + 1);
			JJV.push_back(vid + 1);
			SSV.push_back(weight);

			IIV.push_back(vid + 1);
			JJV.push_back(i + 1);
			SSV.push_back(weight);

			totalweight[i] -= weight;
			totalweight[vid] -= weight;
		}
	}
	for(unsigned int i = 0; i < nv; i ++){
		IIV.push_back(i + 1);
		JJV.push_back(i + 1);
		SSV.push_back(totalweight[i]);
	}
	mexPrintf("\n");
}

void generate_sym_meshlp_matrix_geod(TMesh& mesh, double h, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV, vector<double>& AAV)
{
	unsigned int nv = mesh.v_count();
	unsigned int nf = mesh.f_count();

	//compute area for each vertexs
	AAV.clear();
	AAV.resize(nv, 0);
		
	unsigned int vid0, vid1, vid2;
	double farea;
	for(unsigned int f = 0; f < nf; f ++){
		vid0 = mesh.facet(f).vert(0);
		vid1 = mesh.facet(f).vert(1);
		vid2 = mesh.facet(f).vert(2);
		VECTOR3 vv = cross( mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord(), 
								  mesh.vertex(vid2).coord() - mesh.vertex(vid0).coord() );
		farea = norm(vv) / 2.0;

		AAV[vid0] += farea / 3;
		AAV[vid1] += farea / 3;
		AAV[vid2] += farea / 3;
	}
	
	//for geodesic computation
	vector<double> points;   
	vector<unsigned int> faces;
	for(unsigned int i = 0; i < nv; i ++){
		points.push_back((mesh.vertex(i).coord())[0]);
		points.push_back((mesh.vertex(i).coord())[1]);
		points.push_back((mesh.vertex(i).coord())[2]);
	}
	for(unsigned int i = 0; i < nf; i ++){
		faces.push_back(mesh.facet(i).vert(0));
		faces.push_back(mesh.facet(i).vert(1));
		faces.push_back(mesh.facet(i).vert(2));
	}

	geodesic::Mesh geod_mesh;
	geod_mesh.initialize_mesh_data(points, faces);    //create internal mesh data structure including edges
	geodesic::GeodesicAlgorithmExact algorithm(&geod_mesh); //create exact algorithm for the mesh

	
	//compute the laplacian matrix
	vector<pair<unsigned int, double> >vgdists;
	vector<double> totalweight;
	totalweight.resize(nv, 0);

	double hh = h * h;
	for(unsigned int i = 0; i < nv; i ++){
		//cout<<"i: "<<i<<"\t \r";
		mexPrintf("i: %d\r", i);

		vgdists.clear();		
		compute_one2part_Geodesic_vdist(i, geod_mesh, algorithm, vgdists, h * rho);
		//compute_one2part_Euclidean_vdist(i, mesh, vgdists, h * rho);

		for(unsigned int j = 0; j < vgdists.size(); j ++){
			unsigned int vid = vgdists[j].first;
			if( vid <= i ){
				continue;
			}

			double weight = exp(-vgdists[j].second * vgdists[j].second / hh) * ( 4.0 / (M_PI * hh * hh) );
			//cout<<"vid: "<<vid<<" dist: "<<vgdists[j].second<<" weight: "<<weight<<" vareas: "<<vareas[vid]<<endl;

			weight *=  AAV[vid] * AAV[i];

			IIV.push_back(i + 1);
			JJV.push_back(vid + 1);
			SSV.push_back(weight);

			IIV.push_back(vid + 1);
			JJV.push_back(i + 1);
			SSV.push_back(weight);

			totalweight[i] -= weight;
			totalweight[vid] -= weight;
		}
	}
	for(unsigned int i = 0; i < nv; i ++){
		IIV.push_back(i + 1);
		JJV.push_back(i + 1);
		SSV.push_back(totalweight[i]);
	}
	mexPrintf("\n");
}


void generate_meshlp_matrix(TMesh& mesh, double h, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV)
{
	unsigned int nv = mesh.v_count();
	unsigned int nf = mesh.f_count();

	//compute area for each vertexs
	vector<double> vareas;
	vareas.resize(nv, 0);
		
	unsigned int vid0, vid1, vid2;
	double farea;
	for(unsigned int f = 0; f < nf; f ++){
		vid0 = mesh.facet(f).vert(0);
		vid1 = mesh.facet(f).vert(1);
		vid2 = mesh.facet(f).vert(2);
		VECTOR3 vv = cross( mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord(), 
								  mesh.vertex(vid2).coord() - mesh.vertex(vid0).coord() );
		farea = norm(vv) / 2.0;

		vareas[vid0] += farea / 3;
		vareas[vid1] += farea / 3;
		vareas[vid2] += farea / 3;
	}
	
	//compute the laplacian matrix

	vector<pair<unsigned int, double> >vgdists;
	double totalweight;

	double hh = h * h;
	for(unsigned int i = 0; i < nv; i ++){
		//cout<<"i: "<<i<<"\t \r";
		mexPrintf("i: %d\r", i);

		vgdists.clear();		
		compute_one2part_Euclidean_vdist(i, mesh, vgdists, h * rho);
		
		totalweight = 0;
		for(unsigned int j = 0; j < vgdists.size(); j ++){
			unsigned int vid = vgdists[j].first;
			if( vid == i ){
				continue;
			}

			double weight = exp(-vgdists[j].second * vgdists[j].second / hh) * ( 4.0 / (M_PI * hh * hh) );
			//cout<<"vid: "<<vid<<" dist: "<<vgdists[j].second<<" weight: "<<weight<<" vareas: "<<vareas[vid]<<endl;

			weight *=  vareas[vid];

			IIV.push_back(i + 1);
			JJV.push_back(vid + 1);
			SSV.push_back(weight);

			totalweight -= weight;
		}
		
		IIV.push_back(i + 1);
		JJV.push_back(i + 1);
		SSV.push_back(totalweight);
	}
	
	mexPrintf("\n");
}


void generate_meshlp_matrix_geod(TMesh& mesh, double h, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV)
{
	unsigned int nv = mesh.v_count();
	unsigned int nf = mesh.f_count();

	//compute area for each vertexs
	vector<double> vareas;
	vareas.resize(nv, 0);
		
	unsigned int vid0, vid1, vid2;
	double farea;
	for(unsigned int f = 0; f < nf; f ++){
		vid0 = mesh.facet(f).vert(0);
		vid1 = mesh.facet(f).vert(1);
		vid2 = mesh.facet(f).vert(2);
		VECTOR3 vv = cross( mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord(), 
								  mesh.vertex(vid2).coord() - mesh.vertex(vid0).coord() );
		farea = norm(vv) / 2.0;

		vareas[vid0] += farea / 3;
		vareas[vid1] += farea / 3;
		vareas[vid2] += farea / 3;
	}
	
	//for geodesic computation
	vector<double> points;   
	vector<unsigned int> faces;
	for(unsigned int i = 0; i < nv; i ++){
		points.push_back((mesh.vertex(i).coord())[0]);
		points.push_back((mesh.vertex(i).coord())[1]);
		points.push_back((mesh.vertex(i).coord())[2]);
	}
	for(unsigned int i = 0; i < nf; i ++){
		faces.push_back(mesh.facet(i).vert(0));
		faces.push_back(mesh.facet(i).vert(1));
		faces.push_back(mesh.facet(i).vert(2));
	}

	geodesic::Mesh geod_mesh;
	geod_mesh.initialize_mesh_data(points, faces);    //create internal mesh data structure including edges
	geodesic::GeodesicAlgorithmExact algorithm(&geod_mesh); //create exact algorithm for the mesh

	//compute the laplacian matrix
	vector<pair<unsigned int, double> >vgdists;
	double totalweight;

	double hh = h * h;
	for(unsigned int i = 0; i < nv; i ++){
		//cout<<"i: "<<i<<"\t \r";
		mexPrintf("i: %d\r", i);

		vgdists.clear();		
		compute_one2part_Geodesic_vdist(i, geod_mesh, algorithm, vgdists, h * rho);
		//compute_one2part_Euclidean_vdist(i, mesh, vgdists, h * rho);
		
		totalweight = 0;
		for(unsigned int j = 0; j < vgdists.size(); j ++){
			unsigned int vid = vgdists[j].first;
			if( vid == i ){
				continue;
			}

			double weight = exp(-vgdists[j].second * vgdists[j].second / hh) * ( 4.0 / (M_PI * hh * hh) );
			//cout<<"vid: "<<vid<<" dist: "<<vgdists[j].second<<" weight: "<<weight<<" vareas: "<<vareas[vid]<<endl;

			weight *=  vareas[vid];

			IIV.push_back(i + 1);
			JJV.push_back(vid + 1);
			SSV.push_back(weight);

			totalweight -= weight;
		}
		
		IIV.push_back(i + 1);
		JJV.push_back(i + 1);
		SSV.push_back(totalweight);
	}
	mexPrintf("\n");
}

void generate_meshlp_matrix_adp(TMesh& mesh, double hs, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV)
{
	unsigned int nv = mesh.v_count();
	unsigned int nf = mesh.f_count();

	//compute area for each vertexs
	vector<double> vareas;
	vareas.resize(nv, 0);
		
	unsigned int vid0, vid1, vid2;
	double farea;
	for(unsigned int f = 0; f < nf; f ++){
		vid0 = mesh.facet(f).vert(0);
		vid1 = mesh.facet(f).vert(1);
		vid2 = mesh.facet(f).vert(2);
		VECTOR3 vv = cross( mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord(), 
								  mesh.vertex(vid2).coord() - mesh.vertex(vid0).coord() );
		farea = norm(vv) / 2.0;

		vareas[vid0] += farea / 3;
		vareas[vid1] += farea / 3;
		vareas[vid2] += farea / 3;
	}


	//compute the laplacian matrix

	vector<pair<unsigned int, double> >vgdists;
	double totalweight;

	for(unsigned int i = 0; i < nv; i ++){
		//cout<<"i: "<<i<<"\t \r";
		//mexPrintf("i: %d\r", i);

		double h = 0;
		for(unsigned int j = 0; j < mesh.vertex(i).n_verts(); j ++){
    		unsigned int k = mesh.vertex(i).vert(j);
    		h +=  sqrt( fabs(dot(mesh.vertex(i).coord() - mesh.vertex(k).coord(),
         		      			 mesh.vertex(i).coord() - mesh.vertex(k).coord())) );
		}

		if(mesh.vertex(i).n_verts() > 0){
			h = hs * h / mesh.vertex(i).n_verts();
		}

		mexPrintf("i: %d h: %f\r", i, h);

		double hh = h * h;

		vgdists.clear();		
		compute_one2part_Euclidean_vdist(i, mesh, vgdists, h * rho);
		
		totalweight = 0;
		for(unsigned int j = 0; j < vgdists.size(); j ++){
			unsigned int vid = vgdists[j].first;
			if( vid == i ){
				continue;
			}

			double weight = exp(-vgdists[j].second * vgdists[j].second / hh) * ( 4.0 / (M_PI * hh * hh) );
			//cout<<"vid: "<<vid<<" dist: "<<vgdists[j].second<<" weight: "<<weight<<" vareas: "<<vareas[vid]<<endl;

			weight *=  vareas[vid];

			IIV.push_back(i + 1);
			JJV.push_back(vid + 1);
			SSV.push_back(weight);

			totalweight -= weight;
		}
		
		IIV.push_back(i + 1);
		JJV.push_back(i + 1);
		SSV.push_back(totalweight);
	}
	
	mexPrintf("\n");
}

void generate_meshlp_matrix_adp_geod(TMesh& mesh, double hs, double rho, vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV)
{
	unsigned int nv = mesh.v_count();
	unsigned int nf = mesh.f_count();

	//compute area for each vertexs
	vector<double> vareas;
	vareas.resize(nv, 0);
		
	unsigned int vid0, vid1, vid2;
	double farea;
	for(unsigned int f = 0; f < nf; f ++){
		vid0 = mesh.facet(f).vert(0);
		vid1 = mesh.facet(f).vert(1);
		vid2 = mesh.facet(f).vert(2);
		VECTOR3 vv = cross( mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord(), 
								  mesh.vertex(vid2).coord() - mesh.vertex(vid0).coord() );
		farea = norm(vv) / 2.0;

		vareas[vid0] += farea / 3;
		vareas[vid1] += farea / 3;
		vareas[vid2] += farea / 3;
	}

	//for geodesic computation
	vector<double> points;   
	vector<unsigned int> faces;
	for(unsigned int i = 0; i < nv; i ++){
		points.push_back((mesh.vertex(i).coord())[0]);
		points.push_back((mesh.vertex(i).coord())[1]);
		points.push_back((mesh.vertex(i).coord())[2]);
	}
	for(unsigned int i = 0; i < nf; i ++){
		faces.push_back(mesh.facet(i).vert(0));
		faces.push_back(mesh.facet(i).vert(1));
		faces.push_back(mesh.facet(i).vert(2));
	}

	//for geodesics computation
	geodesic::Mesh geod_mesh;
	geod_mesh.initialize_mesh_data(points, faces);    //create internal mesh data structure including edges
	geodesic::GeodesicAlgorithmExact algorithm(&geod_mesh); //create exact algorithm for the mesh

	//compute the laplacian matrix
	vector<pair<unsigned int, double> >vgdists;
	double totalweight;

	for(unsigned int i = 0; i < nv; i ++){
		//cout<<"i: "<<i<<"\t \r";
		//mexPrintf("i: %d\r", i);

		double h = 0;
		for(unsigned int j = 0; j < mesh.vertex(i).n_verts(); j ++){
    		unsigned int k = mesh.vertex(i).vert(j);
    		h +=  sqrt( fabs(dot(mesh.vertex(i).coord() - mesh.vertex(k).coord(),
         		      			 mesh.vertex(i).coord() - mesh.vertex(k).coord())) );
		}

		if(mesh.vertex(i).n_verts() > 0){
			h = hs * h / mesh.vertex(i).n_verts();
		}

		mexPrintf("i: %d h: %f\r", i, h);

		double hh = h * h;

		vgdists.clear();		
		compute_one2part_Geodesic_vdist(i, geod_mesh, algorithm, vgdists, h * rho);
		//compute_one2part_Euclidean_vdist(i, mesh, vgdists, h * rho);
		
		totalweight = 0;
		for(unsigned int j = 0; j < vgdists.size(); j ++){
			unsigned int vid = vgdists[j].first;
			if( vid == i ){
				continue;
			}

			double weight = exp(-vgdists[j].second * vgdists[j].second / hh) * ( 4.0 / (M_PI * hh * hh) );
			//cout<<"vid: "<<vid<<" dist: "<<vgdists[j].second<<" weight: "<<weight<<" vareas: "<<vareas[vid]<<endl;

			weight *=  vareas[vid];

			IIV.push_back(i + 1);
			JJV.push_back(vid + 1);
			SSV.push_back(weight);

			totalweight -= weight;
		}
		
		IIV.push_back(i + 1);
		JJV.push_back(i + 1);
		SSV.push_back(totalweight);
	}
	
	mexPrintf("\n");
}


void generate_Xu_Meyer_laplace_matrix(TMesh& mesh,  vector<unsigned int>& IIV,  vector<unsigned int>& JJV,  vector<double>& SSV, vector<double>& DDV)
{
	unsigned int nv = mesh.v_count();
	unsigned int nf = mesh.f_count();
		
	//compute the center and area for each facets
	vector<double> fareas;
	fareas.resize(3 * nf);
	vector<double> fcots;
	fcots.resize(3 * nf);

	unsigned int vid0, vid1, vid2;
	for(unsigned int f = 0; f < nf; f ++){
		vid0 = mesh.facet(f).vert(0);
		vid1 = mesh.facet(f).vert(1);
		vid2 = mesh.facet(f).vert(2);
		VECTOR3 vv = cross( mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord(), 
								  mesh.vertex(vid2).coord() - mesh.vertex(vid0).coord() );
		double area = norm(vv);
		
		double l0 = sqrt( fabs(dot(mesh.vertex(vid1).coord() - mesh.vertex(vid2).coord(), mesh.vertex(vid1).coord() - mesh.vertex(vid2).coord())) );
		double l1 = sqrt( fabs(dot(mesh.vertex(vid0).coord() - mesh.vertex(vid2).coord(), mesh.vertex(vid0).coord() - mesh.vertex(vid2).coord())) );
		double l2 = sqrt( fabs(dot(mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord(), mesh.vertex(vid1).coord() - mesh.vertex(vid0).coord())) );
		
		assert(l0 > MYNZERO && l1 > MYNZERO && l2 > MYNZERO);

		double h0 = area / l0;
		double h1 = area / l1;
		double h2 = area / l2;

		fcots[3 * f] = sqrt(fabs(l1 * l1 - h2 * h2)) / h2;
		if( (l1 * l1 + l2 * l2 - l0 * l0) < 0){
			fcots[3 * f] = -fcots[3 * f];
		}
		fcots[3 * f + 1] = sqrt(fabs(l2 * l2 - h0 * h0)) / h0;
		if( (l0 * l0 + l2 * l2 - l1 * l1) < 0){
			fcots[3 * f + 1] = -fcots[3 * f + 1];
		}
		fcots[3 * f + 2] = sqrt(fabs(l0 * l0 - h1 * h1)) / h1;
		if( (l0 * l0 + l1 * l1 - l2 * l2) < 0 ){
			fcots[3 * f + 2] = -fcots[3 * f + 2];
		}

		//--------------------------------------------------
		// if(mode == XU){
		// 	fareas[3 * f] = 1.0 / 8.0 * (l1 * l1 * fcots[3 * f + 1] + l2 * l2 * fcots[3 * f + 2]);
		// 	fareas[3 * f + 1] = 1.0 / 8.0 * (l0 * l0 * fcots[3 * f] + l2 * l2 * fcots[3 * f + 2]);
		// 	fareas[3 * f + 2] = 1.0 / 8.0 * (l1 * l1 * fcots[3 * f + 1] + l0 * l0 * fcots[3 * f]);
		// }	
		// else
		//-------------------------------------------------- 
		{
			if( fcots[3 * f] >= 0 && fcots[3 * f + 1] >= 0 && fcots[3 * f + 2] >= 0 ){
				fareas[3 * f] = 1.0 / 8.0 * (l1 * l1 * fcots[3 * f + 1] + l2 * l2 * fcots[3 * f + 2]);
				fareas[3 * f + 1] = 1.0 / 8.0 * (l0 * l0 * fcots[3 * f] + l2 * l2 * fcots[3 * f + 2]);
				fareas[3 * f + 2] = 1.0 / 8.0 * (l1 * l1 * fcots[3 * f + 1] + l0 * l0 * fcots[3 * f]);
			}
			else if(fcots[3 * f] < 0){
				fareas[3 * f] = area / 4;
				fareas[3 * f + 1] = area / 8;
				fareas[3 * f + 2] = area / 8;
			}
			else if(fcots[3 * f + 1] < 0){
				fareas[3 * f] = area / 8;
				fareas[3 * f + 1] = area / 4;
				fareas[3 * f + 2] = area / 8;
			}
			else{
				fareas[3 * f] = area / 8;
				fareas[3 * f + 1] = area / 8;
				fareas[3 * f + 2] = area / 4;
			}
		}
	}

	DDV.clear();
	DDV.resize(nv, 0);

	double Amix, totalweight;
	map<unsigned int, double> adj_weight;
	for(unsigned int i = 0; i < nv; i ++){
		//cout<<"i: "<<i<<"\t \r";
		mexPrintf("i: %d\r", i);


		Amix = 0;
		totalweight = 0;
		adj_weight.clear();
		for(unsigned int kk = 0; kk < mesh.vertex(i).n_verts(); kk ++){
			adj_weight[mesh.vertex(i).vert(kk)] = 0;
		}
		for(unsigned int kk = 0; kk < mesh.vertex(i).n_facets(); kk ++){
			unsigned int f = mesh.vertex(i).facet(kk);
			int ind = mesh.facet(f).index(i);
			assert( ind >=0 );

			Amix += fareas[3 * f + ind];

			vid1 = mesh.facet(f).vert((ind + 1) % 3);
			vid2 = mesh.facet(f).vert((ind + 2) % 3);
			//cout<<"vid1: "<<vid1<<" vid2: "<<vid2<<endl;
			map<unsigned int, double>::iterator iter;

			iter = adj_weight.find(vid1);
			assert(iter != adj_weight.end());
			iter->second = iter->second + fcots[3 * f + (ind + 2) % 3];

			iter = adj_weight.find(vid2);
			assert(iter != adj_weight.end());
			iter->second = iter->second + fcots[3 * f + (ind + 1) % 3];


			totalweight -= (fcots[3 * f + (ind + 2) % 3] + fcots[3 * f + (ind + 1) % 3]);

			//--------------------------------------------------
			// matrix[vid1 * nv + i] -= fcots[3 * f + (ind + 2) % 3];
			// matrix[vid2 * nv + i] -= fcots[3 * f + (ind + 1) % 3];
			// matrix[i * nv + i] += fcots[3 * f + (ind + 2) % 3] + fcots[3 * f + (ind + 1) % 3];
			//-------------------------------------------------- 
		}
		
		//cout<<"Amix: "<<Amix<<endl;
		assert(Amix > 0);
		Amix *= 2;
		
		for(map<unsigned int, double>::iterator iter = adj_weight.begin(); iter != adj_weight.end(); iter ++){
			IIV.push_back(i + 1);
			JJV.push_back(iter->first + 1);
			SSV.push_back(iter->second);
		}
		IIV.push_back(i + 1);
		JJV.push_back(i + 1);
		SSV.push_back(totalweight);

		DDV[i] = Amix;
	}	
	mexPrintf("\n");

}














void compute_one2part_Euclidean_vdist(unsigned int vid_start, TMesh& mesh, vector<pair<unsigned int, double> >& vgdists, double maxdist)
{
	for(unsigned int j = 0; j < mesh.v_count(); j ++){
		double d = sqrt( dot(mesh.vertex(vid_start).coord() - mesh.vertex(j).coord(),  
									mesh.vertex(vid_start).coord() - mesh.vertex(j).coord()) );
		if(d <= maxdist){
		 	vgdists.push_back( make_pair(j, d) );
		}
	}
	

	//--------------------------------------------------
	// unsigned int nv = mesh.v_count();
	// queue<unsigned int> vert_queue;
	// //mark a vertex which indicates if it is visited 
	// vector<bool> visited_vector;
	// visited_vector.resize(nv, false);
	// vert_queue.push(vid_start);
	// visited_vector[vid_start] = true;
	// //terminate the while loop if exhaust the vert_heap or compute 
	// //the geodesic distance to all other vertices which are connected 
	// //to the vertex vid
	// while( (!vert_queue.empty()) ){
	// 	//exact the first element from the heap  
	// 	unsigned int vid = vert_queue.front();
	// 	vert_queue.pop();
	// 	double d =  sqrt( dot(mesh.vertex(vid_start).coord() - mesh.vertex(vid).coord(), 
	// 				 mesh.vertex(vid_start).coord() - mesh.vertex(vid).coord()) );
	// 	if(d > maxdist){
	// 		continue;
	// 	}
	// 	vgdists.push_back(make_pair(vid, d));
	// 	//for each adjacent vertex on mesh, modify its geodesic distance and move its position 
	// 	//in the heap if necessary. 
	// 	for(unsigned int k = 0; k < mesh.vertex(vid).n_verts(); k ++){
	// 		unsigned int kk = mesh.vertex(vid).vert(k);
	// 		if( !visited_vector[kk]){
	// 			vert_queue.push(kk);
	// 			visited_vector[kk] = true;
	// 		}
	// 	}//for k
	// }//while
	// //sort(vgdists.begin(), vgdists.end());
	//-------------------------------------------------- 
}


void compute_one2part_Geodesic_vdist(unsigned int vid_start, geodesic::Mesh& geod_mesh, geodesic::GeodesicAlgorithmExact& algorithm, vector<pair<unsigned int, double> >& vgdists, double maxdist)
{
	geodesic::SurfacePoint source(&geod_mesh.vertices()[vid_start]);      //create source 
	vector<geodesic::SurfacePoint> all_sources(1,source);  //in general, there could be multiple sources, but now we have only one
	algorithm.propagate(all_sources, maxdist);   //cover the whole mesh

	double distance;
	for(unsigned int i = 0; i < geod_mesh.vertices().size(); ++ i){
   	geodesic::SurfacePoint p(&geod_mesh.vertices()[i]);
		unsigned int best_source = algorithm.best_source(p,distance);      //for a given surface point, find closets source and distance to this source
		if(distance <=  maxdist){
			vgdists.push_back( make_pair(i, distance) );
			//cout << distance <<endl;    //print geodesic distance for every vertex
		}
	}
}

