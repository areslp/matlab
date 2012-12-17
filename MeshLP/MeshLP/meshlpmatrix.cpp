#include <mex.h>
#include "comp_meshlpmatrix.h"
#include "meshlpmatrix.h"
#include <string>

#include "engine.h"

bool generate_sym_meshlp_matrix_matlab(char* filename, unsigned int htype, double hs, double rho, double& h, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV, vector<double>& AAV)
{
	TMesh tmesh;
	if( !(tmesh.ReadOffFile(filename)) ){
   	mexPrintf("Failed to open file %s\n", filename);
  	}

   double maxs, mins, avers;                                                                         
	tmesh.MeshSize(maxs, mins, avers);
	mexPrintf("avers: %f\n", avers);
	if(htype == 0){
		h = avers * hs;
	}
	else{
		h = hs;
	}

	generate_sym_meshlp_matrix(tmesh, h, rho, IIV,  JJV,  SSV, AAV);
	mexPrintf("h: %f\n", h);
}

bool generate_sym_meshlp_matrix_geod_matlab(char* filename, unsigned int htype, double hs, double rho, double& h, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV, vector<double>& AAV)
{
	TMesh tmesh;
	if( !(tmesh.ReadOffFile(filename)) ){
   	mexPrintf("Failed to open file %s\n", filename);
  	}

   double maxs, mins, avers;                                                                         
	tmesh.MeshSize(maxs, mins, avers);
	mexPrintf("avers: %f\n", avers);
	if(htype == 0){
		h = avers * hs;
	}
	else{
		h = hs;
	}

	generate_sym_meshlp_matrix_geod(tmesh, h, rho, IIV,  JJV,  SSV, AAV);
	mexPrintf("h: %f\n", h);
}

//bool generate_sym_meshlp_matrix_geod_matlab_fastmarching(char* filename, unsigned int htype, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV, vector<double>& AAV)
//{
	//TMesh tmesh;
	//if( !(tmesh.ReadOffFile(filename)) ){
      //mexPrintf("Failed to open file %s\n", filename);
     //}
	//double h;
   //double maxs, mins, avers;                                                                         
	//tmesh.MeshSize(maxs, mins, avers);
	//mexPrintf("avers: %f\n", avers);
	//if(htype == 0){
		//h = avers * hs;
	//}
	//else{
		//h = hs;
	//}
	//Engine* engine;
   //if (!(engine = engOpen("\0"))) {
		//cerr<<"Can't start MATLAB engine\n";
		//exit(1);
	//}
	//engEvalString(engine, "path(path, './model')");
	//engEvalString(engine, "path(path, './toolbox_fast_marching')");
	//engEvalString(engine, "path(path, './toolbox_fast_marching/toolbox')");
	//char str[1024];
	//string str_path(filename);
	//string str_file = str_path.substr( str_path.find_last_of( '/' ) +1 );
	////cout<<str_file.c_str()<<endl;
	//sprintf(str, "[vertex,faces] = read_mesh('%s');", str_file.c_str());
	////sprintf(str, "[vertex,faces] = read_mesh('sphere_uniform_1000.off');");

	//engEvalString(engine, str);
	//engEvalString(engine, "options.name = 'dump'");
	//engEvalString(engine, "options.end_points = []");
	//mexPrintf("h: %f\n", h);
	//generate_sym_meshlp_matrix_geod_fastmarching(tmesh, engine, h, rho, IIV,  JJV,  SSV, AAV);
	//engClose(engine);
//}

bool generate_meshlp_matrix_matlab(char* filename, unsigned int htype, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV)
{
	TMesh tmesh;
	if( !(tmesh.ReadOffFile(filename)) ){
   	mexPrintf("Failed to open file %s\n", filename);
  	}

	double h;
   double maxs, mins, avers;                                                                         
	tmesh.MeshSize(maxs, mins, avers);
	mexPrintf("avers: %f\n", avers);
	if(htype == 0){
		h = avers * hs;
	}
	else{
		h = hs;
	}

	mexPrintf("h: %f\n", h);
	generate_meshlp_matrix(tmesh, h, rho, IIV,  JJV,  SSV);
}

bool generate_meshlp_matrix_geod_matlab(char* filename, unsigned int htype, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV)
{
	TMesh tmesh;
	if( !(tmesh.ReadOffFile(filename)) ){
   	mexPrintf("Failed to open file %s\n", filename);
  	}

	double h;
   double maxs, mins, avers;                                                                         
	tmesh.MeshSize(maxs, mins, avers);
	mexPrintf("avers: %f\n", avers);
	if(htype == 0){
		h = avers * hs;
	}
	else{
		h = hs;
	}

	mexPrintf("h: %f\n", h);
	generate_meshlp_matrix_geod(tmesh, h, rho, IIV,  JJV,  SSV);
}

bool generate_meshlp_matrix_adp_matlab(char* filename, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV)
{
	TMesh tmesh;
	if( !(tmesh.ReadOffFile(filename)) ){
   	mexPrintf("Failed to open file %s\n", filename);
  	}
	generate_meshlp_matrix_adp(tmesh, hs, rho, IIV,  JJV,  SSV);
}

bool generate_meshlp_matrix_adp_geod_matlab(char* filename, double hs, double rho, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV)
{
	TMesh tmesh;
	if( !(tmesh.ReadOffFile(filename)) ){
   	mexPrintf("Failed to open file %s\n", filename);
  	}
	generate_meshlp_matrix_adp_geod(tmesh, hs, rho, IIV,  JJV,  SSV);
}

void generate_Xu_Meyer_laplace_matrix_matlab(char* filename, vector<unsigned int>& IIV, vector<unsigned int>& JJV, vector<double>& SSV, vector<double>& DDV)
{
	TMesh tmesh;
	if( !(tmesh.ReadOffFile(filename)) ){
   	mexPrintf("Failed to open file %s\n", filename);
  	}
	generate_Xu_Meyer_laplace_matrix(tmesh,  IIV,  JJV,  SSV, DDV);
}


