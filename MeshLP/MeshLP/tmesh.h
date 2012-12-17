#ifndef __TMESH_H__
#define __TMESH_H__

#include "datastructure.h"

//-------------------------------------------------------------------
//VTMesh: vertex class for triangle mesh
//-------------------
#define VTMESH_FLAG_BOUNDARY  0x1
#define VTMESH_FLAG_NMANIFOLD 0x2
#define VTMESH_FLAG_SELECTED 	0x4

class VTMesh : public VertexBase{
public:
  //Constructor
  VTMesh():VertexBase(){}
  VTMesh(double x, double y, double z):VertexBase(x, y, z){}
  VTMesh(VECTOR3 co):VertexBase(co){}
  VTMesh(const VTMesh& v):VertexBase(v){
    for(unsigned int i = 0; i < v.n_verts(); i ++)
      add_vert( v.vert(i) );
    for(unsigned int i = 0; i < v.n_facets(); i ++)
      add_facet( v.facet(i) );
  }

  //Assignment operator
  VTMesh& operator = (const VTMesh& v){
    if (this != &v){
      this->VTMesh::~VTMesh();
      new (this) VTMesh(v);
    }
    return *this;
  }
                          
  //Functions for accessing data
  unsigned int n_verts() const { return _verts.size();  }
  unsigned int vert(unsigned int i) const { return _verts[i]; }
  unsigned int add_vert(unsigned int vid){ _verts.push_back(vid); return _verts.size() - 1; }
  void remove_vert(unsigned int vid){
    for(unsigned int i = 0; i < _verts.size(); i ++){
      if(vid == _verts[i]){
        _verts[i] = _verts[_verts.size() - 1];
        _verts.pop_back();
        break;
      }
    }//for i
  }
  unsigned int add_unique_vert(unsigned int vid){
    for(unsigned int i = 0; i < _verts.size(); i ++)
      if(vid == _verts[i]) return i;
    _verts.push_back(vid);
    return _verts.size() - 1;
  }
  unsigned int n_facets() const { return _facets.size();  }
  unsigned int facet(unsigned int i)  const { return _facets[i];  }
  unsigned int add_facet(unsigned int fid){ _facets.push_back(fid); return _facets.size() - 1; }
  void remove_facet(unsigned int fid){
    for(unsigned int i = 0; i < _facets.size(); i ++){
      if( fid == _facets[i] ){
        _facets[i] = _facets[_facets.size() - 1];
        _facets.pop_back();
        break;
      }
    }//for i
  }                                                  
  unsigned int add_unique_facet(unsigned int fid){
    for(unsigned int i = 0; i < _facets.size(); i ++)
      if(fid == _facets[i])  return i;;
    _facets.push_back(fid);
    return _facets.size() - 1; 
  }
  
  //Functions for some simple operations
  int index(unsigned int vid){
    for(unsigned int i = 0; i < _verts.size(); i ++)
      if(vid == _verts[i]) return i;
    return -1;
  }

private:
  vector<unsigned int> _facets;
  vector<unsigned int> _verts;
                  
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
//FTMesh: facet class for triangle mesh
//-------------------
#define FTMESH_FLAG_NMANIFOLD 0x01
#define FTMESH_FLAG_REVERSE   0x02
#define FTMESH_FLAG_ORIENTATE 0x10

#define FTMESH_FLAG_BOUNDARY  0x20

#define FTMESH_FLAG_VISITED   0x04
#define FTMESH_FLAG_IN        0x08
#define FTMESH_FLAG_SELECT    0x40
class FTMesh : public TFacetBase{
public:
  //Constructor
  FTMesh():TFacetBase(){}
  FTMesh(unsigned int v0, unsigned int v1, unsigned int v2):TFacetBase(v0, v1, v2){}
  FTMesh(const FTMesh& f):TFacetBase(f){}
  //Assignment operator
  FTMesh& operator = (const FTMesh& f){
    if(this != &f){
      this->FTMesh::~FTMesh();
      new (this) FTMesh(f);
    }
    return *this;
  }
private:
};
//-------------------------------------------------------------------

class TMesh{
public:
  //Constructor
  TMesh(){};
  //Destructor
  virtual ~TMesh(){};

  //Functions for Creation
  bool ReadOffFile(char *filename, bool wcolor);
  bool ReadOffFile(char *filename);
  bool ReadObjFile(char *filename);
  void GenerateMeshTopo();

  bool CheckMeshTopo();
  void PrintMeshTopo();
  void OrientateFacets();
  void MarkNonManifoldness();

  //Functions for output 
  void OutMeshOffFile(char *filename);
  void OutMeshOffFile(FILE *fp, double r = 0, double g = 1, double b = 0, double a = 0.2);
                          
  
  //Function for accessing data
  unsigned int v_count() const { return _vertices.size(); }
  unsigned int f_count() const { return _facets.size(); }

  unsigned int add_vertex(const VTMesh& v){ _vertices.push_back(v); return (_vertices.size() - 1);}
  VTMesh& vertex(unsigned int i ) { return _vertices[i]; }
  unsigned int add_facet(const FTMesh& f){ _facets.push_back(f); return (_facets.size() - 1);}
  FTMesh& facet(unsigned int i) { return _facets[i]; }
		
  void copy_name(char* name) { sprintf(name, "%s", _name); }	
  void set_name(char* name){ sprintf(_name, "%s", name );}
 
  //Functions for bounding box
  void GetBBox();
  bool IsOutBBox(VECTOR3 pt);                 
  void MeshSize(double& maxs, double& mins, double& aves);
  VECTOR3 pmin(){ return _pmin; }
  VECTOR3 pmax(){ return _pmax; }

  
  //Clear
  void clear(){ _vertices.clear();  _facets.clear(); }
            
private:
  char _name[256];
  vector<VTMesh> _vertices;
  vector<FTMesh> _facets;
  VECTOR3 _pmin;
  VECTOR3 _pmax;
};

#endif //_TMESH_H__

