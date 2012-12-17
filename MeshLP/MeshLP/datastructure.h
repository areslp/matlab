#ifndef __DATASTRUCTURE_H__
#define __DATASTRUCTURE_H__


#include <vector>
#include <map>
#include <set>
#include <fstream>
#include <float.h>
#include <assert.h>
#include "matrix.h"

#define MYNZERO 1e-10


//====================================================================
//Begin: Base elements for curve and mesh (manifold and non-manifold)
//====================================================================

//-------------------------------------------------------------------
//VertexBase
//-------------------
class VertexBase{
public:
  //Constructor
  VertexBase(): _flag(0){}
  VertexBase(double x, double y, double z): _flag(0)
  {  _co[0] = x; _co[1] = y; _co[2] = z; }
  VertexBase(VECTOR3 co): _flag(0)
  {  _co = co; }
  
  VertexBase(const VertexBase& vb){
    _co = vb.coord(); 
    _flag = vb.flag();
  }
  //Overloaded operator
  VertexBase& operator = (const VertexBase& vb){
    _co = vb.coord(); 
    _flag = vb.flag();
    return *this;
  }
  //Function for accessing data
  VECTOR3 coord() const { return _co; }
  void set_coord(VECTOR3 co) { _co = co; }

  bool check_flag(unsigned char f) { return _flag&f; }
  void set_flag(unsigned char f) { _flag |= f; }
  void un_set_flag(unsigned char f) { _flag &= (~f); }
  unsigned char flag() const { return _flag; }
  void copy_flag(unsigned char f) { _flag = f; }
  
private:
  VECTOR3 _co;  
  unsigned char _flag;
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
//EdgeBase
//-------------------
class EdgeBase{
public:
  //Constructor
  EdgeBase(): _flag(0) {}
  EdgeBase(unsigned  int vs, unsigned int ve): _flag(0) { set_vs_ve(vs, ve); }
  EdgeBase(const EdgeBase& eb){
    _vs = eb.vs(); _ve = eb.ve();
    _flag = eb.flag();
  }
  //Overloaded operator
  EdgeBase& operator = (const EdgeBase& eb){
    _vs = eb.vs(); _ve = eb.ve();
    _flag = eb.flag();
    return *this;
  }
  
  //Functions for accessing data
  unsigned int vs() const { return _vs; }
  unsigned int ve() const { return _ve; }
  void set_vs_ve(unsigned  int vs, unsigned int ve) { 
    if(vs < ve) { _vs = vs; _ve = ve; }
    else { _vs = ve; _ve = vs; }
  }
  

  bool check_flag(unsigned char f) { return _flag&f; }
  void set_flag(unsigned char f) { _flag |= f; }
  void un_set_flag(unsigned char f) { _flag &= (~f); }
  unsigned char flag() const { return _flag; }
  void copy_flag(unsigned char f) { _flag = f; }

  //Functions for some simple operations
  inline int index(unsigned int vid){
    if(_vs == vid) return 0;
    else if(_ve == vid) return 1;
    else return -1;
  }
  unsigned int opposite_vertex(unsigned int vid) {
    if( vid == _vs ) return _ve;
    else {
      assert(vid == _ve);
      return _vs;
    }
  }

private:   
  unsigned int _vs, _ve;
  unsigned char _flag;
};
//-------------------------------------------------------------------

//-------------------------------------------------------------------
//TFacetBase
//-------------------
//Base of triangle facet 
//-------------------------------------------------------------------
class TFacetBase{
public:
  //Constructors
  TFacetBase():_flag(0){
    _facets[0] = -1; _facets[1] = -1; _facets[2] = -1;
  }
  TFacetBase(unsigned int v0, unsigned int v1, unsigned int v2):_flag(0){
    _verts[0] = v0; _verts[1] = v1; _verts[2] = v2; 
    _facets[0] = -1; _facets[1] = -1; _facets[2] = -1;
  }
  TFacetBase(const TFacetBase& fb){
    _facets[0] = fb.facet(0); _facets[1] = fb.facet(1); _facets[2] = fb.facet(2);
    _verts[0] = fb.vert(0); _verts[1] = fb.vert(1); _verts[2] = fb.vert(2); 
    _flag = fb.flag();
  }
  //Assignment operators
  TFacetBase& operator = (const TFacetBase& fb){
    _facets[0] = fb.facet(0); _facets[1] = fb.facet(1); _facets[2] = fb.facet(2);
    _verts[0] = fb.vert(0); _verts[1] = fb.vert(1); _verts[2] = fb.vert(2); 
    _flag = fb.flag();
    return *this;
  }
  
  //Functions for accessing data
  int facet(int i) const {return _facets[i];}
  void set_facet(int i, int f){ _facets[i] = f;}
  unsigned int vert(int i) const {return _verts[i];}
  void set_vert(int i, unsigned int v){ _verts[i] = v;}
  
  bool check_flag(unsigned short f) { return _flag&f; }
  void set_flag(unsigned short f) { _flag |= f; }
  void un_set_flag(unsigned short f) { _flag &= (~f); }
  unsigned short flag() const { return _flag; }
  void copy_flag(unsigned short f) { _flag = f; }

  //Functions for some simple operations
  int index(unsigned int v){ 
    if(_verts[0] == v) return 0;
    else if(_verts[1] == v) return 1;
    else if(_verts[2] == v) return 2;
    else return -1; 
  }
  
  int findex(int f){ 
    if(_facets[0] == f) return 0;
    else if(_facets[1] == f) return 1;
    else if(_facets[2] == f) return 2;
    else return -1; 
  }
  
  //using this function instead of the above for "bad"(such as non manifold) mesh. 
  int findex(unsigned int vid0, unsigned int vid1){
    if( index(vid0) < 0 || index(vid1) < 0 ) return -1;
    if(_verts[0] != vid0 && _verts[0] != vid1) return 0;
    if(_verts[1] != vid0 && _verts[1] != vid1) return 1;
    if(_verts[2] != vid0 && _verts[2] != vid1) return 2;
    
    cerr<<"Warning: invalid facet or parameter of function findex"<<endl;
    return -1;
  }

private:
  unsigned int _verts[3];
  int _facets[3];
  unsigned short _flag;
};

//-------------------------------------------------------------------
//PFaceBase: 
//-------------------
//Base of polygon facet 
//-------------------------------------------------------------------
class PFacetBase{
public:
  //Constructors
  PFacetBase():_flag(0){}
  PFacetBase(const PFacetBase& fb){
    copy_flag( fb.flag() );
    for(unsigned int i = 0; i < fb.n_verts(); i ++)
      add_vert( fb.vert(i) );
    for(unsigned int i = 0; i < fb.n_facets(); i ++)
      add_facet( fb.facet(i) );
  }
  //Overloaded operators
  PFacetBase& operator = (const PFacetBase& fb){
    copy_flag( fb.flag() );
    for(unsigned int i = 0; i < fb.n_verts(); i ++)
      add_vert( fb.vert(i) );
    for(unsigned int i = 0; i < fb.n_facets(); i ++)
      add_facet( fb.facet(i) );
    return *this;
  }
  //
  void clear(){ _facets.clear(); _verts.clear(); }
  
  //Functions for accessing data
  bool check_flag(unsigned short f) { return _flag&f; }
  void set_flag(unsigned short f) { _flag |= f; }
  void un_set_flag(unsigned short f) { _flag &= (~f); }
  unsigned short flag() const { return _flag; }  
  void copy_flag(unsigned short f) { _flag = f; }
  
               
  unsigned int n_verts() const {return _verts.size();}  
  unsigned int vert(unsigned int i) const {return _verts[i];}
  void set_vert(unsigned int i, unsigned int vid){ _verts[i] = vid; }
  unsigned int add_vert(unsigned int v){ _verts.push_back(v); return _verts.size() - 1;}
  
  unsigned int n_facets() const {return _facets.size();}
  unsigned int facet(unsigned int eind) const {return _facets[eind];}
  void set_facet(unsigned int eind, unsigned int fid) { _facets[eind] = fid; }
  unsigned int add_facet(unsigned int fid){ _facets.push_back(fid); return _facets.size() - 1;}
  
  //Functions for some simple operations
  int index(unsigned int v){ 
    for(unsigned int  i = 0; i < _verts.size(); i ++)
      if( v == _verts[i])  return i;
    return -1;
  }
  int findex(unsigned int f){
    for(unsigned int i = 0; i < n_facets(); i ++){
      if( f == _facets[i] ) return i;    
    }
    return -1;  
  }

private:
  vector<unsigned int> _verts;
  vector<unsigned int> _facets;    //_facets[i] contains the neighboring polygons about edge _verts[i]--_verts[i + 1]
  unsigned short _flag;
};

//-------------------------------------------------------------------
//NPFaceBase: 
//-------------------
//Base of polygon facet for non-manifold 
class NPFacetBase{
public:
  //Constuctors
  NPFacetBase():_flag(0){}
  NPFacetBase(const NPFacetBase& fb){
    copy_flag( fb.flag() );
    for(unsigned int i = 0; i < fb.n_verts(); i ++)
      add_vert( fb.vert(i) );
    _facetss.resize( fb.n_facetss() );
    for(unsigned int i = 0; i < fb.n_facetss(); i ++){
      for(unsigned int j = 0; j < fb.n_ith_facets(i); j ++){
        add_facet( i, fb.facet(i, j) );
      }//for j
    }//for i  
  }
  //Overloaded operators
  NPFacetBase& operator = (const NPFacetBase& fb){
    copy_flag( fb.flag() );
    for(unsigned int i = 0; i < fb.n_verts(); i ++)
      add_vert( fb.vert(i) );
    _facetss.resize( fb.n_facetss() );
    for(unsigned int i = 0; i < fb.n_facetss(); i ++){
      for(unsigned int j = 0; j < fb.n_ith_facets(i); j ++){
        add_facet( i, fb.facet(i, j) );
      }//for j
    }//for i  
    return *this;
  }
  //
  void clear(){ _facetss.clear(); _verts.clear(); }
  
  //Functions for accessing data
  bool check_flag(unsigned short f) { return _flag&f; }
  void set_flag(unsigned short f) { _flag |= f; }
  void un_set_flag(unsigned short f) { _flag &= (~f); }
  unsigned short flag() const { return _flag; }  
  void copy_flag(unsigned short f) { _flag = f; }
  
               
  unsigned int n_verts() const {return _verts.size();}  
  unsigned int vert(unsigned int i) const {return _verts[i];}
  void set_vert(unsigned int i, unsigned int vid){ _verts[i] = vid; }
  unsigned int add_vert(unsigned int v){ _verts.push_back(v); return _verts.size() - 1;}
  
  unsigned int n_facetss() const {return _facetss.size();}
  unsigned int n_ith_facets(unsigned int i) const { return _facetss[i].size(); }
  void resize_facetss(unsigned int n) { _facetss.resize(n); } 
  
  unsigned int facet(unsigned int eind, unsigned int find) const {return (_facetss[eind])[find];}
  void set_facet(unsigned int eind, unsigned int find, unsigned int fid) { (_facetss[eind])[find] = fid; }
  unsigned int add_facet(unsigned int eind, unsigned int fid){ _facetss[eind].push_back(fid); return _facetss[eind].size() - 1;}

  //Functions for some simple operations
  int index(unsigned int v){ 
    for(unsigned int  i = 0; i < _verts.size(); i ++)
      if( v == _verts[i])  return i;
    return -1;
  }
  
  int findex(unsigned int f){
    for(unsigned int i = 0; i < n_facetss(); i ++){
      for(unsigned int j = 0; j < n_ith_facets(i); j ++)
        if( f == facet(i, j) ) return i;    
    }
    return -1;  
  }
  bool findex(unsigned int f, unsigned int& eind, unsigned int& find){ 
    for(unsigned int i = 0; i < n_facetss(); i ++){
      for(unsigned int j = 0; j < n_ith_facets(i); j ++){
        if( f == facet(i, j) ){
          eind = i; find = j; return true;
        }    
      }//for j
    }//for i
    return false;  
  }

  void remove_facet(unsigned int eind, unsigned int fid){
    for(unsigned int i = 0; i < n_ith_facets(eind); i ++){
      if(fid == facet(eind, i)){
        (_facetss[eind])[i] = facet(eind, n_ith_facets(eind) - 1);
        _facetss[eind].pop_back();
      }    
    }//for i  
  }  

private:
  vector<unsigned int> _verts;
  vector< vector<unsigned int> > _facetss;    //_facetss[i] contains the neighboring polygons about edge _verts[i]--_verts[i + 1]
  unsigned short _flag;
};
//-------------------------------------------------------------------

//====================================================================
//End: Base elements for curve and mesh (manifold and non-manifold)
//====================================================================

enum pmode{
	EIGEN, LLF, MEANC
};
enum mmode{
	MEYER, XU, OUR	
};

enum nmode{
	DENSE, SPARSE	
};

enum stype{
	PLANE, SPHERE, SKINSURF
};


#endif //__DATASTRUCTURE_H__
