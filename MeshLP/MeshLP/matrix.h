///////////////////////////////////////////////////////////////////
//   Purpose:    vector and matrix class and operations
//   Paper:      Hierarchical Polygon Tiling with Coverage Masks
//   Created by: Matt Camuto
//
//   Modification:
//        Time            Programmer
//        07-09-99        R. Wenger
//        08-20-99        K. Boner
//        07-15-00        J. Gao
//		  08-09-03		  J. Sun
///////////////////////////////////////////////////////////////////


#ifndef MATRIX_H
#define MATRIX_H

///////////////////////////////////////////////////////////////////
// 2d, 3d and 4d vector and matrix classes.
// Supports vector and matrix addition and subtraction, 
// matrix multiplication, vector dot products, cross products, 
// normalization, etc.
///////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
using namespace std;

typedef unsigned int INT_32;

extern void init_shift();
extern int is_bit_set(INT_32 value, int pos);
extern void set_bit(INT_32 &value, int pos);


// 16 bit structure used for bit mask
class INT_64 {

private :
	INT_32 vec[2];

public :
    INT_64()                                    // constructor
      { Zero();};
    INT_64(const INT_32 x0, const INT_32 x1) // constructor
      { vec[0] = x0; vec[1] = x1;};
		~INT_64() {};
    int Dimension() const
      { return 2; };
    INT_32 & operator [](const int i)           // index i'th element
      { return(vec[i]); };
    INT_32 operator ()(const int i) const        // return i'th element
      { return(vec[i]); };
    const INT_64 & operator =(const INT_64 & v0)  // copy vector v0
      { vec[0] = v0(0); vec[1] = v0(1);
        return(*this); };
    void Zero()                                  // make zero vector
      { vec[0] = vec[1] = 0; };
	void Fullset() {
        vec[0] = ~0;
        vec[1] = ~0;
	};
	void combine(INT_64 &v0)
	{
		vec[0] = vec[0] | v0[0];
		vec[1] = vec[1] | v0[1];
	}
	void combine2(INT_32 &v0, INT_32 &v1)
	{
		vec[0] = vec[0] | v0;
		vec[1] = vec[1] | v1;
	}

	void set_bit(int pos);
	int is_set(int pos);
	int is_zero() { 
		if (vec[0] || vec[1])
			return 0; 
		return 1;
	};
	int is_full() {
        INT_32 full = ~0;
        if ((vec[0]==full)&&(vec[1]==full))
            return 1;        
        return 0;
	};
	void print(FILE *fp) const {fprintf(fp, "%d, %d", vec[0], vec[1]); };
};


// 2d vector used for Point
class VECTOR2 {

private :
    double vec[2];

public :
    VECTOR2()                                    // constructor
		{ Zero(); };
    VECTOR2(const double x0, const double x1) // constructor
		{ vec[0] = x0; vec[1] = x1;};
    int Dimension() const
		{ return 2; };
    double & operator [](const int i)           // index i'th element
		{ return(vec[i]); };
    double operator ()(const int i) const        // return i'th element
		{ return(vec[i]); };
    const VECTOR2 & operator =(const VECTOR2 & v0)  // copy vector v0
		{ vec[0] = v0(0); vec[1] = v0(1); 
        return(*this); };
			
	const VECTOR2  operator *(const double t)		//number mulitple overload
		{ return VECTOR2(vec[0] * t, vec[1] * t); }
	const double  operator *(const VECTOR2 & v0)		//dot product overload
		{ return vec[0] * v0.vec[0] + vec[1] * v0.vec[1]; }
	const VECTOR2  operator +(const VECTOR2 & v0)	//add overload
		{ return VECTOR2(vec[0] + v0.vec[0], vec[1] + v0.vec[1]); }
//		{ return VECTOR2(vec[0], vec[1]); }
	const VECTOR2  operator - (const VECTOR2 & v0)	//sub overload
		{ return VECTOR2(vec[0] - v0.vec[0], vec[1] - v0.vec[1]); }
//		{ return VECTOR2(vec[0], vec[1]); }

	void Normalize();

    void Zero()                                  // make zero vector
		{ vec[0] = vec[1] = 0.0; };
	void Set(const double x0, const double x1)
		{ vec[0] = x0; vec[1] = x1; };
};

// Line class
typedef struct line_2d {
	VECTOR2 start;
	VECTOR2 end;
}Line;

// 3d vector class
class VECTOR3 {

private :
public :
    double vec[3];

public :
		friend class VECTOR4;
    VECTOR3()                                    // constructor
		{	Zero(); };
	VECTOR3(const double x0, const double x1, const double x2) // constructor
		{	vec[0] = x0; vec[1] = x1; vec[2] = x2; };
	void Set(const double x0, const double x1, const double x2)
		{	vec[0] = x0; vec[1] = x1; vec[2] = x2; };
    int Dimension() const
		{	return 3; };
    double & operator [](const int i)           // index i'th element
		{	return(vec[i]); };
    double operator ()(const int i) const        // return i'th element
		{	return(vec[i]); };
    const VECTOR3 & operator =(const VECTOR3 & v0)  // copy vector v0
		{	vec[0] = v0(0); vec[1] = v0(1); vec[2] = v0(2); 
			return(*this); };
    void Zero()                                  // make zero vector
		{	vec[0] = vec[1] = vec[2] = 0.0; };
	double GetMax();
	void Clamp();
	const double Magnitude(){
		return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]); 
	}
	
	void Normalize()
	{	  
    double norm = Magnitude();
    if(norm != 0){
  		norm = 1.0 / norm;		
	  	vec[0] = vec[0] * norm; vec[1] = vec[1] * norm; vec[2] = vec[2] * norm; 
    }
	}

	const VECTOR3 Hat()
	{	double norm = Magnitude(); 
		if(norm != 0){
			norm = 1.0 / norm;
			return VECTOR3(vec[0] * norm, vec[1] * norm, vec[2] * norm);
		}
		return VECTOR3(vec[0], vec[1], vec[2]);
	}
/*	
	const VECTOR3 operator *(const double t) //number multiply 
	{	return VECTOR3( vec[0] * t, vec[1] * t, vec[2] * t );	}
	const double operator *(const VECTOR3 & v0) //dot product (v' * v)
	{	return vec[0] * v0(0) + vec[1] * v0(1) + vec[2] * v0(2);	}
	const VECTOR3 operator ^(const VECTOR3 & v0) //cross product 
	{	return VECTOR3(	vec[1] * v0(2) - vec[2] * v0(1),  
		vec[2] * v0(0) - vec[0] * v0(2),
		vec[0] * v0(1) - vec[1] * v0(0) ); }
	const VECTOR3 operator +(const VECTOR3 & v0) //add overload
	{	return VECTOR3( vec[0] + v0(0), vec[1] + v0(1), vec[2] + v0(2) );	}
	const VECTOR3 operator -(const VECTOR3 & v0) //sub overload
	{	return VECTOR3( vec[0] - v0(0), vec[1] - v0(1), vec[2] - v0(2) );	}
*/	
	
	const VECTOR3& operator +=(const VECTOR3 & v0) //+= overload
	{	vec[0] += v0(0); vec[1] += v0(1); vec[2] += v0(2);
		return (*this); }
    const VECTOR3& operator -=(const VECTOR3 & v0) //+= overload
	{   vec[0] -= v0(0); vec[1] -= v0(1); vec[2] -= v0(2);
		return (*this); }

};

// 4d vector class
class VECTOR4 {

private :
    double vec[4];

public :
    VECTOR4()                                    // constructor
		{	Zero(); };
    VECTOR4(VECTOR3 v)    // constructor
		{	vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; vec[3] = 1.0; };
	VECTOR4(VECTOR3 v, double x)
		{   vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; vec[3] = x; };
    VECTOR4(const double x0, const double x1, 
	    const double x2, const double x3)    // constructor
		{	vec[0] = x0; vec[1] = x1; vec[2] = x2; vec[3] = x3; };

    int Dimension() const
		{	return 4; };
    double & operator [](const int i)            // index i'th element
		{	return(vec[i]); };
    double operator ()(const int i) const        // return i'th element
		{	return(vec[i]); };
    const VECTOR4 & operator =(const VECTOR4 & v0)  // copy vector v0
		{	vec[0] = v0(0); vec[1] = v0(1); vec[2] = v0(2); vec[3] = v0(3);
			return(*this); };
    const VECTOR4 & operator =(const VECTOR3 & v0)  // copy vector v0
		{	vec[0] = v0(0); vec[1] = v0(1); vec[2] = v0(2); vec[3] = 1.0;
			return(*this); };
	
	VECTOR3 get_vector3()
		{	VECTOR3 temp(vec[0]/vec[3], vec[1]/vec[3], vec[2]/vec[3]); return temp;};
    void Zero()                                  // make zero vector
		{	vec[0] = vec[1] = vec[2] = vec[3] = 0.0; };
    void Normalize();                            // normalize vector
};

// 3d matrix
class MATRIX3 {

private :
    VECTOR3 mat[3];       // a vector represents each matrix row

public :

    MATRIX3()                                    // constructor
		{	Identity(); };
    MATRIX3(const VECTOR3 & v0, const VECTOR3 & v1, const VECTOR3 & v2)
		{	mat[0] = v0; mat[1] = v1; mat[2] = v2; };  // constructor
    int Dimension() const
		{	return 3; };
    VECTOR3 & operator [](const int i)           // index row i
		{	return(mat[i]); };
    // Note: reference row i, column j of MATRIX3 m0 as m0[i][j] (not m0[i,j])
    VECTOR3 operator()(const int i) const        // return row i
		{	return(mat[i]); };
    double operator ()(const int i, const int j) const   
		{	return(mat[i](j)); };                    // return element (i,j)
    MATRIX3 & operator =(const MATRIX3 & m0)     // copy matrix m0
		{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2); 
			return(*this); };
    void Identity();                             // set to identity
    void Zero();
    
/*   
   const MATRIX3 operator *(const double t) //number multiply 
    {	return MATRIX3( mat[0] * t, mat[1] * t, mat[2] * t );	}

   const MATRIX3 operator +(const MATRIX3 & m0) //add overload
    {	return MATRIX3( mat[0] + m0(0), mat[1] + m0(1), mat[2] + m0(2) );}
    const MATRIX3 operator -(const MATRIX3 & m0) //sub overload
    {	return MATRIX3( mat[0] - m0(0), mat[1] - m0(1), mat[2] - m0(2) );}
*/
    const void operator +=(const MATRIX3 & m0) //+= overload
	{	mat[0] += m0(0); mat[1] += m0(1); mat[2] += m0(2);}
    const void operator -=(const MATRIX3 & m0) //-= overload
    {	mat[0] -= m0(0); mat[1] -= m0(1); mat[2] -= m0(2);}
 
    
};

// 4d matrix
class MATRIX4 {

private :
    VECTOR4 mat[4];       // a vector represents each matrix row

public :		

    MATRIX4()                                    // constructor
		{	Identity(); };
    MATRIX4(const VECTOR4 & v0, const VECTOR4 & v1, 
			const VECTOR4 & v2, const VECTOR4 & v3) // constructor
		{	mat[0] = v0; mat[1] = v1; mat[2] = v2; mat[3] = v3; };  
    int Dimension() const
		{	return 4; };
    VECTOR4 & operator [](int i)                 // index row i
		{	return(mat[i]); };
    // Note: reference row i, column j of MATRIX4 m0 as m0[i][j] (not m0[i,j])
    VECTOR4 operator()(const int i) const        // return row i
		{	return(mat[i]); };
    double operator ()(const int i, const int j) const
		{	return(mat[i](j)); };                    // return element (i,j)
    MATRIX4 & operator =(const MATRIX4 & m0)     // copy matrix m0
		{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2); mat[3] = m0(3);
			return(*this); };
    MATRIX4 & operator =(const MATRIX3 & m0)     // copy matrix m0
		{	mat[0] = m0(0); mat[1] = m0(1); mat[2] = m0(2); 
			VECTOR4 temp(0.0,0.0,0.0,1.0);
			mat[3] = temp;
			return(*this); };
    void Identity();                             // set to identity
};


//************************
// INT_64 operations
//************************

inline INT_64 operator &(const INT_64 & v0, const INT_64 & v1)
// return v0 & v1
{	return(INT_64(v0(0)&v1(0), v0(1)&v1(1))); };

inline INT_64 operator |(const INT_64 & v0, const INT_64 & v1)
// return v0 | v1
{	return(INT_64(v0(0)|v1(0), v0(1)|v1(1))); };

inline INT_64 operator ~(const INT_64 & v0)
// return ~v0 
{	return(INT_64(~v0(0), ~v0(1))); };


//************************
// VECTOR2 operations
//************************
inline double dot(const VECTOR2 & v0, const VECTOR2 & v1)   
// return dot product of v0 and v1
{	return(v0(0)*v1(0) + v0(1)*v1(1)); };

inline VECTOR2 operator +(const VECTOR2 & v0, const VECTOR2 & v1)
// return v0 + v1
{	return(VECTOR2(v0(0) + v1(0), v0(1) + v1(1))); };

inline VECTOR2 operator -(const VECTOR2 & v0, const VECTOR2 & v1)
// return v0 - v1
{	return(VECTOR2(v0(0) - v1(0), v0(1) - v1(1))); };

inline VECTOR2 operator *(double x0, const VECTOR2 & v0)
// return x0*v0
{	return(VECTOR2(x0*v0(0), x0*v0(1))); };

inline VECTOR2 operator *(const VECTOR2 & v0, double x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); };
inline double norm(const VECTOR2& v)
{   return sqrt(v(0)*v(0) + v(1)*v(1)); }

inline double norm2(const VECTOR2& v)
{   return v(0)*v(0) + v(1)*v(1);}


//************************
// VECTOR3 operations
//************************

inline double dot(const VECTOR3 & v0, const VECTOR3 & v1)   
// return dot product of v0 and v1
{	return(v0(0)*v1(0) + v0(1)*v1(1) + v0(2)*v1(2)); };

inline VECTOR3 cross(const VECTOR3 & v0, const VECTOR3 & v1)
// return cross product of v0 and v1
{	return(VECTOR3(v0(1)*v1(2) - v0(2)*v1(1), 
			v0(2)*v1(0) - v0(0)*v1(2),
			v0(0)*v1(1) - v0(1)*v1(0))); };

inline VECTOR3 operator +(const VECTOR3 & v0, const VECTOR3 & v1)
// return v0 + v1
{	return(VECTOR3(v0(0) + v1(0), v0(1) + v1(1), v0(2) + v1(2))); };

inline VECTOR3 operator -(const VECTOR3 & v0, const VECTOR3 & v1)
// return v0 - v1
{	return(VECTOR3(v0(0) - v1(0), v0(1) - v1(1), v0(2) - v1(2))); };

inline VECTOR3 operator *(double x0, const VECTOR3 & v0)
// return x0*v0
{	return(VECTOR3(x0*v0(0), x0*v0(1), x0*v0(2))); };

inline VECTOR3 operator *(const VECTOR3 & v0, double x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); };

inline VECTOR3 operator /(const VECTOR3 & v0, double x0)
// return v0*x0 (= x0*v0)
{	double x = 1.0 / x0; return VECTOR3( x * v0(0), x * v0(1), x * v0(2) ); };


inline double norm(const VECTOR3& v)
{   return sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2)); }

inline double norm2(const VECTOR3& v)
{   return v(0)*v(0) + v(1)*v(1) + v(2)*v(2);}

inline VECTOR4 get_vector4(VECTOR3 vec)
{	VECTOR4 temp(vec[0], vec[1], vec[2], 1.0); return temp;};


//************************
// VECTOR4 operations
//************************

inline double dot(const VECTOR4 & v0, const VECTOR4 & v1)   
// return dot product of v0 and v1
{	return(v0(0)*v1(0) + v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)); };

inline VECTOR4 operator +(const VECTOR4 & v0, const VECTOR4 & v1)
// return v0 + v1
{	return(VECTOR4(v0(0)+v1(0), v0(1)+v1(1), v0(2)+v1(2), v0(3)+v1(3))); };

inline VECTOR4 operator -(const VECTOR4 & v0, const VECTOR4 & v1)
// return v0 - v1
{	return(VECTOR4(v0(0)-v1(0), v0(1)-v1(1), v0(2)-v1(2), v0(3)-v1(3))); };

inline VECTOR4 operator *(double x0, const VECTOR4 & v0)
// return x0*v0
{	return(VECTOR4(x0*v0(0), x0*v0(1), x0*v0(2), x0*v0(3))); };

inline VECTOR4 operator *(const VECTOR4 & v0, double x0)
// return v0*x0 (= x0*v0)
{	return(x0*v0); };

inline double norm(const VECTOR4& v)
{	return sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2) + v(3)*v(3));	}

inline double norm2(const VECTOR4& v)
{   return v(0)*v(0) + v(1)*v(1) + v(2)*v(2) + v(3)*v(3);}

//************************
// MATRIX3 operations
//************************

MATRIX3 operator +(const MATRIX3 & m0, const MATRIX3 & m1); // return m0 + m1
MATRIX3 operator -(const MATRIX3 & m0, const MATRIX3 & m1); // return m0 - m1
MATRIX3 operator *(const MATRIX3 & m0, const MATRIX3 & m1); // return m0 * m1
MATRIX3 operator *(const double x0, const MATRIX3 & m0);    // return x0 * m0
MATRIX3 operator *(const MATRIX3 & m0, const double x0);    // return m0 * x0
VECTOR3 operator *(const MATRIX3 & m0, const VECTOR3 & v0); // return m0 * v0
VECTOR3 operator *(const VECTOR3 & v0, const MATRIX3 & m0); // return v0 * m0

//************************
// MATRIX4 operations
//************************

MATRIX4 operator +(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 + m1
MATRIX4 operator -(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 - m1
MATRIX4 operator *(const MATRIX4 & m0, const MATRIX4 & m1); // return m0 * m1
MATRIX4 operator *(const double x0, const MATRIX4 & m0);    // return x0 * m0
MATRIX4 operator *(const MATRIX4 & m0, const double x0);    // return m0 * x0
VECTOR4 operator *(const MATRIX4 & m0, const VECTOR4 & v0); // return m0 * v0
VECTOR4 operator *(const VECTOR4 & v0, const MATRIX4 & m0); // return v0 * m0
VECTOR3 operator *(const MATRIX4 & m0, const VECTOR3 & v0); // return m0 * v0
VECTOR3 operator *(const VECTOR3 & v0, const MATRIX4 & m0); // return v0 * m0


MATRIX4 inverse(const MATRIX4 & m);  // return inverse of m; return 0 matrix if
                                     // m is singular
MATRIX4 rotate_matrix(int type, double angle); // type: 1:x, 2:y, 3:z
MATRIX4 translate_matrix(double dx, double dy, double dz);
MATRIX4 scale_matrix(double sx, double sy, double sz);

//******************************
// output procedures, operators
//******************************

ostream & operator<<(ostream & os, const VECTOR3 & v);  // output vector
ostream & operator<<(ostream & os, const VECTOR4 & v);  // output vector
void print_matrix(const MATRIX3 & mat,                  // print matrix
		  const int indent = 0,
		  const int cout_width = 6, const int precision = 1);
void print_matrix(const MATRIX4 & mat,                  // print matrix
		  const int indent = 0,
		  const int cout_width = 6, const int precision = 1);

#endif
