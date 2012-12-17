///////////////////////////////////////////////////////////////////
//   Purpose:    vector and matrix class and operations
//   Created by: Matt Camuto
//
//   Modification:
//        Time            Programmer
//        07-09-99        R. Wenger
//        07-15-00        J. Gao
//		  08-09-03		  J. Sun
///////////////////////////////////////////////////////////////////

#include "matrix.h"

#include <iomanip>
#include <math.h>
using namespace std;
int Shift[32];

void init_shift() {
	Shift[0] = 1;
	for (int i=1; i<32; i++) {
		Shift[i] = 1<<i; 
	}
}

int is_bit_set(INT_32 value, int pos){
	if ((pos<1)||(pos>32)) 
		return 0;
	if (value & Shift[pos-1]) {
		return 1;
	}
	return 0;		
}

void set_bit(INT_32 &value, int pos){
	if ((pos<1)||(pos>32)) 
		return;
	value |= Shift[pos-1];
}

//*********************************************************
//
//                INT64 class definitions
//
//*********************************************************

// Set the posth [1, 64] bit as 1
void INT_64::set_bit(int pos){
	if ((pos<1)||(pos>64)) 
		return;
	int flag = (pos>32)?1:0;	
	int move = (pos-1)%32;
	//vec[flag] |= Shift[move];
	if (move==0) vec[flag] |=1; 
	else vec[flag] |= 1<<move;
}

// determine whether the posth(from right to left: 1,2,3...) bit is set
int INT_64::is_set(int pos){
	if ((pos<1)||(pos>64)) 
		return 0;
	int flag = (pos>32)?1:0;	
	int move = (pos-1)%32;
	int temp = 0;
	if (move==0) temp |=1; 
	else temp |= 1<<move;
	if (temp & vec[flag])
		return 1;
	return 0;		
}


//*********************************************************
//
//                VECTOR2 class definitions
//
//*********************************************************
void VECTOR2::Normalize()
{
	double d = vec[0] * vec[0] + vec[1] * vec[1];
	d = sqrt(d);
	if(d < 1e-5)
		return;
	vec[0] /= d;
	vec[1] /= d;
}


//*********************************************************
//
//                VECTOR3 class definitions
//
//*********************************************************


// get the maximum value
double VECTOR3::GetMax() 
{
	double maxval = vec[0];
	if (vec[1] > maxval) maxval = vec[1];
	if (vec[2] > maxval) maxval = vec[2];
	return maxval;
}

// make sure all dimension <=1.0
void VECTOR3::Clamp() 
{
	for (int i = 0; i < Dimension(); i++)
		if (vec[i]>1.0) vec[i] = 1.0;
}

//*********************************************************
//
//                VECTOR4 class definitions
//
//*********************************************************

// normalize vector
// Assumes (pre-condition): vector != (0,0,0,0)
void VECTOR4::Normalize()
{
	double norm = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

	if(norm < 1e-5)
		return;

	for (int i = 0; i < Dimension(); i++)
		vec[i] = vec[i]/norm;
}


//*********************************************************
//
//                MATRIX3 class definitions
//
//*********************************************************

// set matrix to identity matrix
void MATRIX3::Identity()
{
	for (int i = 0; i < Dimension(); i++) {
		mat[i].Zero();
		mat[i][i] = 1.0;
	};
}

void MATRIX3::Zero()
{
    for (int i = 0; i < Dimension(); i++) {
        mat[i].Zero();
    };
}

//*********************************************************
//
//                MATRIX4 class definitions
//
//*********************************************************

// set matrix to identity matrix
void MATRIX4::Identity()
{
	for (int i = 0; i < Dimension(); i++) {
		mat[i].Zero();
		mat[i][i] = 1.0;
	};
}

//************************
// MATRIX3 operations
//************************

// return m0 + m1
MATRIX3 operator +(const MATRIX3 & m0, const MATRIX3 & m1)
{
	MATRIX3 result;

	result[0] = m0(0) + m1(0);
	result[1] = m0(1) + m1(1);
	result[2] = m0(2) + m1(2);

	return(result);
}

// return m0 - m1
MATRIX3 operator -(const MATRIX3 & m0, const MATRIX3 & m1)
{
	MATRIX3 result;

	result[0] = m0(0) - m1(0);
	result[1] = m0(1) - m1(1);
	result[2] = m0(2) - m1(2);

	return(result);
}

// return m0 * m1
MATRIX3 operator *(const MATRIX3 & m0, const MATRIX3 & m1)
{
	MATRIX3 result;

	for (int i = 0; i < m0.Dimension(); i++)
		for (int j = 0; j < m0.Dimension(); j++) {
			result[i][j] = 0;
		for (int k = 0; k < m0.Dimension(); k++)
			result[i][j] += m0(i,k) * m1(k,j);
		};

	return(result);
}

// return x0 * m0
MATRIX3 operator *(const double x0, const MATRIX3 & m0)
{
	MATRIX3 result;

	result[0] = x0*m0(0);
	result[1] = x0*m0(1);
	result[2] = x0*m0(2);

	return(result);
}

// return m0 * x0
MATRIX3 operator *(const MATRIX3 & m0, const double x0)
{ return(x0*m0); };

// return m0 * v0
VECTOR3 operator *(const MATRIX3 & m0, const VECTOR3 & v0)
{
	VECTOR3 result;

	result[0] = dot(m0(0),v0);
	result[1] = dot(m0(1),v0);
	result[2] = dot(m0(2),v0);

	return(result);
}

// return v0 * m0
VECTOR3 operator *(const VECTOR3 & v0, const MATRIX3 & m0)
{
	VECTOR3 result;

	result[0] = v0(0)*m0(0,0) + v0(1)*m0(1,0) + v0(2)*m0(2,0);
	result[1] = v0(0)*m0(0,1) + v0(1)*m0(1,1) + v0(2)*m0(2,1);
	result[2] = v0(0)*m0(0,2) + v0(1)*m0(1,2) + v0(2)*m0(2,2);

	return(result);
}

//************************
// MATRIX4 operations
//************************

// get the rotation matrix
// type: 1:x, 2:y, 3:z
MATRIX4 rotate_matrix(int type, double angle) {
    MATRIX4 rot_mat;
	switch(type) {
		case 1: // rotate by x
		{
			double anglex = (angle/180.0f)*3.14159f;
			double sinx = (double)sin((double)anglex);
			double cosx = (double)cos((double)anglex);
			rot_mat[1][1] = rot_mat[2][2] = cosx;
			rot_mat[1][2] = -sinx;
			rot_mat[2][1] = sinx;
            break;
		}
		case 2: // rotate by y
		{
			double angley = (angle/180.0f)*3.14159f;
			double siny = (double)sin((double)angley);
			double cosy = (double)cos((double)angley);
			rot_mat[0][0] = rot_mat[2][2] = cosy;
			rot_mat[0][2] = siny;
			rot_mat[2][0] = -siny;
			break;
		}
		case 3: // rotate by z
		{
			double anglez = (angle/180.0f)*3.14159f;
			double sinz = (double)sin((double)anglez);
			double cosz = (double)cos((double)anglez);
			rot_mat[0][0] = rot_mat[1][1] = cosz;
			rot_mat[0][1] = -sinz;
			rot_mat[1][0] = sinz;
			break;
		}
	}
	return rot_mat;
}

// get translation matrix
MATRIX4 translate_matrix(double tx, double ty, double tz){
	MATRIX4 translation_mat;
	translation_mat[0][3] = tx;
	translation_mat[1][3] = ty;
	translation_mat[2][3] = tz;
	return translation_mat;
}

// get scaling matrix
MATRIX4 scale_matrix(double sx, double sy, double sz){
	MATRIX4 scale_mat;
	scale_mat[0][0] = sx;
	scale_mat[1][1] = sy;
	scale_mat[2][2] = sz;
	return scale_mat;
}

// return m0 + m1
MATRIX4 operator +(const MATRIX4 & m0, const MATRIX4 & m1)
{
	MATRIX4 result;

	result[0] = m0(0) + m1(0);
	result[1] = m0(1) + m1(1);
	result[2] = m0(2) + m1(2);
	result[3] = m0(3) + m1(3);

	return(result);
}

// return m0 - m1
MATRIX4 operator -(const MATRIX4 & m0, const MATRIX4 & m1)
{
	MATRIX4 result;

	result[0] = m0(0) - m1(0);
	result[1] = m0(1) - m1(1);
	result[2] = m0(2) - m1(2);
	result[3] = m0(3) - m1(3);

	return(result);
}

// return m0 * m1
MATRIX4 operator *(const MATRIX4 & m0, const MATRIX4 & m1)
{
	MATRIX4 result;

	for (int i = 0; i < m0.Dimension(); i++)
		for (int j = 0; j < m0.Dimension(); j++) {
			result[i][j] = 0;
		for (int k = 0; k < m0.Dimension(); k++)
			result[i][j] += m0(i,k) * m1(k,j);
		};

	return(result);
}

// return x0 * m0
MATRIX4 operator *(const double x0, const MATRIX4 & m0)
{
	MATRIX4 result;

	result[0] = x0*m0(0);
	result[1] = x0*m0(1);
	result[2] = x0*m0(2);
	result[3] = x0*m0(3);

	return(result);
}

// return m0 * x0
MATRIX4 operator *(const MATRIX4 & m0, const double x0)
{ return(x0*m0); };

// return m0 * v0
VECTOR4 operator *(const MATRIX4 & m0, const VECTOR4 & v0)
{
	VECTOR4 result;

	result[0] = dot(m0(0),v0);
	result[1] = dot(m0(1),v0);
	result[2] = dot(m0(2),v0);
	result[3] = dot(m0(3),v0);

	return(result);
}


// return m0 * v0
VECTOR3 operator *(const MATRIX4 & m0, const VECTOR3 & v0)
{
	VECTOR4 v(v0);
	VECTOR3 result;	

	double temp = dot(m0(3),v);
	result[0] = dot(m0(0),v)/temp;
	result[1] = dot(m0(1),v)/temp;
	result[2] = dot(m0(2),v)/temp;
	return(result);
}


VECTOR4 operator *(const VECTOR4 & v0, const MATRIX4 & m0)
// return v0 * m0
{
	VECTOR4 result;

	result[0] = v0(0)*m0(0,0) + v0(1)*m0(1,0) + v0(2)*m0(2,0) + v0(3)*m0(3,0);
	result[1] = v0(0)*m0(0,1) + v0(1)*m0(1,1) + v0(2)*m0(2,1) + v0(3)*m0(3,1);
	result[2] = v0(0)*m0(0,2) + v0(1)*m0(1,2) + v0(2)*m0(2,2) + v0(3)*m0(3,2);
	result[3] = v0(0)*m0(0,3) + v0(1)*m0(1,3) + v0(2)*m0(2,3) + v0(3)*m0(3,3);

	return(result);
}

VECTOR3 operator *(const VECTOR3 & v0, const MATRIX4 & m0)
// return v0 * m0
{
	VECTOR3 result;
	double temp = v0(0)*m0(0,3) + v0(1)*m0(1,3) + v0(2)*m0(2,3) + m0(3,3);

	result[0] = (v0(0)*m0(0,0) + v0(1)*m0(1,0) + v0(2)*m0(2,0) + m0(3,0))/temp;
	result[1] = (v0(0)*m0(0,1) + v0(1)*m0(1,1) + v0(2)*m0(2,1) + m0(3,1))/temp;
	result[2] = (v0(0)*m0(0,2) + v0(1)*m0(1,2) + v0(2)*m0(2,2) + m0(3,2))/temp;

	return(result);
}

//Code was taken from the original 'edge' library written by dave ebert.
MATRIX4 inverse(const MATRIX4 & m) {
	register int lp,i,j,k;
	static double wrk[4][8];
	static double a, b;
	MATRIX4 result;
  
	for( i=0; i<4; i++ )	/* Set up matrices */
	{
		for( j=0; j<4; j++ )
		{
			wrk[i][j]=(double)m(i,j);
			wrk[i][j+4]=0.0;
			result[i][j] = 0.0;
		}
		wrk[i][i+4]=1.0;
    }
  
	for( lp=0; lp<4; lp++ )	/* Loop over all rows */
	{
		a=0.0;
		j=(-1);
		for( i=lp; i<4; i++ )	/* Find largest non-zero element */
		{
			b=wrk[i][lp];
			if( b< 0.0 )
				b=(-b);
			if( b>a )
			{
				a=b;
				j=i;
			}
		}
		if( j!=lp )			/* If not on diagonal, put it there */
		{
			if( j<0 )		/* Singular if none found */
				return(result);
			else			/* Exchange rows from lp to end */
				for( k=lp; k<8; k++ )
			{	
				a=wrk[j][k];
				wrk[j][k]=wrk[lp][k];
				wrk[lp][k]=a;
			}
		}
		a=wrk[lp][lp];		/* Normalize working row */
		for( i=lp; i<8; i++ )
			wrk[lp][i]/=a;
      
		for( i=lp+1; i<8; i++ )  /* Adjust rest of work space */
		{
			b=wrk[lp][i];
			for( j=0; j<4; j++ )	/* One column at a time */
				if( j!=lp )
					wrk[j][i]-=wrk[j][lp]*b;
		}
    }

	for( i=0; i<4; i++ )	/* Return result matrix */
		for( j=0; j<4; j++ )
			result[i][j]=(double)wrk[i][j+4];
	return(result);
}

//*******************************
// output procedures/operators
//*******************************

// output vector
ostream & operator<<(ostream & os, const VECTOR3 & v)
{
	os << "( " << v(0) << " " << v(1) << " " << v(2) << " )";
	return(os);
}

// output vector
ostream & operator<<(ostream & os, const VECTOR4 & v)
{
	os << "( " << v(0) << " " << v(1) << " " << v(2) << " " << v(3) << " )";
	return(os);
}

// routine to print out matrix.  Useful for debugging.
//   call at the beginning of a new line for proper formatting
// mat = matrix
// indent = # of indented spaces
// cout_width = output field width
// precision = output precision
void print_matrix(const MATRIX3 & mat,
		  const int indent, const int cout_width, const int precision)
{
	cout.precision(precision);
	cout.flags(ios::fixed);
	for (int i = 0; i < mat.Dimension(); i++) {
        int j;
		for (j = 0; j < indent; j++)
			cout << " ";
		cout << "( ";
		for (j = 0; j < mat.Dimension(); j++)
			cout << setw(cout_width) << mat(i,j) << " ";
		cout << ")" << endl;
	};
}

// routine to print out matrix.  Useful for debugging.
//   call at the beginning of a new line for proper formatting
// mat = matrix
// indent = # of indented spaces
// cout_width = output field width
// precision = output precision
void print_matrix(const MATRIX4 & mat,
		  const int indent, const int cout_width, const int precision)
{
  cout.precision(precision);
  cout.flags(ios::fixed);
  for (int i = 0; i < mat.Dimension(); i++) {
    int j;
    for (j = 0; j < indent; j++)
      cout << " ";
    cout << "( ";
    for (j = 0; j < mat.Dimension(); j++)
      cout << setw(cout_width) << mat(i,j) << " ";
    cout << ")" << endl;
  }
}

