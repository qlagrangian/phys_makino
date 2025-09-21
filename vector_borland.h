#ifndef  STARLAB_VECTOR_H
#  define  STARLAB_VECTOR_H
#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

#ifndef real
# define real double
#endif
/*-----------------------------------------------------------------------------
 *  vector  --  a class for 3-dimensional vectors
 *-----------------------------------------------------------------------------
 */
class vector
    {
    private:
        real element[VLEN];

    public:
//	Default: initialize to zero.
        vector(real c = 0)
	    {for (int i = 0; i<VLEN;i++)element[i] = c;}

//  []: the return type is declared as a reference (&), so that it can be used
//  on the left-hand side of an asignment, as well as on the right-hand side,
//  i.e.  v[1] = 3.14  and  x = v[2]  are both allowed and work as expected.

        real & operator [] (int i)       {return element[i];}

	inline void print() {for(int i = 0; i<VLEN;i++)
			     cout << element[i] << " ";
			     cout << "\n";}

//	Unary -

        const vector operator - ()
	    {vector v;for(int i = 0; i<VLEN;i++)v.element[i] = - element[i];
	          return v;}

        friend const vector operator + (const vector &, const vector & );
        friend const vector operator - (const vector &, const vector & );
        friend vector operator * (real, const vector & );
        friend real operator * (const vector &, const vector & ); 
        friend vector operator * (const vector &, real);
        friend vector operator / (const vector &, real);

//	Vector +=, -=, *=, /=

        vector& operator += (const vector& b)
	    {for(int i = 0; i<VLEN;i++)element[i] += b.element[i];       
	     return *this;}

	vector& operator -= (const vector& b)
	    {for(int i = 0; i<VLEN;i++)element[i] -= b.element[i];
	     return *this;}

	vector& operator *= (const real b)
	    {for(int i = 0; i<VLEN;i++)element[i] *= b;
	     return *this;}

	vector& operator /= (const real b)
	    {for(int i = 0; i<VLEN;i++)element[i] /= b;
	     return *this;}
	
	//      Input / Output
	
        friend ostream & operator << (ostream & , const vector & );
	
	friend istream & operator >> (istream & , vector & );
    };

inline  ostream & operator << (ostream & s, const vector & v)
{
    for(int i = 0; i<VLEN;i++)s << v.element[i] << "  " ;
    return s;
}
inline  istream & operator >> (istream & s, vector & v)
{for(int i = 0; i<VLEN;i++)s >> v.element[i]; return s;}

inline  const vector operator + (const vector &v1, const vector & v2)
{
    vector v3;
    for(int i = 0; i<VLEN;i++)v3.element[i] = v1.element[i]+ v2.element[i];
    return v3;
}
inline  const vector operator - (const vector &v1, const vector & v2)
{
    vector v3;
    for(int i = 0; i<VLEN;i++)v3.element[i] = v1.element[i]- v2.element[i];
    return v3;
}

inline  vector operator * (real b, const vector & v)
{
    vector v3;
    for(int i = 0; i<VLEN;i++)v3.element[i] = b* v.element[i];
    return v3;
}
inline  vector operator * (const vector & v, real b)
{
    vector v3;
    for(int i = 0; i<VLEN;i++)v3.element[i] = b* v.element[i];
    return v3;
}

inline  real operator * (const vector & v1, const vector & v2)
{
    real x = 0;
    for(int i = 0; i<VLEN;i++) x+= v1.element[i]* v2.element[i];
    return x;
}


inline  vector operator / (const vector & v, real b)
{
    vector v3;
    for(int i = 0; i<VLEN;i++)v3.element[i] =  v.element[i]/b;
    return v3;
}

typedef vector (*vfunc_ptr)(vector &); 
typedef vector (vfunc)(vector &); 

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\
//  |  the end of:  |         /|\         |  inc/vector.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\
 
