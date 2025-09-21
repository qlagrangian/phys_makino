#ifndef  NBODY_H
#  define  NBODY_H

#include  "vector.h"

class  body
{
private:
    vector  pos;            // position (3-D Cartesian vector)
    vector  vel;            // velocity: (d/dt) pos
    vector  acc;            // acceleration: (d/dt) vel
    vector  acc_old;        // previous acceleration: (d/dt) vel
    real phi;               //potential
    real mass;              //mass
public:
    body(){pos = vel = acc = 0; }
    void push_pos(double h){pos += h*vel;}
    void push_vel(double h){vel += h*acc;}
    void set_pos(const vector & p ){pos = p;}
    void set_vel(const vector & v ){vel = v;}
    void set_mass(double m ){mass = m;}
    void clear_acc_and_phi(){acc = 0; phi = 0;}
    void trapezoidal_iterate(body & prev_val, double h);
     
    friend void calc_and_acc_force_between_two_particles(body * p1,
							 body * p2,
							 double eps2);
    double kinetic_energy(){return mass*vel*vel*0.5;}
    double potential_energy(){return phi*mass;}

    real & phasevar (int i) {
	if(i<NDIM){
	    return pos[i];
	}else{
	    return vel[i];
	}
    }

    friend ostream & operator << (ostream & , const body & );
    
    friend istream & operator >> (istream & , body & );


    void dump()	{
	cerr << "p= " <<pos << " v=" << vel << "m =" << mass <<"\n";
    }
	    
};

inline  ostream & operator << (ostream & s, const body & b)
{
    s << b.pos << b.vel << b.mass;
    return s;
}
inline  istream & operator >> (istream & s, body & b)
{
    s>> b.pos >> b.vel >> b.mass;
    return s;
}



typedef struct nbody{
    int n; // number of particles;
    real eps2;
    body * bp; // pointer to body array
};


#endif NBODY_H

