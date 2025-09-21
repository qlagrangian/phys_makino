#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  <stdiostream.h>

main()
{
    double x, v,k, t,h ;
    k = 1; x = 1; v = 0;
    cerr << "Enter h:";
    cin >> h;
    h *= 2*M_PI;
    cout.precision(16);
    cout  << " x= " << x << " v= " << v << " h = " <<  h <<"\n";
    for(t = h; t < 2*M_PI + h/2 ; t+= h ){
	double kx1,kx2,kx3,kx4;
	double kv1,kv2,kv3,kv4;
	kx1 = v*h;
	kv1 = -k*x*h;
	kx2 = (v+kv1/2)*h;
	kv2 = -k*(x+kx1/2)*h;
	kx3 = (v+kv2/2)*h;
	kv3 = -k*(x+kx2/2)*h;
	kx4 = (v+kv3)*h;
	kv4 = -k*(x+kx3)*h;
	x += (kx1+2*kx2+2*kx3+kx4)/6;
	v += (kv1+2*kv2+2*kv3+kv4)/6;
    }
    cout.precision(16);
    cout << "T= " << t << " x= " << x << " v= " << v << "\n";
    cout  << "           err_x= " << x-cos(t)
          << " err_v= " << v+sin(t) << "\n";
    cout  << "  t= " << t << " sin = " << sin(t) << "\n";
}
