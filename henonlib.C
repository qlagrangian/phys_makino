//
// HENONLIB.C
//
// tools for HENONx.C
//

const double EPS = 1e-10;


typedef struct plot_control{
    int xvarid;
    int yvarid;
    int secvarid;
    int mode;
}PLOT_CONTROL;

struct plot_control pctr;

void initgraph(const double xmax, const double h0, const int reuse=0)
{
    if (reuse == 0){
	cerr << "Enter var to plot (x, y, and sos):";
	cin >> pctr.xvarid >>  pctr.yvarid >> pctr.secvarid ;
	if ( pctr.secvarid < 0){
	    pctr.mode = 0;
	}else{
	    pctr.mode = 1;
	}
	char  pgdevname[255];
	cerr << "Enter devname";
	cin >> pgdevname ;
	
	if(cpgopen(pgdevname) != 1) exit(EXIT_FAILURE);
	cpgask(0);
    }else{
	cpgeras();
    }
    cpgenv(-xmax, xmax, -xmax, xmax, 1, 0);
    const char * name[4]={"q1","p1","q2","p2"};
    char label[255];
    sprintf(label, "H0 = %g", h0);
    cpglab(name[pctr.xvarid],name[pctr.yvarid], label);
}


const vector dxdt( vector & x)
{
    vector d;
    d[0] = x[1];
    d[1] = -x[0] - 2* x[0]*x[2];
    d[2] = x[3];
    d[3] = -x[2] - x[0]*x[0] + x[2]*x[2];
    return d;
}


double energy( vector &x)
{
    return 0.5*(x[1]*x[1]+x[3]*x[3]) + 0.5*(x[0]*x[0]+x[2]*x[2])
	+ x[0]*x[0]*x[2]-x[2]*x[2]*x[2]/3;
}

#ifndef WRONG_PLOT
int plot( vector x,
	   vector xold,
	   vfunc_ptr f,
	   const double h)
{

    if (pctr.mode == 0 ){
	cpgpt1(x[pctr.xvarid],x[pctr.yvarid],-1);
	return 0;
    }else    if (x[pctr.secvarid]*xold[pctr.secvarid] <= 0) {
	double hfrac = h;
	double hfold;
	do{
	    hfold = hfrac;
	    hfrac *= -xold[pctr.secvarid]/
		(x[pctr.secvarid]-xold[pctr.secvarid]);
	    x = xold;
	    int first_call = 1;
	    yoshida6b(x,f,hfrac,first_call);
	}while ((fabs((hfrac-hfold)/h)> EPS)&&(hfrac/h> EPS));
	cpgpt1(x[pctr.xvarid],x[pctr.yvarid],-1);
	return 1;
    }
    return 0;


}
void plot( vector x,  vector xold)
{

    if (pctr.mode == 0 ){
	cpgpt1(x[pctr.xvarid],x[pctr.yvarid],-1);
    }else    if (x[pctr.secvarid]*xold[pctr.secvarid] <= 0) {
	vector xinterp = (-x*xold[pctr.secvarid]+xold*x[pctr.secvarid])
	    /(x[pctr.secvarid]-xold[pctr.secvarid]);
	cpgpt1(xinterp[pctr.xvarid],xinterp[pctr.yvarid],-1);
	cout << xinterp[pctr.xvarid] << " " <<xinterp[pctr.yvarid] << "\n"; 
    }	


}
#else
void plot( vector x,  vector xold)
{

    if (pctr.mode == 0 ){
	cpgpt1(x[pctr.xvarid],x[pctr.yvarid],-1);
    }else    if (x[pctr.secvarid]*xold[pctr.secvarid] <= 0) {
	cpgpt1(x[pctr.xvarid],x[pctr.yvarid],-1);
    }	


}
#endif
