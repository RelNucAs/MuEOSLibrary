#include <math.h> //problem in exp of long double when using <cmath>?
#include <algorithm>

#include "glag.hpp"
#include "gleg.hpp"
#include "parameters.hpp"
#include "FD_functions.hpp"
#include "complete_FG.hpp"

using namespace parameters;
using namespace aparicio;

void extremes(double eta, double* s1, double* s2, double *s3) {
//.......Given the parameter eta the subroutine returns the integration extremes for GFF
          double csi, xa, xb, xc;

	  csi = (1./sigma)*log(1.+exp(sigma*(static_cast<long double>(eta-d))));
	  
	  xa = (a1+b1*csi+c1*csi*csi)/(1.+c1*csi);
          xb = (a2+b2*csi+c2*d2*csi*csi)/(1.+e2*csi+c2*csi*csi);
          xc = (a3+b3*csi+c3*d3*csi*csi)/(1.+e3*csi+c3*csi*csi);

          *s1 = xa-xb;
          *s2 = xa;
          *s3 = xa+xc;

	return;
}

void extremes_ed(double eta, double* s1, double* s2, double *s3) {
//.......Given the parameter eta the subroutine returns the integration extremes for the first eta derivative of GFF
          double csi, xa, xb, xc;

          csi = (1./sigmae)*log(1.+exp(sigmae*static_cast<long double>(eta-de)));

          xa = (a1e+b1e*csi+c1e*csi*csi)/(1.+c1e*csi);
          xb = (a2e+b2e*csi+c2e*d2e*csi*csi)/(1.+e2e*csi+c2e*csi*csi);
          xc = (a3e+b3e*csi+c3e*d3e*csi*csi)/(1.+e3e*csi+c3e*csi*csi);

          *s1 = xa-xb;
          *s2 = xa;
          *s3 = xa+xc;
	
	return;
}


double gleg_integration(double f[]){
	double r = 0.;
	for (int i=0;i<ngle;i++){
		r += f[i]*wgle[i];
	}
	return r;
}

double glag_integration(double f[]){
        double r = 0.;
        for (int i=0;i<ngla;i++){
                r += f[i]*wgla[i];
        }
        return r;
}


double compute_res(double eta, double theta, float k) {
        double r1,r2,r3,r4;
        double s1 = 0., s2 = 0., s3 = 0.;
       	double f1[ngle], f2[ngle], f3[ngle], f4[ngla];

        extremes(eta,&s1,&s2,&s3);

	fermi_func1(eta,theta,k,0.,sqrt(s1),f1);
        r1 = gleg_integration(f1);

	fermi_func2(eta,theta,k,s1,s2,f2);
        r2 = gleg_integration(f2);
	
	fermi_func2(eta,theta,k,s2,s3,f3);
        r3 = gleg_integration(f3);
	
	fermi_func3(eta,theta,k,s3,f4);
        r4 = glag_integration(f4);

	//printf("r1 = %.3e, r2 = %.3e, r3 = %.3e, r4 = %.3e\n", r1, r2, r3, r4);
	//for (int i=0;i<ngle;i++) printf("%.3e\n", f4[i]);
	return r1+r2+r3+r4;
}


double compute_res_ed(double eta, double theta, float k) {
        double r1,r2,r3,r4;
        double s1 = 0., s2 = 0., s3 = 0.;
        double f1[ngle], f2[ngle], f3[ngle], f4[ngla];

        extremes_ed(eta,&s1,&s2,&s3);

        fermi_ed1(eta,theta,k,0.,sqrt(s1),f1);
        r1 = gleg_integration(f1);

        fermi_ed2(eta,theta,k,s1,s2,f2);
        r2 = gleg_integration(f2);

        fermi_ed2(eta,theta,k,s2,s3,f3);
        r3 = gleg_integration(f3);

        fermi_ed3(eta,theta,k,s3,f4);
        r4 = glag_integration(f4);

        return r1+r2+r3+r4;
}
        
