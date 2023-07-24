#include <math.h>
#include <algorithm>
#include <iostream>

#include "fermi_integrals.hpp"
#include "../num_tools/num_tools.hpp"
#include "../num_tools/integration/integration.hpp"

namespace aparicio {
//.......Define parameters for GFF and theta derivatives
  double const d=3.3609, sigma=0.091186;
  double const a1=6.7774, b1=1.1418, c1=2.9826;
  double const a2=3.7601, b2=0.093719, c2=0.021063, d2=31.084, e2=1.0056;
  double const a3=7.5669, b3=1.1695, c3=0.75416, d3=6.6558, e3=-0.12819;

//.......Define parameters for first eta derivative
  double const de=4.99551, sigmae=9.11856e-2;
  double const a1e=6.77740, b1e=1.14180, c1e=2.98255;
  double const a2e=3.76010, b2e=9.37188e-2, c2e=2.10635e-2, d2e=3.95015e1, e2e=1.00557;
  double const a3e=7.56690, b3e=1.16953, c3e=7.54162, d3e=7.564734, e3e=-1.28190e-1;
}

using namespace aparicio;

void fermi_func1(double eta, double theta, float k, double a, double b, double f[]){
	double x, tmp;
	for (int i=0;i<ngle;i++) {
		x = (b-a)*xgle[i]+a;

		tmp = x*x-eta;
		tmp = safe_FDexp(tmp);

		f[i] = (b-a)*2.*pow(x,(2.*k+1.))*pow((1.+x*x*theta/2.),0.5)*tmp;
		//printf("x = %.10e, f = %.10e\n", x, f[i]);
	}
	//printf("\n");
	return;
}

void fermi_func2(double eta, double theta, float k, double a, double b, double f[]){
        double x, tmp;
        for (int i=0;i<ngle;i++) {
                x = (b-a)*xgle[i]+a;

                tmp = x-eta;
		tmp = safe_FDexp(tmp);

                f[i] = (b-a)*pow(x,k)*pow((1.+x*theta/2.),0.5)*tmp;
		//printf("x = %.10e, f = %.10e\n", x, f[i]);
        }
	//printf("\n");
        return;
}
         


void fermi_func3(double eta, double theta, float k, double s3, double f[]){
	double x, tmp;
	for (int i=0;i<ngla;i++) {
		x = xgla[i];

		tmp = -x-s3+eta;
		tmp = safe_FDexp(tmp);

		f[i] = exp(eta-s3)*pow(x+s3,k)*pow((1.+(x+s3)*theta/2.),0.5) * tmp; /// (1.+exp(-x-s3+eta));
	}

	if (alpha != 0.) {
		for (int i=0;i<ngla;i++) f[i] = f[i]/pow(xgla[i],alpha);
	}
	return; 
}


void fermi_ed1(double eta, double theta, float k, double a, double b, double f[]) {
        double x, tmp;
        for (int i=0;i<ngle;i++) {
                x = (b-a)*xgle[i]+a;

                tmp = x*x-eta;
                tmp = std::min(tmp,120.);
        	tmp = exp(tmp)/pow(1.+exp(tmp),2.);

                f[i] = (b-a)*2.*pow(x,(2.*k+1.))*pow((1.+x*x*theta/2.),0.5)*tmp;
        }
        return;
}

void fermi_ed2(double eta, double theta, float k, double a, double b, double f[]) {
        double x, tmp;
        for (int i=0;i<ngle;i++) {
                x = (b-a)*xgle[i]+a;

        	tmp = x-eta;
		tmp = std::min(tmp,120.);
		tmp = exp(tmp)/pow(1.+exp(tmp),2.);
		
		f[i] = (b-a)*pow(x,k)*pow((1.+theta*x/2.),0.5)*tmp;
	}
	return;
}


void fermi_ed3(double eta, double theta, float k, double s3, double f[]) {
        double x, tmp, tmp1, tmp2;
        for (int i=0;i<ngla;i++) {
                x = xgla[i];

		tmp1 = 2.*x+s3-eta;
		tmp1 = std::min(tmp1,120.);
		tmp2 = x+s3-eta;
		tmp2 = std::min(tmp2,120.);
		tmp  = exp(tmp1)/pow((1.+exp(tmp2)),2.);

		f[i] = pow((x+s3),k)*pow((1.+(x+s3)*theta/2.),0.5)*tmp;
	}

	if (alpha != 0.) {
		for (int i=0;i<ngla;i++) f[i] = f[i]/pow(xgla[i],alpha);
	}
	return;
}

void extremes(double eta, double* s1, double* s2, double *s3) {
//.......Given the parameter eta the subroutine returns the integration extremes for GFF
          double csi, xa, xb, xc;

          // @WARNING: this works with <math.h>, does not work with <cmath>
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
          
          // @WARNING: this works with <math.h>, does not work with <cmath>
          csi = (1./sigmae)*log(1.+exp(sigmae*static_cast<long double>(eta-de)));

          xa = (a1e+b1e*csi+c1e*csi*csi)/(1.+c1e*csi);
          xb = (a2e+b2e*csi+c2e*d2e*csi*csi)/(1.+e2e*csi+c2e*csi*csi);
          xc = (a3e+b3e*csi+c3e*d3e*csi*csi)/(1.+e3e*csi+c3e*csi*csi);

          *s1 = xa-xb;
          *s2 = xa;
          *s3 = xa+xc;

        return;
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


GFDs compute_all_FD(double eta, double theta) {
  GFDs out;

  // Generalized Fermi-Dirac integrals for leptons
  out.f12 = compute_res(eta, theta, 0.5); // k = 1/2
  out.f32 = compute_res(eta, theta, 1.5); // k = 3/2
  out.f52 = compute_res(eta, theta, 2.5); // k = 5/2

  const double a_eta = - (eta + 2./theta);

  // Generalized Fermi-Dirac integrals for anti-leptons
  out.f12 = compute_res(a_eta, theta, 0.5); // k = 1/2
  out.f32 = compute_res(a_eta, theta, 1.5); // k = 3/2
  out.f52 = compute_res(a_eta, theta, 2.5); // k = 5/2
  
  return out;
}
