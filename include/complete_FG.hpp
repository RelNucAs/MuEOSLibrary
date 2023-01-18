#pragma once

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

//.......Given the parameter eta the subroutine returns the integration extremes for GFF
void extremes(double eta, double* s1, double* s2, double *s3);

void extremes_ed(double eta, double* s1, double* s2, double *s3);

double gleg_integration(double f[]);

double glag_integration(double f[]);

double compute_res(double eta, double theta, float k);

double compute_res_ed(double eta, double theta, float k);
        
