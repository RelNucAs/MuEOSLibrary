#ifndef FERMI_INTEGRALS_HPP
#define FERMI_INTEGRALS_HPP

/*============================================================================*/

// file: generalized_FD.cpp

struct GFDs {
  double f12;
  double f32;
  double f52;
  double a_f12;
  double a_f32;
  double a_f52;
};
typedef struct GFDs GFDs;

// Calculation of Generalized Fermi-Dirac integrals

void fermi_func1(double eta, double theta, float k, double a, double b, double f[]);

void fermi_func2(double eta, double theta, float k, double a, double b, double f[]);

void fermi_func3(double eta, double theta, float k, double s3, double f[]);

void fermi_ed1(double eta, double theta, float k, double a, double b, double f[]);

void fermi_ed2(double eta, double theta, float k, double a, double b, double f[]);

void fermi_ed3(double eta, double theta, float k, double s3, double f[]);

void extremes(double eta, double* s1, double* s2, double *s3);

void extremes_ed(double eta, double* s1, double* s2, double *s3);

double compute_res(double eta, double theta, float k);

double compute_res_ed(double eta, double theta, float k);

/*============================================================================*/

// file: nonrel_FD.cpp

// Calculation of non-relativistic Fermi-Dirac integrals

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -9/2 */
double Fermi_integral_m92(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -7/2 */
double Fermi_integral_m72(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k=-5/2 */
double Fermi_integral_m52(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -3/2 */
double Fermi_integral_m32(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = -1/2 */
double Fermi_integral_m12(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 0 */
double Fermi_integral_0(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 1/2 */
double Fermi_integral_p12(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 1 */
double Fermi_integral_p1(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 3/2 */
double Fermi_integral_p32(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 2 */
double Fermi_integral_p2(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 5/2 */
double Fermi_integral_p52(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 3 */
double Fermi_integral_p3(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 7/2 */
double Fermi_integral_p72(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 4 */
double Fermi_integral_p4(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 9/2 */
double Fermi_integral_p92(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 5 */
double Fermi_integral_p5(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 11/2 */
double Fermi_integral_p112(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 6 */
double Fermi_integral_p6(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 13/2 */
double Fermi_integral_p132(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 7 */
double Fermi_integral_p7(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 15/2 */
double Fermi_integral_p152(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 8 */
double Fermi_integral_p8(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 17/2 */
double Fermi_integral_p172(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 9 */
double Fermi_integral_p9(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 19/2 */
double Fermi_integral_p192(const double x);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 10 */
double Fermi_integral_p10(const double y);

/* double precision rational minimax approximation of non-relativistic */
/* Fermi-Dirac integral of order k = 21/2 */
double Fermi_integral_p212(const double x);

/*============================================================================*/

#endif //FERMI_INTEGRALS_HPP
