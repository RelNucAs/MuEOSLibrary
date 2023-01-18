#pragma once

struct FD_params { float k; double eta; double theta; double a; double b; };

double safe_FDexp(double x);

void fermi_func1(double eta, double theta, float k, double a, double b, double f[]);

void fermi_func2(double eta, double theta, float k, double a, double b, double f[]);

void fermi_func3(double eta, double theta, float k, double s3, double f[]);

void fermi_ed1(double eta, double theta, float k, double a, double b, double f[]);

void fermi_ed2(double eta, double theta, float k, double a, double b, double f[]);

void fermi_ed3(double eta, double theta, float k, double s3, double f[]);
