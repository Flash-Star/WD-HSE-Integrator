#ifndef TEOS_H
#include "composition.h"

struct Teos_values {
	double P;
	double Gamma;
	double e; // internal energy /g
	double c_V;
	double c_P;
	double chi_T;
	double chi_rho;
	double del_ad;
	double sound;
};

double Teos_P(double rho, double T, Composition_t *comp);

int Teos_all( struct Teos_values *out, double rho, double T,
		Composition_t *comp);
double Teos_rho(double P, double T, Composition_t *comp);

#endif /* TEOS_H */
