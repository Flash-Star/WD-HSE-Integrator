#include <math.h>
#include "cgs.h"
#include "composition.h"
#include "Teos.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
/* interface for Timmes "helmholtz" EIP eos. */

/* this must match the declarations in vector_eos.dek */
#define NROWMAX 1000

struct {
	double temp_row[NROWMAX];
	double  den_row[NROWMAX];
	double abar_row[NROWMAX];
	double zbar_row[NROWMAX];
	double zeff_row[NROWMAX];
} thinp_ ;

struct {
	int jlo_eos;
	int jhi_eos;
} eosvec2_ ;

struct {
	double ptot_row[NROWMAX];
	double  dpt_row[NROWMAX];
	double  dpd_row[NROWMAX];
	double  dpa_row[NROWMAX];
	double  dpz_row[NROWMAX];
	double etot_row[NROWMAX];
	double  det_row[NROWMAX];
	double  ded_row[NROWMAX];
	double  dea_row[NROWMAX];
	double  dez_row[NROWMAX];
	double stot_row[NROWMAX];
	double  dst_row[NROWMAX];
	double  dsd_row[NROWMAX];
	double  dsa_row[NROWMAX];
	double  dsz_row[NROWMAX];
} thtot_ ;

struct {
	double   cp_row[NROWMAX];
	double   cv_row[NROWMAX];
	double gam1_row[NROWMAX];
	double gam2_row[NROWMAX];
	double gam3_row[NROWMAX];
	double   cs_row[NROWMAX];
} thder_ ;

struct {
	double  pcou_row[NROWMAX];
	double  ecou_row[NROWMAX];
	double  scou_row[NROWMAX];
	double plasg_row[NROWMAX];
} thcou_ ;

int helmeos_(void);

double Teos_P(double rho, double T, Composition_t *comp) {
	
	thinp_.temp_row[0] = T;
	thinp_.den_row[0] = rho;
	thinp_.abar_row[0] = 1.0/comp->yi;
	thinp_.zbar_row[0] = comp->ye*thinp_.abar_row[0];

	eosvec2_.jlo_eos=1;
	eosvec2_.jhi_eos=1;

	/* call it */
	helmeos_();
	return thtot_.ptot_row[0];
}

int Teos_all( struct Teos_values *out, double rho, double T,
		Composition_t *comp) {
	/* this fills the common blocks */
	Teos_P(rho,T,comp);
	
	out->P       = thtot_.ptot_row[0];
	out->Gamma   = thcou_.plasg_row[0];
	out->e       = thtot_.etot_row[0];; // internal energy /g
	out->c_V     = thder_.cv_row[0];
	out->c_P     = thder_.cp_row[0];
	out->chi_T   = T/out->P * thtot_.dpt_row[0];
	out->chi_rho = rho/out->P * thtot_.dpd_row[0];
	out->del_ad  = (thder_.gam3_row[0]-1.0)/thder_.gam1_row[0];
	out->sound   = thder_.cs_row[0];

	return 0;
}

struct findrho_params {
	double P;
	double T;
	Composition_t *comp;
};

double findrho_f_Pdiff( double rho, void *params) {
	struct findrho_params *p = params;
	struct Teos_values v;

	Teos_all( &v, rho, p->T, p->comp);
	return (v.P - p->P);
}

double findrho_df_Pdiff( double rho, void *params) {
	struct findrho_params *p = params;
	struct Teos_values v;

	Teos_all( &v, rho, p->T, p->comp);
	return (v.chi_rho*rho/v.P);
}

void findrho_fdf_Pdiff( double rho, void *params, double *y, double *df) {
	struct findrho_params *p = params;
	struct Teos_values v;

	Teos_all( &v, rho, p->T, p->comp);
	*y = v.P - p->P;
	*df = v.chi_rho*rho/v.P;
	return;
}


double Teos_rho(double P, double T, Composition_t *comp) {
	double rhoguess,rhoguessnd,rhoguessd,mue,mu,Prad;
	int status;
	int iter=0, max_iter=100;
	const gsl_root_fsolver_type *type=gsl_root_fsolver_brent;
	gsl_root_fsolver *solver;
	double rhohi,rholo;
	gsl_function F;
	struct findrho_params params = {P,T,comp};

	F.function = &findrho_f_Pdiff;
	F.params = &params;

	/* compute guess */
	Prad=cgs_a/3*T*T*T*T;
        /*
	mue = 2.0/(1+comp->mass_fraction[1]);
        mu = 1/(1/mue+1/comp->mui);
        */
        mue = 1.0/comp->ye;
        mu = 1.0/(comp->ye+comp->yi);
        if (P <= 1e23) { /* non-relativistic */
                //if (pow(T,3.0/2)*6e-9 > P*mu*cgs_amu/cgs_k/T/mue) {
                        /* non-degenerate */
                        rhoguessnd = (P-Prad)*mu*cgs_amu/cgs_k/T;
                //} else {
                        /* degenerate */
                        rhoguessd = mue*pow(P/1.004e13,3.0/5);
               // }
        } else { /* relativistic */
                //if (pow(T,3)*4.6e-24 > P*mu*cgs_amu/cgs_k/T/mue) {
                        /* non-degenerate */
                        rhoguessnd = P*mu*cgs_amu/cgs_k/T;
                //} else {
                        /* degenerate */
                        rhoguessd = mue*pow(P/1.243e15,3.0/4);
                //}
        }
	if (rhoguessnd < rhoguessd ) rhoguess = rhoguessnd;
        else rhoguess = rhoguessd;
//	fprintf(stderr, "%e %e rhoguess = %e\n",P,T, rhoguess);
        rhohi = rhoguess*1e1;
	rholo= rhoguess/1e2;
	
	solver = gsl_root_fsolver_alloc(type);
	gsl_root_fsolver_set(solver,&F,rholo,rhohi);
	do {
		iter++;
		status = gsl_root_fsolver_iterate(solver);
		rholo = gsl_root_fsolver_x_lower(solver);
		rhohi = gsl_root_fsolver_x_upper(solver);
		//fprintf(stderr, "%e %e\n", rholo, rhohi);
		status = gsl_root_test_interval(rholo,rhohi,0,1e-8);
	} while (status == GSL_CONTINUE && iter < max_iter);

	return gsl_root_fsolver_root(solver);
}

