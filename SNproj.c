#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cgs.h"
#include "composition.h"
#include "Teos.h"

#define NVARS 7

#define iRads 0
#define iDens 1
#define iTemp 2
#define iC12 3
#define iO16 4
#define iNe20 5
#define iNe22 6

#define diags 0

/* make a pre-Ia supernova
 * given rhoc, Tc and a file "ref_profile.dat"
 * containing the T, mass fraction profile.
 */

double f_Mcore(int nzones, double *(refData)[]);

int main( int argc, char** argv) {
	char ifname[] = "ref_profile.dat";
        FILE *fin;
        fin = fopen(ifname,"r");
	char header[NVARS+1][20];
        int i;
        for(i=0;i<NVARS+1;i++) {
        	char *temp[20];
                fscanf(fin,"%s",temp);
                strcpy(header[i],temp);
      	}
	int nzones;
        fscanf(fin,"%i",&nzones);	
	double **data;
        data = (double**)malloc(nzones * sizeof(double*));
        for(i=0;i<nzones;i++) {
        	data[i] = (double*)malloc(NVARS * sizeof(double));
        }
        for(i=0;i<nzones;i++) {
        int j;
        for (j=0;j<NVARS;j++) {
        	fscanf(fin,"%lf",&data[i][j]);
           	}
        }
	
	f_Mcore (nzones, data);

	for(i=0;i<nzones;i++) {
        	free(data[i]);
      	}
	fclose(fin);

	return 0;
}


struct root_next_zone_params {
	double P0;
	double rho0;
	double T0;
	double g;
	double delta;
	double T;
	Composition_t *comp;
};

/* function to zero */

double root_f( double rho, void* params) {
	struct root_next_zone_params *p = params;
	double P0=p->P0;
	double T=p->T;
	double rho0=p->rho0;
	double g=p->g;
	double delta=p->delta;
	double Prtx;
	double Pw;
	double f;
	Composition_t *comp = p->comp;

	Prtx = Teos_P(rho,T,comp);
	Pw = P0 + 0.5*delta*g*(rho+rho0);
//	f = (Pw-Prtx)/(P0-Prtx);
	f = Pw - Prtx;
//	f = (Pw - Prtx)/Pw;

	if (diags == 1) {
		printf("delta: %21.15e\n",delta);
		printf("g: %21.15e\n",g);
		printf("rho: %21.15e\n",rho);
		printf("rho0: %21.15e\n",rho0);
		printf("P0: %21.15e\n",P0);
		printf("Prtx: %21.15e\n",Prtx);
		printf("Pw: %21.15e\n",Pw);
		printf("f: %21.15e\n",f);
	}
	return f;
}

double root_df_drho(double rho, void* params) {
	struct root_next_zone_params *p = params;
	struct Teos_values v;
	double g=p->g;
	double P0=p->P0;
	double T=p->T;
	double rho0=p->rho0;
	double delta=p->delta;
	double dfdrho,dpdrho;
	double dpwdrho,dprtxdrho;
	double Prtx;
	double Pw;
	Composition_t *comp = p->comp;

	Teos_all(&v,rho,T,comp);
	dpdrho = (v.chi_rho)*(v.P/rho);
	dfdrho = 0.5*delta*g - dpdrho;

	Prtx = Teos_P(rho,T,comp);
	Pw = P0 + 0.5*delta*g*(rho+rho0);
	dprtxdrho = (v.chi_rho)*(v.P/rho);
	dpwdrho = 0.5*delta*g;
//	dfdrho = (dpwdrho-dprtxdrho)/Pw - (Pw-Prtx)*dpwdrho/(Pw*Pw);
	// Below for: f = (Pw-Prtx)/(P0-Prtx);
	
//	dfdrho = (dpwdrho-dprtxdrho)/(P0-Prtx) - (Pw-Prtx)*(-dprtxdrho)/((P0-Prtx)*(P0-Prtx));	
	if (diags == 1) {
		printf("chi_rho: %21.15e\n",v.chi_rho);
		printf("P: %21.15e\n",v.P);
		printf("rho: %21.15e\n",rho);
		printf("dprtxdrho: %21.15e\n",dprtxdrho);
		printf("dpwdrho: %21.15e\n",dpwdrho);
		printf("delta: %21.15e\n",delta);
		printf("g: %21.15e\n",g);
		printf("dfdrho: %21.15e\n",dfdrho);
	}
	return dfdrho;
}

Composition_t get_zone_comp(int nzones, double *(data)[], int n) {
	Composition_t comp;
	double C12 = data[n][iC12];
	double Ne20 = data[n][iNe20];
	double Ne22 = data[n][iNe22];
	double O16 = data[n][iO16];
        comp.ye = 0.5*(C12 + O16 + Ne20) + 10.0/22.0*Ne22;
        comp.yi = C12/12.0 + O16/16.0 + Ne20/20.0 + Ne22/22.0;
	return comp;
}

int get_sign_double(double d) {
	if (d == 0.0) return 0;
	if (d/fabs(d) < 0.0) {
		return -1;
	} else {
		return 1;
	}
}

/* Build from core out */
double f_Mcore(int nzones, double *(refData)[]) {
	struct root_next_zone_params params;
	struct Teos_values v;
	int max_iter=200;
	double f, dfdrho, deltarho, oldrho;
	double ftol = 1.0e-10;
	double rhoc, Tc, rho, P, T;
	double rcenter, Mrinner, rinner, router;
	double c12, o16, ne20, ne22;
	double deltaP, grhoave, perr;
        Composition_t comp;
	double deltaR = refData[1][iRads]-refData[0][iRads];

	int nz = 0;

	/* start at center */
	/* the central zone is index 0 in refData[] */
	rho = refData[0][iDens];
	rhoc = rho;
	Tc = refData[0][iTemp];
        comp = get_zone_comp(nzones,refData,0);
	P = Teos_P(rhoc,Tc,&comp);
	T = Tc;
//	printf("rhoc: %21.15e\n",rhoc);
//	printf("Pc: %21.15e\n",P);
//	printf("Tc: %21.15e\n",Tc);
//	printf("core ye: %21.15e\n",comp.ye);
//	printf("core yi: %21.15e\n",comp.yi);

	rcenter = 0.5*deltaR;

	// don: print all with pressure diagnostics	
//	printf("# radius density pressure P-P0 0.5g*dr*(rho+rho0) pdiff_erel temperature c12 o16 ne20 ne22\n");
//	printf("%i\n",nzones);
	o16 = 1.0-refData[nz][iC12]-refData[nz][iNe20]-refData[nz][iNe22];	
//	printf("%21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e\n", rcenter, rho, P, 0.0, 0.0, 0.0, T, refData[nz][iC12], o16, refData[nz][iNe20], refData[nz][iNe22]);

//	// don: print all
//	printf("# radius density pressure temperature c12 o16 ne20 ne22\n");
//	printf("%i\n",nzones);
//	o16 = 1.0-refData[nz][iC12]-refData[nz][iNe20]-refData[nz][iNe22];	
//	printf("%21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e\n", rcenter, rho, P, T, refData[nz][iC12], o16, refData[nz][iNe20], refData[nz][iNe22]);
//

	// don: print for flash
	printf("# radius density temperature c12 ne20 ne22\n");
	printf("%i\n",nzones);
	printf("%21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e\n", rcenter, rho, T, refData[nz][iC12], refData[nz][iNe20], refData[nz][iNe22]);

	Mrinner = rhoc*4.0*M_PI/3.0 * deltaR*deltaR*deltaR;
	rcenter = 1.5*deltaR;
	rinner  = deltaR;
	params.rho0 = rhoc;
	nz++;

	do {
		/* set root finding parameters for next zone out */
		params.rho0 = rho;
		params.T0 = T;
		T = refData[nz][iTemp];
		params.T = T;
		params.P0 = P;
		params.delta = deltaR;
		params.g = -cgs_G*Mrinner/rinner/rinner;
                comp = get_zone_comp(nzones,refData,nz);
               	params.comp = &comp;
		
		/* now get the next density */
		int iter=0;
		double Pw,tmp;
		double fl,fh,dfdrl,dfdrh,rl,rh;
		rl = rho;
		rh = rho;
		do { // Newton Solver with density as the independent variable.
			oldrho = rho;
			Pw = params.P0 + 0.5*params.delta*params.g*(rho+params.rho0);
			tmp=0.5*params.delta*params.g*(rho+params.rho0);
			Teos_all(&v,rho,params.T,params.comp);
			P = v.P;
			deltaP = P - params.P0;
                	grhoave = 0.5*(params.g)*deltaR*(params.rho0 + rho);
                	perr = (deltaP - grhoave)/deltaP;
			// Get the value of f
			f = root_f(rho,&params);
			if (iter == 0) {
				fh = f;
				fl = f;
			}
			dfdrho = root_df_drho(rho,&params);
			if (get_sign_double(f) == 1 && (fabs(f) < fabs(fh) || get_sign_double(fh) != 1)) {
				fh = f;
				rh = rho;
				dfdrh = dfdrho;
			} else if (get_sign_double(f) == -1 && (fabs(f) < fabs(fl) || get_sign_double(fl) != -1)) {
				fl = f;
				rl = rho;
				dfdrl = dfdrho;
			} // note if f=0, this is converged!
			// Test if f is low enough in magnitude
			//if (fabs(f) < fabs(ftol*(tmp/(params.P0 + tmp)))) {
	//		if (fabs(f) < ftol) {
			if (fabs(perr) < ftol) {
				// Iteration has converged
				if (diags == 1) {
					printf("Iteration converged in zone: %i\n",nz);
					printf("(f)/Pw: %21.15e\n",fabs(f)/Pw);
					printf("|perr|: %21.15e\n",fabs(perr));
					printf("rho: %21.15e\n",rho);
	//				printf("(f): %21.15e\n",fabs(f));
				}
				break;
			} else {
				// Iteration has not converged, prepare to keep going
				deltarho = -f/dfdrho;
				rho = fmax(0.9*rho,fmin(1.1*rho,rho+deltarho));	
				iter++;
				if (diags == 1) {
					printf("Iteration not converged in zone: %i\n",nz);
					printf("f: %21.15e\n",f);
					printf("new rho/prev rho: %21.15e\n",rho/oldrho);
					printf("deltarho: %21.15e\n",deltarho);
					printf("Next iteration: %i\n",iter);
				}
			}
//			printf("iteration: %i\n",iter);
		} while (iter < max_iter);
		if (iter == max_iter) {	// Newton Solver did not converge, so try something more bruteforce. 
			 printf("Newton Solver did not converge in zone: %i\n",nz); //exit(1); 	
			// Just seek down till you get a pressure balance to tolerance ftol
		//	double fl,fh,dfdrl, dfdrh, rl, rh; // don: use the values from +/- f in the preceding loop.
		//	int sign_dfdrl, sign_dfdrh;
			int jter = 0;
			int max_jter = 1;
			double lowfactor = 0.99999;
		//	rh = params.rho0;
		//	rl = lowfactor*rh;
		//	dfdrh = root_df_drho(rh,&params);
		//	fh = root_f(rh,&params);
		//	sign_dfdrh = get_sign_double(dfdrh);
			do {
				//dfdrl = root_df_drho(rl,&params);
				//sign_dfdrl = get_sign_double(dfdrl);			
				fl = root_f(rl,&params);
				//if (sign_dfdrl == -sign_dfdrh) {
				//if (fabs(fl) > fabs(fh)) {
				if (get_sign_double(fl) == -1 && get_sign_double(fh) == 1) {
					// Stop and subdivide
					int kter = 0;
					int max_kter = 10000000;
					int converged = 0;
					double target;
					double deltaP, grhoave, perr;
			//		double temp_dfdrho;
		//			int sign_temp_dfdrho;
					do {
						rho = 0.5*(rl + rh);
						f = root_f(rho,&params);
						tmp=0.5*params.delta*params.g*(rho+params.rho0);	
//						target = fabs(ftol*(tmp/(params.P0 + tmp)));
//						printf("Target: %21.15e\n",target);
						//target = ftol*1.0e-3;
						Teos_all(&v,rho,params.T,params.comp);
						P = v.P;
						deltaP = P - params.P0;
                				grhoave = 0.5*(params.g)*deltaR*(params.rho0 + rho);
                				perr = (deltaP - grhoave)/deltaP;
						//if (fabs(f) < ftol && fabs(perr) < ftol) {
						//printf("perr: %21.15e\n",perr);
						//printf("ftol: %21.15e\n",ftol);
						if (fabs(perr) < ftol) {
							printf("Subdivision |f| and perr converged to ftol, f: %21.15e, perr: %21.15e\n",f,perr);
							printf("rho: %21.15e\n",rho);
							converged = 1;
							break;
						} else {
//							temp_dfdrho = root_df_drho(rho,&params);
//							sign_temp_dfdrho = get_sign_double(temp_dfdrho);
//							if (sign_temp_dfdrho == sign_dfdrl) {
//								rl = rho;
//								printf("kter: %i, Setting RL=RHO: %21.15e\n",kter,rl);
//							} else if (sign_temp_dfdrho == sign_dfdrh) {
//								rh = rho;
//								printf("kter: %i, Setting RH=RHO: %21.15e\n",kter,rh);
//							} else {
//								printf("kter: %i, Concern: dfdrho zero but f hasn't reached ftol!\n",kter);
//								exit(1);						
//							}	
							// Do something clever here based on fl, fh, and f...
							if (get_sign_double(f) == 1) {
								rh = rho;
								fh = f;
								dfdrh = root_df_drho(rh,&params);
							} else if (get_sign_double(f) == -1) {
								rl = rho;
								fl = f;
								dfdrl = root_df_drho(rl,&params);
							} // note: if get_sign_double(f) = 0, it is converged!	
							kter++;
						}
					} while (kter < max_kter);		
					if (converged == 0) {
						printf("max_kter too low! f=%21.15e\n",f);
					}
					break;
				} else {
					// Keep iterating
					printf("jter: %i, NOTE: fl/fh signs wrong!!\n",jter);
					printf("fl: %21.15e, fh: %21.15e\n",fl,fh);
					rl = lowfactor*rl;
					jter++;
				}			
			} while (jter < max_jter);
			if (jter == max_jter) {printf("jter iteration did not converge in zone: %i\n",nz); exit(1);}
		}

		// update the pressure with the current rho
		Teos_all(&v,rho,T,&comp);
		P = v.P;

		o16 = 1.0-refData[nz][iC12]-refData[nz][iNe20]-refData[nz][iNe22];	
		
		// don: print all with pressure diagnostics
		deltaP = P - params.P0;
		grhoave = 0.5*(params.g)*deltaR*(params.rho0 + rho);
		perr = (deltaP - grhoave)/deltaP;
//		printf("%21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e\n", rcenter, rho, P, deltaP, grhoave, perr, T, refData[nz][iC12], o16, refData[nz][iNe20], refData[nz][iNe22]);
		
		// don: print all
		// printf("%21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e\n", rcenter, rho, P, T, refData[nz][iC12], o16, refData[nz][iNe20], refData[nz][iNe22]);

		// don: print for flash
		printf("%21.15e  %21.15e  %21.15e  %21.15e  %21.15e  %21.15e\n", rcenter, rho, T, refData[nz][iC12], refData[nz][iNe20], refData[nz][iNe22]);
		
		router = rinner+deltaR;

		/* step to next cell */
		Mrinner = Mrinner + rho*4.0*M_PI/3.0*(router*router+router*rinner+rinner*rinner)*(router-rinner);
		rinner = router;
		rcenter = rcenter+deltaR;
		nz++;

	} while (nz < nzones);
	
	return Mrinner;
}


