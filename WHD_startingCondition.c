#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "WHD_startingCondition.h"

int readStartingCondition(StartingCondition *sc,char* filename){	
		
	const char par_names[][20]={"mu_x","mu_y","mu_z","sigma_logx","sigma_logy","sigma_logz"};
	
	double** par_values = (double**) malloc(6*sizeof(double*));
	par_values[0] = &(sc->mu_x);
	par_values[1] = &(sc->mu_y);
	par_values[2] = &(sc->mu_z);
	par_values[3] = &(sc->sigma_logx);
	par_values[4] = &(sc->sigma_logy);
	par_values[5] = &(sc->sigma_logz);
	
	double** par_log_values = (double**) malloc(6*sizeof(double*));
	par_log_values[0] = &(sc->mu_logx);
	par_log_values[1] = &(sc->mu_logy);
	par_log_values[2] = &(sc->mu_logz);
	par_log_values[3] = &(sc->sigma_logx);
	par_log_values[4] = &(sc->sigma_logy);
	par_log_values[5] = &(sc->sigma_logz);

	double** par_values_min = (double**) malloc(6*sizeof(double*));
	par_values_min[0] = &(sc->loop_x_min);
	par_values_min[1] = &(sc->loop_y_min);
	par_values_min[2] = &(sc->loop_z_min);
	par_values_min[3] = &(sc->sigma_logx_min);
	par_values_min[4] = &(sc->sigma_logy_min);
	par_values_min[5] = &(sc->sigma_logz_min);

	double** par_values_max = (double**) malloc(6*sizeof(double*));
	par_values_max[0] = &(sc->loop_x_max);
	par_values_max[1] = &(sc->loop_y_max);
	par_values_max[2] = &(sc->loop_z_max);
	par_values_max[3] = &(sc->sigma_logx_max);
	par_values_max[4] = &(sc->sigma_logy_max);
	par_values_max[5] = &(sc->sigma_logz_max);

	int** par_values_nb = (int**) malloc(6*sizeof(int*));
	par_values_nb[0] = &(sc->loop_x_nb);
	par_values_nb[1] = &(sc->loop_y_nb);
	par_values_nb[2] = &(sc->loop_z_nb);
	par_values_nb[3] = &(sc->sigma_logx_nb);
	par_values_nb[4] = &(sc->sigma_logy_nb);
	par_values_nb[5] = &(sc->sigma_logz_nb);
	
	FILE *partxt;
	partxt = fopen(filename, "r");
	if (partxt == NULL){
	    printf("Error opening file %s!\n",filename);
	    return(-1);
  	}
	
	char line[200];
	char param[100];		
	double buf;

	for(int j=0;j<6;j++){
		sc->count[j]=0;
		rewind(partxt);
		short int test=0;
		while(fgets(line, sizeof line, partxt) != NULL){
		  sscanf(line,"%s : %lf",param,&buf);
			if(strcmp(param,par_names[j])==0){
			  //printf("%s : %f\n",param,buf);
				*(par_values[j])= buf;
				if(j<3){
					*(par_log_values[j])= log(buf);
					} else {
					*(par_log_values[j])= buf;
				}
				*(par_values_min[j])= -1.0;
				*(par_values_max[j])= -1.0;
				*(par_values_nb[j])= 0;
				//printf("%f\n",*par_values[j]);
				//printf("%f\n",*par_log_values[j]);
				test=1;
				break;
			} 
			
		}
		if(test==0){
			printf("ERROR, missing parameter: %s in %s!\n",par_names[j], filename);
			fclose(partxt);
			free(par_values);
			free(par_log_values);
			free(par_values_min);
			free(par_values_max);
			free(par_values_nb);
			return(-1);
		}
	}
	fclose(partxt);	
	free(par_values);
	free(par_log_values);
	free(par_values_min);
	free(par_values_max);
	free(par_values_nb);
	return(0);
}
/*initialize incrementation to next starting condition values if optional arguments were given*/

void printStartingX(StartingCondition sc) {
  printf("mu_x = %g, sigma_logx = %g\n",sc.mu_x,sc.sigma_logx);
}
void printStartingY(StartingCondition sc) {
  printf("mu_y = %g, sigma_logy = %g\n",sc.mu_y,sc.sigma_logy);
}
void printStartingZ(StartingCondition sc) {
  printf("mu_z = %g, sigma_logz = %g\n",sc.mu_z,sc.sigma_logz);
}

void printStarting(StartingCondition sc) {
  printStartingX(sc);
  printStartingY(sc);
  printStartingZ(sc);  
}

int nextStartingCondition(StartingCondition *sc){
  double** par_values = (double**) malloc(6*sizeof(double*));
	par_values[0] = &(sc->mu_x);
	par_values[1] = &(sc->mu_y);
	par_values[2] = &(sc->mu_z);
	par_values[3] = &(sc->sigma_logx);
	par_values[4] = &(sc->sigma_logy);
	par_values[5] = &(sc->sigma_logz);

  double** par_log_values = (double**) malloc(6*sizeof(double*));
	par_log_values[0] = &(sc->mu_logx);
	par_log_values[1] = &(sc->mu_logy);
	par_log_values[2] = &(sc->mu_logz);
	par_log_values[3] = &(sc->sigma_logx);
	par_log_values[4] = &(sc->sigma_logy);
	par_log_values[5] = &(sc->sigma_logz);
  
  double** par_values_min = (double**) malloc(6*sizeof(double*));
	par_values_min[0] = &(sc->loop_x_min);
	par_values_min[1] = &(sc->loop_y_min);
	par_values_min[2] = &(sc->loop_z_min);
	par_values_min[3] = &(sc->sigma_logx_min);
	par_values_min[4] = &(sc->sigma_logy_min);
	par_values_min[5] = &(sc->sigma_logz_min);

  int** par_values_nb = (int**) malloc(6*sizeof(int*));
	par_values_nb[0] = &(sc->loop_x_nb);
	par_values_nb[1] = &(sc->loop_y_nb);
	par_values_nb[2] = &(sc->loop_z_nb);
	par_values_nb[3] = &(sc->sigma_logx_nb);
	par_values_nb[4] = &(sc->sigma_logy_nb);
	par_values_nb[5] = &(sc->sigma_logz_nb);

  double** par_values_incr = (double**) malloc(6*sizeof(double*));
	par_values_incr[0] = &(sc->loop_x_incr);
	par_values_incr[1] = &(sc->loop_y_incr);
	par_values_incr[2] = &(sc->loop_z_incr);
	par_values_incr[3] = &(sc->sigma_logx_incr);
	par_values_incr[4] = &(sc->sigma_logy_incr);
	par_values_incr[5] = &(sc->sigma_logz_incr);

  for(int i=0;i<6;i++) {
    if(*(par_values_nb[i])>1) {
      if(sc->count[i] < *(par_values_nb[i]) - 1) {
	if(sc->logscale[i]==0){
	*(par_values[i]) += *(par_values_incr[i]);
        *(par_log_values[i]) = log(*(par_values[i]));
	} else {
	*(par_log_values[i]) += *(par_values_incr[i]);
	*(par_values[i]) = exp(*(par_log_values[i]));
	}
	sc->count[i]++;

	//printf("%g %g\n",p->alpha,p->gamma);
	free(par_values);
        free(par_log_values);
	free(par_values_min);
	free(par_values_incr);
	free(par_values_nb);
	return 1;
      } else {
	if(sc->logscale[i]==0){
	*(par_values[i]) = *(par_values_min[i]);
	*(par_log_values[i]) = log(*(par_values_min[i]));
	} else {
	*(par_values[i]) = exp(*(par_values_min[i]));
	*(par_log_values[i]) = *(par_values_min[i]);
	}
	sc->count[i]=0;
      }
    }
  }
  free(par_values);
  free(par_log_values);
  free(par_values_min);
  free(par_values_incr);
  free(par_values_nb);
  return 0;
}

int setLoopXYZ(char who,double min,double max,int nb,int logscale,StartingCondition *sc) {
	double incr = (max-min)/((double) nb - 1.0);
	switch(who){
		case 'x':
		sc->mu_x = logscale==1?exp(min):min;
		sc->mu_logx = logscale==1?min:log(min);
		sc->loop_x_min = min;
		sc->loop_x_max = max;
		sc->loop_x_nb = nb;
		sc->loop_x_incr = incr;
		sc->logscale[0] = logscale==1?1:0;
		break;

		case 'y':
		sc->mu_y = logscale==1?exp(min):min;
		sc->mu_logy = logscale==1?min:log(min);
		sc->loop_y_min = min;
		sc->loop_y_max = max;
		sc->loop_y_nb = nb;
		sc->loop_y_incr = incr;
		sc->logscale[1] = logscale==1?1:0;
		break;

		case 'z':
		sc->mu_z = logscale==1?exp(min):min;
		sc->mu_logz = logscale==1?min:log(min);
		sc->loop_z_min = min;
		sc->loop_z_max = max;
		sc->loop_z_nb = nb;
		sc->loop_z_incr = incr;
		sc->logscale[2] = logscale==1?1:0;
		break;

		default:
		printf("Error, cannot set Loop for variable %c",who);
		break;
	}

	
}

int setLoopSigmaXYZ(int who_s,double min_s,double max_s,int nb_s,StartingCondition *sc) {
	double incr_s = (max_s-min_s)/((double) nb_s - 1.0);
	switch(who_s){
		case 'x':
		sc->sigma_logx = min_s;
		sc->sigma_logx_min = min_s;
		sc->sigma_logx_max = max_s;
		sc->sigma_logx_nb = nb_s;
		sc->sigma_logx_incr = incr_s;
		break;

		case 'y':
		sc->sigma_logy = min_s;
		sc->sigma_logy_min = min_s;
		sc->sigma_logy_max = max_s;
		sc->sigma_logy_nb = nb_s;
		sc->sigma_logy_incr = incr_s;
		break;

		case 'z':
		sc->sigma_logz = min_s;
		sc->sigma_logz_min = min_s;
		sc->sigma_logz_max = max_s;
		sc->sigma_logz_nb = nb_s;
		sc->sigma_logz_incr = incr_s;
		
		break;

		default:
		printf("Error, cannot set Loop for variable Sigma_log%c",who_s);
		break;
	}

}

double getXini_gsl(StartingCondition sc,gsl_rng * r){
	if(sc.sigma_logx<=0.0) return(sc.mu_x);
	double res = sc.mu_logx + gsl_ran_gaussian(r, sc.sigma_logx);
	return (exp(res));
}

double getYini_gsl(StartingCondition sc,gsl_rng * r){
	if(sc.sigma_logy<=0.0) return(sc.mu_y);
	double res = sc.mu_logy + gsl_ran_gaussian(r, sc.sigma_logy);
	return (exp(res));
}

double getZini_gsl(StartingCondition sc,gsl_rng * r){
	if(sc.sigma_logz<=0.0) return(sc.mu_z);
	double res = sc.mu_logz + gsl_ran_gaussian(r, sc.sigma_logz);
	return (exp(res));
}
