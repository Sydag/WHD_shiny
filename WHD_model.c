#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_cdf.h>

//#include <R.h>
//#include <Rmath.h>

#include "WHD_model.h"

void initializeParameter(Parameter* p) {
  p->par_names = (char**) malloc(14*sizeof(char*));
  for(int i=0;i<14;i++) p->par_names[i] = (char*) malloc(20*sizeof(char));    
  strcpy(p->par_names[0],"alpha");
  strcpy(p->par_names[1],"beta");
  strcpy(p->par_names[2],"gamma");
  strcpy(p->par_names[3],"delta");
  strcpy(p->par_names[4],"eta");
  strcpy(p->par_names[5],"ksi");
  strcpy(p->par_names[6],"phi");
  strcpy(p->par_names[7],"psi");
  strcpy(p->par_names[8],"omega");
  strcpy(p->par_names[9],"u");
  strcpy(p->par_names[10],"v");
  strcpy(p->par_names[11],"w");
  strcpy(p->par_names[12],"l");
  strcpy(p->par_names[13],"zcrit");
  
  p->par_values = (double**) malloc(14*sizeof(double*));
  p->par_values[0] =&(p->alpha);
  p->par_values[1] =&(p->beta);
  p->par_values[2]=&(p->gamma);
  p->par_values[3]=&(p->delta);
  p->par_values[4]=&(p->eta);
  p->par_values[5]=&(p->ksi);
  p->par_values[6]=&(p->phi);
  p->par_values[7]=&(p->psi);
  p->par_values[8]=&(p->omega);
  p->par_values[9]=&(p->u);
  p->par_values[10]=&(p->v);
  p->par_values[11]=&(p->w);
  p->par_values[12]=&(p->l);
  p->par_values[13]=&(p->zcrit);

  p->par_min_values = (double**) malloc(14*sizeof(double*));
  p->par_min_values[0] =&(p->alpha_min);
  p->par_min_values[1] =&(p->beta_min);
  p->par_min_values[2]=&(p->gamma_min);
  p->par_min_values[3]=&(p->delta_min);
  p->par_min_values[4]=&(p->eta_min);
  p->par_min_values[5]=&(p->ksi_min);
  p->par_min_values[6]=&(p->phi_min);
  p->par_min_values[7]=&(p->psi_min);
  p->par_min_values[8]=&(p->omega_min);
  p->par_min_values[9]=&(p->u_min);
  p->par_min_values[10]=&(p->v_min);
  p->par_min_values[11]=&(p->w_min);
  p->par_min_values[12]=&(p->l_min);
  p->par_min_values[13]=&(p->zcrit_min);
  
  p->par_max_values = (double**) malloc(14*sizeof(double*));
  p->par_max_values[0] =&(p->alpha_max);
  p->par_max_values[1] =&(p->beta_max);
  p->par_max_values[2]=&(p->gamma_max);
  p->par_max_values[3]=&(p->delta_max);
  p->par_max_values[4]=&(p->eta_max);
  p->par_max_values[5]=&(p->ksi_max);
  p->par_max_values[6]=&(p->phi_max);
  p->par_max_values[7]=&(p->psi_max);
  p->par_max_values[8]=&(p->omega_max);
  p->par_max_values[9]=&(p->u_max);
  p->par_max_values[10]=&(p->v_max);
  p->par_max_values[11]=&(p->w_max);
  p->par_max_values[12]=&(p->l_max);
  p->par_max_values[13]=&(p->zcrit_max);

  p->par_nb_values = (int**) malloc(14*sizeof(int*));
  p->par_nb_values[0] =&(p->alpha_nb);
  p->par_nb_values[1] =&(p->beta_nb);
  p->par_nb_values[2]=&(p->gamma_nb);
  p->par_nb_values[3]=&(p->delta_nb);
  p->par_nb_values[4]=&(p->eta_nb);
  p->par_nb_values[5]=&(p->ksi_nb);
  p->par_nb_values[6]=&(p->phi_nb);
  p->par_nb_values[7]=&(p->psi_nb);
  p->par_nb_values[8]=&(p->omega_nb);
  p->par_nb_values[9]=&(p->u_nb);
  p->par_nb_values[10]=&(p->v_nb);
  p->par_nb_values[11]=&(p->w_nb);
  p->par_nb_values[12]=&(p->l_nb);
  p->par_nb_values[13]=&(p->zcrit_nb);
  
  p->par_incr_values = (double**) malloc(14*sizeof(double*));
  p->par_incr_values[0] =&(p->alpha_incr);
  p->par_incr_values[1] =&(p->beta_incr);
  p->par_incr_values[2]=&(p->gamma_incr);
  p->par_incr_values[3]=&(p->delta_incr);
  p->par_incr_values[4]=&(p->eta_incr);
  p->par_incr_values[5]=&(p->ksi_incr);
  p->par_incr_values[6]=&(p->phi_incr);
  p->par_incr_values[7]=&(p->psi_incr);
  p->par_incr_values[8]=&(p->omega_incr);
  p->par_incr_values[9]=&(p->u_incr);
  p->par_incr_values[10]=&(p->v_incr);
  p->par_incr_values[11]=&(p->w_incr);
  p->par_incr_values[12]=&(p->l_incr);
  p->par_incr_values[13]=&(p->zcrit_incr);

  p->count = (int*) malloc(14*sizeof(int));
  
  for(int i=0;i<14;i++) {
    p->count[i];
    *(p->par_nb_values[i])=0;
  }

  p->initialized = 1;
}

void freeParameter(Parameter* p) {
  for(int i=0;i<14;i++) free(p->par_names[i]);
  free(p->par_names);
  free(p->par_values);
  free(p->par_min_values);
  free(p->par_max_values);
  free(p->par_nb_values);
  free(p->par_incr_values);
  free(p->count);    
}

int readParameter(Parameter *p,char* filename){
  initializeParameter(p);
  
  FILE *partxt;
  partxt = fopen(filename, "r");
  if (partxt == NULL){
    printf("Error opening file %s!\n",filename);
    return(-1);
  }

  char line[200];
  char param[100];		
  double buf;

  for(int j=0;j<14;j++){
    p->count[j]=0;
    rewind(partxt);
    short int test=0;
    while(fgets(line, sizeof line, partxt) != NULL){
      sscanf(line,"%s : %lf",param,&buf);
      if(strcmp(param,p->par_names[j])==0){
	printf("%s = %f\n",param,buf); //vérif
	*(p->par_values[j])= buf;
	*(p->par_min_values[j])= buf;
	*(p->par_max_values[j])= buf;
	*(p->par_nb_values[j])= 1;
	*(p->par_incr_values[j])= 0.0;
	test=1;
	break;
      } 			
    }
    if(test==0){
      printf("ERROR, missing parameter: %s in %s!\n",p->par_names[j], filename);
      fclose(partxt);	
      return(-1);
    }
  }
  fclose(partxt);

  //printf("%f %f %f\n",p->alpha,p->beta,p->psi);
  return(0);
}

int  copyParameter(Parameter from,Parameter* to) {
  initializeParameter(to);
  for(int j=0;j<14;j++){
    to->count[j]=0;    
    *(to->par_values[j])     = *(from.par_values[j]);
    *(to->par_min_values[j]) = *(from.par_min_values[j]);
    *(to->par_max_values[j]) = *(from.par_max_values[j]);
    *(to->par_nb_values[j])  = *(from.par_nb_values[j]);
    *(to->par_incr_values[j])= *(from.par_incr_values[j]);
  }
  return 0;
}

int readParameterRange(Parameter* p,char* filename){
  if(p->initialized != 1) {
    printf("ERROR Parameter should be initialized before range is read!\n");
    return(-1);
  }
  
  FILE *partxt = fopen(filename, "r");
  if (partxt == NULL){
    printf("Error opening file %s!\n",filename);
    return(-1);
  }
    
  char line[200];
  char parname[100];		
  double bufmin,bufmax;
  int bufnb,cmpt=0;	

  while(fgets(line, sizeof line, partxt) != NULL){
    if(sscanf(line,"%s : %lf %lf %d",parname,&bufmin,&bufmax,&bufnb)==4) {
      for(int j=0;j<14;j++) {
	if(strcmp(parname,p->par_names[j])==0){
	  *(p->par_min_values[j]) = bufmin;
	  *(p->par_values[j])     = bufmin;
	  *(p->par_max_values[j]) = bufmax;
	  *(p->par_nb_values[j])  = bufnb;
	  if(bufnb>1) {
	    cmpt++;
	    *(p->par_incr_values[j])  = (bufmax-bufmin)/((double) bufnb - 1.0);
	  } else {
	    *(p->par_incr_values[j])  = 0.0;
	  }
	  printf("%s varied from %g to %g by increment of %g\n",parname,*(p->par_min_values[j]),*(p->par_max_values[j]),*(p->par_incr_values[j]));
	  break;
	}
      }
    }
  }

  fclose(partxt);	
  return(cmpt);
}

int nextParameter(Parameter* p) {
  for(int i=0;i<14;i++) {
    if(*(p->par_nb_values[i])>1) {
      if(p->count[i] < *(p->par_nb_values[i]) - 1) {
	*(p->par_values[i]) += *(p->par_incr_values[i]);
	p->count[i]++;
	//printf("%g %g\n",p->alpha,p->gamma);
	return 1;
      } else {
	*(p->par_values[i]) = *(p->par_min_values[i]);
	p->count[i]=0;
      }
    }
  }
  return 0;
}

void printVariableParameter(FILE* f,Parameter p) {
  for(int i=0;i<14;i++) {
    if(*(p.par_nb_values[i])>1) {
      if(!f) {
	printf("%s;",p.par_names[i]);
      } else {	
	fprintf(f,"%s;",p.par_names[i]);
      }
    }
  }
  fflush(f);
}

void printParameterState(FILE* f,Parameter p) {
  for(int i=0;i<14;i++) {
    if(*(p.par_nb_values[i])>1) {
      if(!f) {
	printf("%f;",*(p.par_values[i]));
      } else {
	fprintf(f,"%f;",*(p.par_values[i]));
      }
    }
  }
  fflush(f);
}

void printParameter(FILE* f,Parameter p) {
  for(int i=0;i<14;i++) {
    if(!f) {
      printf("%s=%f\n",p.par_names[i],*(p.par_values[i]));
    } else {
      fprintf(f,"%s=%f\n",p.par_names[i],*(p.par_values[i]));
    }    
  }
  fflush(f);
}

/*functionnal responses and their derivatives*/
double F(double x, Parameter p){
  //double res = (p.gamma+p.alpha*pow(x,p.u))/(1.0+p.beta*pow(x,p.v));
  double res = p.gamma + (p.alpha) * pow(x,p.u)/(1.0+(p.beta)*pow(x,p.v));
  return (res);
}

double G(double y, double z, Parameter p){
  double res = 1.0/(1.0 +pow(y,p.w) + (p.psi)*pow(z,p.l));
  return (res);
}

double dFdx(double x,Parameter p){
  //yields nan if x=0 and u<1 or v<1
  double D = p.beta*pow(x,p.v)+1.0;
  double B = p.alpha*pow(x,(p.u-1.0))/D;
  double res = B*(p.u-p.beta*p.v*pow(x,p.v)/D);
  //printf("%g %g %g\n",D,B,res);
  return (res);
}

double dGdy(double y,double z, Parameter p){
  double D = p.psi*pow(z,p.l)+pow(y,p.w)+1.0;
  double res = -p.w*pow(y,p.w-1.0)/(D*D);
  return (res);
}

double dGdz(double y,double z, Parameter p){
  if(p.psi<=0) return 0.0;
  //yields nan if z=0 and l<1
  double D = p.psi*pow(z,p.l)+pow(y,p.w)+1.0;
  double res = -p.l*p.psi*pow(z,p.l-1.0)/(D*D);
  return (res);
}

/*differential equations*/
double dxdt(double x, double y,double z, Parameter p){
	(void)(z);
	double res =  x*(1.0-x-p.delta*y);
	return res;
}
double dlogxdt(double x, double y,double z, Parameter p){
	(void)(z);
	double res =  1.0-x-p.delta*y;
	return res;
}
double dydt(double x, double y, double z, Parameter p){
	double res = F(x,p)*G(y,z,p)-p.phi*y;
	return res;
}
double dlogydt(double x, double y, double z, Parameter p){
	double res = F(x,p)*G(y,z,p)/y-p.phi;
	return res;
}
double dzdt(double x, double y, double z, Parameter p){
	double res = p.omega*x+p.eta*y-p.ksi*z;
	return res;
}
double dlogzdt(double x, double y, double z, Parameter p){
        double res = (p.omega*x+p.eta*y)/z-p.ksi;
	return res;
}
double d2xdt2(double x, double y, double z, Parameter p){
  double dx = dxdt(x,y,z,p);
  double dy = dydt(x,y,z,p);

  double res = dx*(1.0-x-p.delta*y)-x*(dx+p.delta*dy);
  return res;
}
double d2ydt2(double x, double y, double z, Parameter p){
  double f = F(x,p);
  double g = G(y,z,p);
  double df = dFdx(x,p);
  double dgdy = dGdy(y,z,p);
  double dgdz = dGdz(y,z,p);
  double dx = dxdt(x,y,z,p);
  double dy = dydt(x,y,z,p);
  double dz = dzdt(x,y,z,p);

  double res = dx*df*g + f*(dy*dgdy+dz*dgdz) - dy*p.phi;
  return (res);
}
double d2zdt2(double x, double y, double z, Parameter p){
  double dx = dxdt(x,y,z,p);
  double dy = dydt(x,y,z,p);
  double dz = dzdt(x,y,z,p);

  double res = p.omega*dx+p.eta*dy-p.ksi*dz;
  return (res);
}


//functions for numerical integration
// 1. dx/dt , dy/t, dz/dt
int dxyz_dt (double t, const double y[], double f[],
      void *param)
{
  (void)(t); /* avoid unused parameter warning */
  Parameter* p = (Parameter*) param;
  f[0] = dxdt(y[0],y[1],y[2],*p);
  f[1] = dydt(y[0],y[1],y[2],*p);
  f[2] = dzdt(y[0],y[1],y[2],*p);

  return GSL_SUCCESS;
}

/*jacobian matrix*/
int compute_jacob_dxyz_dt(const double y[], double* dfdy, int complete_system,Parameter p) {
  //in either of following conditions dFdx, dGdy or dGdz will be nan
  if(y[0]<=0.0 && (p.u<1.0 || p.v<1.0)) return -1;
  if(y[1]<=0.0 && p.w<1.0) return -1;
  if(complete_system==1 && p.psi>0 && y[2]<=0.0 && p.l<1.0) return -1;

  //x
  dfdy[0] = 1.0-p.delta*y[1]-2.0*y[0];
  dfdy[1] = -p.delta*y[0];
  if(complete_system!=0) dfdy[2] = 0.0;
  //y  
  if(complete_system!=0) {
    dfdy[3] = dFdx(y[0],p)*G(y[1],y[2],p);
    dfdy[4] = F(y[0],p)*dGdy(y[1],y[2],p)-p.phi;
    dfdy[5] = F(y[0],p)*dGdz(y[1],y[2],p);
  } else {
    dfdy[2] = dFdx(y[0],p)*G(y[1],0.0,p);
    dfdy[3] = F(y[0],p)*dGdy(y[1],0.0,p)-p.phi;    
  }
  //z
  if(complete_system!=0) {
    dfdy[6] = p.omega;
    dfdy[7] = p.eta;
    dfdy[8] = -p.ksi;
  }
  
  return 0;
}

int jacob_dxyz_dt (double t, const double y[], double *dfdy,
     double dfdt[], void *param)
{
  (void)(t); /* avoid unused parameter warning */
  Parameter* p = (Parameter*) param;

  if(compute_jacob_dxyz_dt(y,dfdy,1,*p)!=0) return GSL_FAILURE;
  
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 3, 3);
  gsl_matrix * m = &dfdy_mat.matrix;
  
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  return GSL_SUCCESS;
}

/*integrating along z instead of t*/
double dxdz(double x, double y,double z, Parameter p) {
	double res = dxdt(x,y,z,p)/dzdt(x,y,z,p);
	return (res);
}

double dydz(double x, double y,double z, Parameter p) {
	double res = dydt(x,y,z,p)/dzdt(x,y,z,p);
	return (res);
}
double dtdz(double x, double y,double z, Parameter p) {
	double res = 1.0/dzdt(x,y,z,p);
	return (res);
}

//functions for numerical integration
// 1. dx/dt , dy/t, dz/dt
int dxyt_dz (double z, const double y[], double f[],
      void *param)
{
  
  Parameter* p = (Parameter*) param;
  f[0] = dxdt(y[0],y[1],z,*p);
  f[1] = dydt(y[0],y[1],z,*p);
  f[2] = 1.0/dzdt(y[0],y[1],z,*p);
  f[0]*=f[2];
  f[1]*=f[2];
  
  return GSL_SUCCESS;
}

/*jacobian matrix*/
int jacob_dxyt_dz (double z, const double y[], double *dfdy,
     double dfdz[], void *param)
{
  Parameter* p = (Parameter*) param;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 3, 3);
  gsl_matrix * m = &dfdy_mat.matrix;

  double dx = dxdt(y[0],y[1],z,*p);
  double dy = dydt(y[0],y[1],z,*p);
  double dz = dzdt(y[0],y[1],z,*p);

  double ddxdx = 1.0-2.0*y[0]-p->delta*y[1];
  double ddxdy = -p->delta*y[0];
  double ddxdz = 0.0; 

  double ddydx = dFdx(y[0],*p)*G(y[1],z,*p);
  double ddydy = F(y[0],*p)*dGdy(y[1],z,*p)-p->phi;
  double ddydz = F(y[0],*p)*dGdz(y[1],z,*p);

  double ddzdx = p->omega;
  double ddzdy = p->eta;
  double ddzdz = -p->ksi;
  
  double dz2 = pow(dz,2.0);

  gsl_matrix_set (m, 0, 0, (ddxdx*dz-dx*ddzdx)/dz2);
  gsl_matrix_set (m, 0, 1, (ddxdy*dz-dx*ddzdy)/dz2);
  gsl_matrix_set (m, 0, 2, 0.0);	

  gsl_matrix_set (m, 1, 0, (ddydx*dz-dy*ddzdx)/dz2);
  gsl_matrix_set (m, 1, 1, (ddydy*dz-dy*ddzdy)/dz2);
  gsl_matrix_set (m, 1, 2, 0.0);

  gsl_matrix_set (m, 2, 0, -ddzdx/dz2);
  gsl_matrix_set (m, 2, 1, -ddzdy/dz2);
  gsl_matrix_set (m, 2, 2, 0.0);

  dfdz[0] = (ddxdz*dz-dx*ddzdz)/dz2;
  dfdz[1] = (ddydz*dz-dy*ddzdz)/dz2;
  dfdz[2] = -ddzdz/dz2;
  return GSL_SUCCESS;
}

// 3. dlogx/dt , dy/t, dz/dt
int dlogxyz_dt (double t, const double y[], double f[],
      void *param)
{
  (void)(t); /* avoid unused parameter warning */
  Parameter* p = (Parameter*) param;
  double x = exp(y[0]);
  f[0] = dlogxdt(x,y[1],y[2],*p);
  f[1] = dydt(x,y[1],y[2],*p);
  f[2] = dzdt(x,y[1],y[2],*p);

  return GSL_SUCCESS;
}

int jacob_dlogxyz_dt (double t, const double y[], double *dfdy,
     double dfdt[], void *param)
{
  Parameter* p = (Parameter*) param;
  gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 3, 3);
  gsl_matrix * m = &dfdy_mat.matrix;

  double x  = exp(y[0]);
  
  double ddxdx = 1.0-2.0*x-p->delta*y[1];
  double ddxdy = -p->delta*y[0];
  double ddxdz = 0.0; 

  double ddydx = dFdx(x,*p)*G(y[1],y[2],*p);
  double ddydy = F(x,*p)*dGdy(y[1],y[2],*p)-p->phi;
  double ddydz = F(x,*p)*dGdz(y[1],y[2],*p);

  double ddzdx = p->omega;
  double ddzdy = p->eta;
  double ddzdz = -p->ksi;
  
  gsl_matrix_set (m, 0, 0, ddxdx*x);
  gsl_matrix_set (m, 0, 1, ddxdy);
  gsl_matrix_set (m, 0, 2, 0.0);	

  gsl_matrix_set (m, 1, 0, ddydx*x);
  gsl_matrix_set (m, 1, 1, ddydy);
  gsl_matrix_set (m, 1, 2, ddydz);

  gsl_matrix_set (m, 2, 0, ddzdx*x);
  gsl_matrix_set (m, 2, 1, ddzdy);
  gsl_matrix_set (m, 2, 2, ddzdz);
  
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  dfdt[2] = 0.0;
  return GSL_SUCCESS;
}

/*equilibrium x value are roots of H*/
double H(double x, void *vp){
  Parameter* p = (Parameter*) vp;
  
  double y = (1.0-x)/p->delta;
  
  double A = (p->alpha)*pow(x,p->u);
  double B = 1.0+pow(x,p->v)*(p->beta);
  double C = 1.0 + pow(y,p->w);
  if(p->psi>0) {
    //problem here when ksi is zero : equilibrium cannot be computed
    //if psi>0 and ksi=0
    double z = ((p->omega)*x+(p->eta)*y)/(p->ksi);
    C += (p->psi)*pow(z,p->l);
  }
  double D = A + (p->gamma - p->phi*y*C)*B;
  return(D);
}

double Hlogx(double lx, void *vp){
  return H(exp(lx),vp);
}

double dH(double x, void *vp){
  Parameter* p = (Parameter*) vp;

  double y = (1.0-x)/p->delta;
  double dy = -1.0/p->delta;

  double dA = (p->alpha)*(p->u)*pow(x,p->u-1.0);
  double B  = 1.0+pow(x,p->v)*(p->beta);
  double dB = (p->v)*pow(x,p->v-1.0)*(p->beta);
  double C  = 1.0 + pow(y,p->w);
  double dC = (p->w)*dy*pow(y,p->w-1.0);  
  if(p->psi>0) {
    //problem here when ksi is zero : equilibrium cannot be computed
    //if psi>0 and ksi=0
    double z  = (p->omega*x+p->eta*y)/p->ksi;
    double dz = (p->omega+p->eta*dy)/p->ksi;
    C  += p->psi*pow(z,p->l);
    dC += p->psi*dz*(p->l)*pow(z,p->l-1.0);
  }
  
  //double D = A + (p->gamma - p->phi*y*C)*B;
  double dD = dA + dB*(p->gamma-p->phi*y*C) - B*(p->phi)*(dy*C+y*dC);
  
  return(dD);
}

double dHlogx(double lx, void *vp){
  double x = exp(lx);
  return x*dH(x,vp);
}

void HdH(double x, void *vp, double* h,double *dh){
  Parameter* p = (Parameter*) vp;  

  double y = (1.0-x)/p->delta;
  double dy = -1.0/p->delta;
  double C = 1.0 + pow(y,p->w);
  double dC = p->w*dy*pow(y,p->w-1.0);
  
  if(p->psi>0) {
    //problem here when ksi is zero : equilibrium cannot be computed
    //if psi>0 and ksi=0
    double z  = (p->omega*x+p->eta*y)/p->ksi;
    double dz = (p->omega+p->eta*dy)/p->ksi;
    C  += p->psi*pow(z,p->l);
    dC += p->psi*dz*(p->l)*pow(z,p->l-1.0);
  }
  
  double A  = p->alpha*pow(x,p->u);
  double dA = p->alpha*p->u*pow(x,p->u-1.0);
  double B  = pow(x,p->v)*p->beta+1.0;
  double dB = p->v*pow(x,p->v-1.0)*p->beta;
  
  double D  = A + (p->gamma-p->phi*y*C)*B;
  double dD = dA + dB*(p->gamma-p->phi*y*C) - B*p->phi*(dy*C+y*dC);
  
  *h = D;
  *dh = dD;
}

void HdHlogx(double lx, void *vp, double* h,double *dh){
  double x = exp(lx);
  HdH(x,vp,h,dh);
  *dh *= x;
}

/* gamma value for a given equilirium x */
double P_gamma_eq(double x, void* vp) {
  Parameter p = *((Parameter*) vp);
  
  double y = (1.0-x)/p.delta;
  double A = 1.0+pow(y,p.w);
  if(p.psi>0) {
    double z = ((p.omega)*x+(p.eta)*y)/p.ksi;
    A += (p.psi)*pow(z,p.l);
  }
  double B = F(x,p) - (p.gamma);
  
  double geq = (p.phi)*y*A - B;
  return(geq);
}

double P_gamma_leq(double lx, void* vp) {
  return P_gamma_eq(exp(lx),vp);
}

double dP_gamma_eq(double x, void* vp) {
  Parameter p = *((Parameter *) vp);

  double y = (1.0-x)/(p.delta);
  double dy = -1.0/(p.delta);
  double A = 1.0+pow(y,p.w);
  double dA = (p.w)*dy*pow(y,p.w-1.0);
  if(p.psi>0) {
    double z = ((p.omega)*x+(p.eta)*y)/p.ksi;
    double dz = ((p.omega) + dy*(p.eta))/p.ksi;
    A += (p.psi)*pow(z,p.l);
    dA += (p.psi)*(p.l)*dz*pow(z,(p.l)-1.0);
  }
  double B = F(x,p)- (p.gamma);
  double dB = dFdx(x,p);

  double dgeq = (p.phi)*(dy*A+y*dA)-dB;
  return(dgeq); //why -dgeq
}

double dP_gamma_leq(double lx, void* vp) {
  return dP_gamma_eq(exp(lx),vp);
}


double P_alpha_eq(double x, void* vp) {

  Parameter p = *((Parameter*) vp);

  double y = (1.0-x)/p.delta;
  double A = 1.0+pow(y,p.w);
  if(p.psi>0) {
    double z = ((p.omega)*x+(p.eta)*y)/p.ksi;
    A += (p.psi)*pow(z,p.l);
  }
  A *= (p.phi)*y;
  A -= p.gamma;
  
  double B = (1.0+p.beta*pow(x,p.v))/pow(x,p.u);
  
  double a_eq = A*B;
  return(a_eq);
}

double P_alpha_leq(double lx, void* vp) {
  return P_alpha_eq(exp(lx),vp);
}

double dP_alpha_eq(double x, void* vp) {
  Parameter p = *((Parameter*) vp);
  

  double y  = (1.0-x)/p.delta;
  double dy = -1.0/p.delta;
  
  double A = 1.0+pow(y,p.w);  
  double z;
  if(p.psi>0) {
    z = ((p.omega)*x+(p.eta)*y)/p.ksi;
    A += (p.psi)*pow(z,p.l);
  }
  
  double dA = p.phi*dy*A + p.w*p.phi*dy*pow(y,p.w);
  if(p.psi>0) {
    double dz = (p.omega+p.eta*dy)/p.ksi;
    dA += p.l*p.psi*p.phi*y*dz*pow(z,p.l-1);
  }

  A *= (p.phi)*y;
  A -= p.gamma;

  double B = (1.0+p.beta*pow(x,p.v))/pow(x,p.u);
  double dB = p.beta*(p.v-p.u)*pow(x,p.v-p.u-1.0) - p.u*pow(x,-p.u-1.0);

  double da_eq = dB*A + B*dA;

  return(da_eq);
}

double dP_alpha_leq(double lx, void* vp) {
  return dP_alpha_eq(exp(lx),vp);
}

//critical gamma above which clearance is stable
double gamma_clearance(Parameter p) {
  double gc =  1 + 1/pow(p.delta,p.w);
  if(p.psi>0 && p.ksi>0) gc += p.psi*pow((p.eta/(p.delta*p.ksi)),p.l);
  gc *= p.phi/p.delta;
  return gc;
}

//critical gamma above which the host dies with no infection
double gamma_survival(Parameter p) {
  double yeq = p.ksi*p.zcrit/p.eta;
  double gc = p.phi * yeq * ( 1.0 + pow(yeq,p.w) + p.psi*pow(p.zcrit,p.l));
  return gc;
}

//critical gamma below which an infection will always kill
//because equilibrium z is above zcrit
double gamma_kill(Parameter p) {
  double x = p.eta/p.delta;
  x = (p.ksi*p.zcrit-x)/(p.omega-x); //x
  double y = (1.0-x)/p.delta;
  if(x<0 || y<0) return(NAN);
  double A = p.phi*y*(1.0+pow(y,p.w)+p.psi*pow(p.zcrit,p.l));
  double B = p.alpha*pow(x,p.u)/(1.0+p.beta*pow(x,p.v));
  return(A-B);
}

//critical alpha  below which an infection will always kill
//because equilibrium z is above zcrit
double alpha_kill(Parameter p) {
  double x = p.eta/p.delta;
  x = (p.ksi*p.zcrit-x)/(p.omega-x); //x
  double y = (1.0-x)/p.delta;
  if(x<0 || y<0) return(NAN);
  double A = p.phi*y*(1.0+pow(y,p.w)+p.psi*pow(p.zcrit,p.l))-p.gamma;
  double B = pow(x,p.u)/(1.0+p.beta*pow(x,p.v));
  return(A/B);
}

//find range of gamma values in which the system is bistable
//maxg is the maximum number of cirtical values (shoud be 2)
//curv??
int gamma_crit (int maxg,double* g_crit,int* curv,Parameter p) {
  for(int i=0;i<maxg;i++) g_crit[i]=-1.0;
  
  if(p.ksi==0 && p.psi>0.0) {
    #ifdef RMOD
    error("ksi = ",p.ksi," and psi=",p.psi,". Calculation of critical gamma stopped because division by zero will occur.\n");
    #else
    printf("ksi = %g and psi = %g. Calculation of critical gamma stopped because division by zero will occur.\n",p.ksi,p.psi);
    #endif    
    return(-1);
  }

  void gsl_permisive_error_handler(const char * reason,
				      const char * file,
				      int line,
				      int gsl_errno) {
    printf("\nError while computing critical gamma: %s\n", gsl_strerror (gsl_errno));
    printf("%s\n\n",reason);
  }
  gsl_error_handler_t*  old_handler = gsl_set_error_handler (&gsl_permisive_error_handler);

  int status;
  int iter = 0, max_iter = 100;

  gsl_function FdP;
  FdP.function = &dP_gamma_leq;
  FdP.params = &p;

  
  //find zeros of dP
  const gsl_root_fsolver_type *Tr = gsl_root_fsolver_brent;
  gsl_root_fsolver *sr = gsl_root_fsolver_alloc (Tr);
  
  double r = -1.0;
  double lxmin = -50.0,lxmax=0.0,lxincr=1e-2;
  double lx_lo,lx_hi;
  int    pos = 0;
  for(double lx=lxmin;lx<lxmax & pos<maxg;lx+=lxincr) {
    if(dP_gamma_leq(lx,&p)*dP_gamma_leq(lx+lxincr,&p)<0) {
      gsl_root_fsolver_set (sr, &FdP,lx,lx+lxincr);
      iter=0;
      do
	{
	  iter++;
	  status = gsl_root_fsolver_iterate (sr);
	  r = gsl_root_fsolver_root (sr);
	  lx_lo = gsl_root_fsolver_x_lower (sr);
	  lx_hi = gsl_root_fsolver_x_upper (sr);
	  status = gsl_root_test_interval (lx_lo, lx_hi,0,0.001);
	  //printf("%d gamma_lo=%g\n",iter,r);
	}
      while (status == GSL_CONTINUE && iter < max_iter);
  
      if (status == GSL_SUCCESS) {
	//printf("%d %g %g\n",pos,r,P_gamma_eq(exp(r),&p));
	g_crit[pos] = P_gamma_eq(exp(r),&p);
	
	double dp_lo = dP_gamma_leq(lx,&p);
	double dp_hi = dP_gamma_leq(lx+lxincr,&p);
	curv[pos]=0;
	if(dp_lo>0 & dp_hi<0) {
	  curv[pos] = 1; //r is a maximum
	}
	if(dp_lo<0 & dp_hi>0) {
	  curv[pos] = -1; //r is a minimum
	}
	pos++;
      }
    }
  }
  gsl_root_fsolver_free (sr);
  gsl_set_error_handler (old_handler);
  //returns the number of roots; should be two if the system is bistable
  return pos;
}

//find range of alpha values in which the system is bistable
//maxa is the maximum number of cirtical values (shoud be 2)
int alpha_crit (int maxa,double* a_crit,int* curv,Parameter p) {
  for(int i=0;i<maxa;i++) a_crit[i]=-1.0;
  
  if(p.ksi==0 && p.psi>0.0) {
    #ifdef RMOD
    error("ksi = ",p.ksi," and psi=",p.psi,". Calculation of critical alpha stopped because division by zero will occur.\n");
    #else
    printf("ksi = %g and psi = %g. Calculation of critical alpha stopped because division by zero will occur.\n",p.ksi,p.psi);
    #endif    
    return(-1);
  }

  void gsl_permisive_error_handler(const char * reason,
				      const char * file,
				      int line,
				      int gsl_errno) {
    printf("\nError while computing critical alpha: %s\n", gsl_strerror (gsl_errno));
    printf("%s\n\n",reason);
  }
  gsl_error_handler_t*  old_handler = gsl_set_error_handler (&gsl_permisive_error_handler);

  int status;
  int iter = 0, max_iter = 100;

  gsl_function FdP;
  FdP.function = &dP_alpha_leq;
  FdP.params = &p;

  
  //find zeros of dP
  const gsl_root_fsolver_type *Tr = gsl_root_fsolver_brent;
  gsl_root_fsolver *sr = gsl_root_fsolver_alloc (Tr);
  
  double r = -1.0;
  double lxmin = -50.0,lxmax=0.0,lxincr=1e-2;
  double lx_lo,lx_hi;
  int    pos = 0;
  for(double lx=lxmin;lx<lxmax & pos<maxa;lx+=lxincr) {
    if(dP_alpha_leq(lx,&p)*dP_alpha_leq(lx+lxincr,&p)<0) {
      gsl_root_fsolver_set (sr, &FdP,lx,lx+lxincr);
      iter=0;
      do
	{
	  iter++;
	  status = gsl_root_fsolver_iterate (sr);
	  r = gsl_root_fsolver_root (sr);
	  lx_lo = gsl_root_fsolver_x_lower (sr);
	  lx_hi = gsl_root_fsolver_x_upper (sr);
	  status = gsl_root_test_interval (lx_lo, lx_hi,0,0.001);
	  //printf("%d alpha_lo=%g\n",iter,r);
	}
      while (status == GSL_CONTINUE && iter < max_iter);
  
      if (status == GSL_SUCCESS) {
	//printf("%d %g %g\n",pos,r,P_alpha_leq(r,&p));
	a_crit[pos] = P_alpha_leq(r,&p);
	
	double dp_lo = dP_alpha_leq(lx,&p);
	double dp_hi = dP_alpha_leq(lx+lxincr,&p);
	curv[pos]=0;
	if(dp_lo>0 & dp_hi<0) {
	  curv[pos] = 1; //r is a maximum
	}
	if(dp_lo<0 & dp_hi>0) {
	  curv[pos] = -1; //r is a minimum
	}
	pos++;
      }
    }
  }
  gsl_root_fsolver_free (sr);
  gsl_set_error_handler (old_handler);
  //returns the number of roots; should be two if the system is bistable
  return pos;
}



//computation of equilibrium y in the absence of any infection
double dy_NoInf(double y, void* vp) {
  Parameter p = *((Parameter*) vp);
   double dy =  1.0 + pow(y,p.w);
   if(p.psi>0 && p.eta>0) {
     double zeq = y*p.eta/p.ksi;
     dy+= p.psi*pow(zeq,p.l);
   }
   dy *= p.phi*y;
   dy -= F(0,p);
   return(dy);
}

int NoInfeq(double *yeq,Parameter p) {  
  if(p.gamma<=0) {
    *yeq = 0.0;
    return(1);
  }

  *yeq = -1.0;
  const gsl_root_fsolver_type *Tr = gsl_root_fsolver_brent;
  gsl_root_fsolver *sr = gsl_root_fsolver_alloc (Tr);
    
  gsl_function Fdy;
  Fdy.function = &dy_NoInf;
  Fdy.params = &p;

  int max_iter = 100,iter=0,status;
  double ymax = p.gamma/p.phi;  
  double ymin = 1.0/p.delta;
  if(ymin!=ymin || ymin>=ymax) {
    ymin = 0.0;
  } else {
    if(dy_NoInf(ymin,&p)*dy_NoInf(0.0,&p)<0) {
      ymax= ymin;
      ymin=0.0;
    }
  }
  double y=ymin;
  if(dy_NoInf(0,&p)*dy_NoInf(ymax,&p)<0) {
    gsl_root_fsolver_set (sr, &Fdy, ymin, ymax);
    iter=0;
    //printf("%d yeq=%g\n",iter,y);
    do
      {
	iter++;
	status = gsl_root_fsolver_iterate (sr);
	y = gsl_root_fsolver_root (sr);
	double y_lo = gsl_root_fsolver_x_lower (sr);
	double y_hi = gsl_root_fsolver_x_upper (sr);
	status = gsl_root_test_interval (y_lo,y_hi,0,0.001);
	//printf("%d yeq=%g\n",iter,y);
      }
    while (status == GSL_CONTINUE && iter < max_iter);
    gsl_root_fsolver_free (sr);
    
    if (status == GSL_SUCCESS) {
      *yeq = y;      
      return(1);
    }
  }
  
  return(0);
}

int setYZ(int n,double xtmp,double* xeq,double* yeq,double* zeq,Parameter p) {
  //invalid value of x
  if(xtmp<0 || xtmp>1.0) return n;
  //is x an equilibrium point??
  if(fabs(H(xtmp,&p))>1e-3) return n;
    
  //is equilibrium point already stored in xeq?
  int test = 0;
  for(int i=0;i<n && test==0;i++) {
    if(xtmp==0 && xeq[i]==0) {
      test=1;
    } else {
      //log transform x values to increse precision
      if(fabs(log(xtmp)-log(xeq[i]))<1e-4) test=1;
    }
  }
  if(test==1) return n; //equilibrium point is not new!

  //sort equilibrium points in ascending order according to xeq
  int i=0;
  while(xeq[i]<xtmp && i<n) i++;
  if(i<n) {
    for(int j=n-1;j>=i;j--) {
      xeq[j+1]=xeq[j];
      yeq[j+1]=yeq[j];
      zeq[j+1]=zeq[j];
    }
  }
  
  //store new equilibrium
  xeq[i] = xtmp;
  yeq[i] = (1.0 - xeq[i])/p.delta;
  if(p.ksi>0) {
    zeq[i] = (p.omega*xeq[i]+p.eta*yeq[i])/p.ksi;
  } else {
    //z has no equilibrium value
    if(xeq[i]<=0.0 && p.eta<=0.0) {
      //pathogens have been cleared and defenses cause no damage
      //z final value is finite and depends on initial x value
      zeq[i] =  NAN;
    } else {
      //z goes to infinity
      zeq[i] = INFINITY;
    }
  }

  //return new number of equilibria
  return ++n;  
}

double** equilibrium(int *nbeq,Parameter p) {
  int nbeq_max = 10;
  int maxiter = 1000;
  
  double *xeq = malloc(nbeq_max*sizeof(double));
  double *yeq = malloc(nbeq_max*sizeof(double));
  double *zeq = malloc(nbeq_max*sizeof(double));

  for(int i=0;i<nbeq_max;i++) {
    xeq[i] = -1.0;
    yeq[i] = -1.0;
    zeq[i] = -1.0;
  }
  double xtmp=-1.0,ytmp=-1.0,ztmp=-1.0;
  *nbeq = 0;
  
  void gsl_permisive_error_handler(const char * reason,
				      const char * file,
				      int line,
				      int gsl_errno) {
    printf("\nError while computing equilibrium: %s\n", gsl_strerror (gsl_errno));
    printf("%s\n\n",reason);
  }
  gsl_error_handler_t*  old_handler = gsl_set_error_handler (&gsl_permisive_error_handler);

  //equilibrium with no infection
  if(NoInfeq(&ytmp,p)>0) {
    yeq[*nbeq] = ytmp;
    xeq[*nbeq] = 0.0;
    if(p.ksi>0) {
      zeq[*nbeq] = p.eta*yeq[*nbeq]/p.ksi;    
    } else {
      //damages are not repaired
      if(p.eta<=0.0) {
      //pathogens have been cleared and defenses cause no damage
      //z final value is finite but depends on initial x value
      //z is therefore set to nan
	zeq[*nbeq] =  NAN;
      } else {
      //pathogens have been cleared but defenses cause damages
      //z goes to infinity
	zeq[*nbeq] = INFINITY;
      }
    }
    (*nbeq)++;
  }

  //equilibrium with xeq>0  
  if(p.omega >0 && p.psi>0  && p.ksi<=0) {
    xeq[*nbeq]=1.0;
    yeq[*nbeq]=0.0;
    zeq[*nbeq]=INFINITY;
    (*nbeq)++;
  } else {
    int status;
    
    //const gsl_root_fdfsolver_type *T = gsl_root_fdfsolver_newton;
    //gsl_root_fdfsolver *s = gsl_root_fdfsolver_alloc (T);
    
    //find roots that may lie in segments [x;x+dx]
    gsl_function F_H;
    F_H.function = &Hlogx;
    F_H.params = &p;

    gsl_function F_dH;
    F_dH.function = &dHlogx;
    F_dH.params = &p;

    const gsl_root_fsolver_type *Tr = gsl_root_fsolver_brent;
    gsl_root_fsolver *sr = gsl_root_fsolver_alloc (Tr);

    int    iter = 0;
    int    nb_lx  = 5000;
    double lx_min=-64.0;
    double lx_max=0.0;
    double lx_incr_old = (lx_max-lx_min)/(double)(nb_lx-1);
    double lx_incr = lx_incr_old;
    double lxtmp;
    
    for(double lx_max_tmp=lx_min+lx_incr; lx_max_tmp<=lx_max & *nbeq<nbeq_max; lx_max_tmp+=lx_incr) {
      
      double lx_min_tmp = lx_max_tmp-lx_incr;
      double h_lx_min_tmp = Hlogx(lx_min_tmp,&p);
      double h_lx_max_tmp = Hlogx(lx_max_tmp,&p);
      lx_incr = lx_incr_old;
	      
      if(h_lx_min_tmp*h_lx_max_tmp > 0) {
	//H has the same sign on both lx_min_tmp and lx_max_tmp
	//check that the sign of it's derivative does not change
	double dh_lx_min_tmp = dHlogx(lx_min_tmp,&p);
	double dh_lx_max_tmp = dHlogx(lx_max_tmp,&p);

	if(dh_lx_min_tmp*dh_lx_max_tmp < 0) {
	  //the sign of the derivative has changed!
	  // two roots may lie between lx_min_tmp and lx_max_tmp
	  gsl_root_fsolver_set (sr, &F_dH, lx_min_tmp, lx_max_tmp);
	  iter=0;
	  do
	    {
	      iter++;
	      status = gsl_root_fsolver_iterate (sr);
	      lxtmp = gsl_root_fsolver_root (sr);
	      double lx_lo = gsl_root_fsolver_x_lower (sr);
	      double lx_hi = gsl_root_fsolver_x_upper (sr);
	      status = gsl_root_test_interval (lx_lo,lx_hi,0,1e-6);
	    }
	  while (status == GSL_CONTINUE && iter < maxiter);
  
	  if (status == GSL_SUCCESS) {
	    lx_incr = lx_max_tmp - lxtmp;
	    lx_max_tmp = lxtmp;	    
	  }
	}
      }
      
      if(h_lx_min_tmp*h_lx_max_tmp < 0) {
	//there's a root somewhere inbetween lx_min_tmp and lx_max_tmp !
	gsl_root_fsolver_set (sr, &F_H, lx_min_tmp, lx_max_tmp);
	iter=0;
	do
	  {
	    iter++;
	    status = gsl_root_fsolver_iterate (sr);
	    lxtmp = gsl_root_fsolver_root (sr);
	    double lx_lo = gsl_root_fsolver_x_lower (sr);
	    double lx_hi = gsl_root_fsolver_x_upper (sr);
	    status = gsl_root_test_interval (lx_lo,lx_hi,0,1e-6);
	  }
	while (status == GSL_CONTINUE && iter < maxiter);
  
	if (status == GSL_SUCCESS) {
	  *nbeq = setYZ(*nbeq,exp(lxtmp),xeq,yeq,zeq,p);
	}

	if((*nbeq) >= nbeq_max) {
#ifdef RMOD
	  warning("WARNING : cannot compute more than ",nbeq_max," equilibria!\n");
#else	
	  printf("WARNING : cannot compute more than %d equilibria!\n",nbeq_max);
#endif
	}	
      }
    }
    if(iter>=maxiter) {
#ifdef RMOD
	  warning("WARNING : number of iterations exceeded ",maxiter," while computing equilibria!\n");
#else	
	  printf("WARNING : number of iterations exceeded %d while computing equilibria!\n",nbeq_max);
#endif
    }
    gsl_root_fsolver_free (sr);
    gsl_set_error_handler (old_handler);
  }

  //store equilibria in an array
  double** eq = malloc(6*sizeof(double*)); //last 2 columns are of real and imaginary parts of dominant eigenvalues
  
  if(*nbeq<nbeq_max) {
    eq[0] = realloc(xeq,(*nbeq)*sizeof(double));
    eq[1] = realloc(yeq,(*nbeq)*sizeof(double));
    eq[2] = realloc(zeq,(*nbeq)*sizeof(double));
  } else {
    eq[0] = xeq;
    eq[1] = yeq;
    eq[2] = zeq;
  }
  eq[3] = malloc((*nbeq)*sizeof(double));
  eq[4] = malloc((*nbeq)*sizeof(double));
  eq[5] = malloc((*nbeq)*sizeof(double));
		 
  for(int i=2;i<(*nbeq);i++) {
    if(dH(eq[0][i-1],&p)*dH(eq[0][i],&p)>0) {
      printf("WARNING one equilibrium is probably missing between x=%g and x=%g\n",eq[0][i-1],eq[0][i]);
    }
  }


  //compute eigenvalues for each equilibrium
  
  //in the particular case where ksi=0, the complete system has no equilibrium point as z has no equilibrium value
  //if psi is zero, stability is then computed on the reduced system, with z excluded
  //if psi>0 and ksi=0 stability cannot be computed.

  //the stability of homeostatic state (equilibrium with no bacteria) may be impossible to compute from the
  //jacobian matrix because dFdx can be infinite when x=0 if u or v are less than one. A solution to this
  //issue is to use (1/x) dx/dt as the dominant eigenvalue.

  if(*nbeq>0) {
    eq[3][0] = 1-eq[1][0]*p.delta;
  }
  
  int size;
  if(p.ksi<=0) {
    if(p.psi<=0) {
      size=2;
    } else {
      size=-1;
    }
  } else {
    size=3;
  }
  if(size>0) {
    gsl_vector_complex *eval = gsl_vector_complex_alloc (size);		 
    gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc (size);

    //option of the eigenvalue routines
    // If compute_t is set to 1, the full Schur form T will be computed by gsl_eigen_nonsymm()
    // If balance is set to 1, a balancing transformation is applied to the matrix prior to computing eigenvalues.
    // The second option appears to be usefull when x approaches zero.
    // https://www.gnu.org/software/gsl/doc/html/eigen.html
    gsl_eigen_nonsymm_params(0,1,w);
  
    
    double* y = malloc(size*sizeof(double));
    double *dfdy = malloc(size*size*sizeof(double));

    for(int i=1;i<*nbeq;i++) {    //stability of equilibrium with no bacteria is not computed here!
      for(int j=0;j<size;j++) y[j] = eq[j][i];
      if(compute_jacob_dxyz_dt(y,dfdy,size==2?0:1,p)==0) {
	//for(int j=0;j<9;j++) printf("%d -- %g\n",j,dfdy[j]);    
	gsl_matrix_view m = gsl_matrix_view_array (dfdy,size,size);
    
	gsl_eigen_nonsymm (&m.matrix, eval, w); //eigenvalues are stored in eval
    
	gsl_vector_view re_ev = gsl_vector_complex_real(eval);
	gsl_vector_view im_ev = gsl_vector_complex_imag(eval);

	//look for eigenvalue with the highest real part
	//and test if some eigen values have non-zero imaginary part
	double max_re_ev = gsl_vector_get(&(re_ev.vector),0);
	double max_im_ev = gsl_vector_get(&(im_ev.vector),0);
	double max_im_ev_glob = max_im_ev;
	for(int j=1;j<size;j++) {
	  double tmp_re = gsl_vector_get(&(re_ev.vector),j);
	  double tmp_im = gsl_vector_get(&(im_ev.vector),j);
	  if(tmp_re>max_re_ev) {
	    max_re_ev = tmp_re;
	    max_im_ev = tmp_im;
	  }
	  if(fabs(tmp_im) > fabs(max_im_ev_glob)) max_im_ev_glob = tmp_im;	  
	}
	eq[3][i] = max_re_ev;
	eq[4][i] = max_im_ev;
	eq[5][i] = max_im_ev_glob;
	//printf("%d -- x=%f y=%f z=%f -- ev = %f + %f.i\n",i,y[0],y[1],y[2],max_re_ev,max_im_ev);    
      } else {
	eq[3][i] = 0.0;
	eq[4][i] = 0.0;
	eq[5][i] = 0.0;	
      }
    }
    free(y);
    free(dfdy);
    gsl_eigen_nonsymm_free(w);
    gsl_vector_complex_free(eval);
  } else {
    for(int i=0;i<*nbeq;i++) {
      eq[3][i]=-1;
      eq[4][i]=-1;
      eq[5][i]=-1;
    }
  }
      
  return eq; 
}

//one infection: returns state at each time t
int dyninf_complete(int nb,double* x, double* y, double* z,double* t,Parameter* p){
  gsl_odeiv2_system sys = { dxyz_dt, jacob_dxyz_dt, 3, p };
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_step * s  = gsl_odeiv2_step_alloc (T, 3);
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-8, 0.0);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (3);  
  
  double* var = (double*) malloc(3*sizeof(double));
  var[0] = x[0];
  var[1] = y[0];
  var[2] = z[0];
  double time = t[0];
  double h = 1e-6;
  int direction = 1;
  if(nb>1 && t[1]<t[0]) {
    direction=-1;
    h*=-1;
  }
  
  int i;
  for(i=1;i<nb;i++) {
    int status=GSL_SUCCESS;
    while (direction*time < direction*t[i]) {
      status = gsl_odeiv2_evolve_apply (e, c, s,
					&sys,
					&time, t[i],
					&h, var);	  
      if ( status!=GSL_SUCCESS ) {
#ifdef RMOD
	warning("GSL failed integrating the system!");
#else
	printf("GSL failed integrating the system!");
#endif
	free(var);
	return(i);
      }
    }
    x[i] = var[0];
    y[i] = var[1];
    z[i] = var[2];
  }
  
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  free(var);
  
  return(i);
}

//one infection: 
int dyninf(double* var,double tmax, int stopWhenDead,int logx,Parameter* p){
  // integration will use runge-kutta order 8 method
  gsl_odeiv2_system sys = { logx==1?dlogxyz_dt:dxyz_dt, logx==1?jacob_dlogxyz_dt:jacob_dxyz_dt, 3, p };
  if(logx==1) var[0] = log(var[0]);
  
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_step * s  = gsl_odeiv2_step_alloc (T, 3);
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-6, 0.0);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (3);  

  int    eventHasOccurred = 0;
  double t = 0;
  double h = 1e-6;

  double  t_old;
  double* var_old = (double*) malloc(3*sizeof(double));
  
  int status=GSL_SUCCESS;
  while (t < tmax) {
    t_old=t;
    for(int i=0;i<3;i++) var_old[i] = var[i];
    
    status = gsl_odeiv2_evolve_apply (e, c, s,
				      &sys,
				      &t, tmax,
				      &h, var);	  
    if ( status!=GSL_SUCCESS ) {
#ifdef RMOD
      warning("GSL failed integrating the system!");
#else
      printf("GSL failed integrating the system!");
#endif
      eventHasOccurred = -1;
      break;
    }

    if(stopWhenDead==1 && var[2]>=p->zcrit) {
      eventHasOccurred = 1;
      for(int i=0;i<3;i++) var[i] = var_old[i];
      t=t_old;
      break; //the host is dead !
    }
  }
  var[3] = t;
  if(logx==1) var[0] = exp(var[0]);
  
  free(var_old);
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);

  return(eventHasOccurred);
}

double xdecrease(double x,double y,double z,Parameter p) {
  double dx =  dxdt(x,y,z,p);
  return(dx);
}

double EVENT_xstart,EVENT_ystart,EVENT_zstart,EVENT_mindist;
double xlocalmin(double x,double y,double z,Parameter p) {
  double xst = x-EVENT_xstart;
  double yst = y-EVENT_ystart;
  double zst = z-EVENT_zstart;  
  double d  = sqrt(xst*xst+yst*yst+zst*zst); //distance to starting point
  if(d<=EVENT_mindist) return -1.0; //system hasn't moved from starting point
  double d2x =  d2xdt2(x,y,z,p);
  if(d2x<0) return -1.0; //secondary derivative is negative: this cannot be a local min!  
  double dx  =  dxdt(x,y,z,p);
  return dx;
}

double xlocalmax(double x,double y,double z,Parameter p) {
  double xst = x-EVENT_xstart;
  double yst = y-EVENT_ystart;
  double zst = z-EVENT_zstart;  
  double d  = sqrt(xst*xst+yst*yst+zst*zst); //distance to starting point
  if(d<=EVENT_mindist) return 1.0; //system hasn't moved from starting point
  double d2x =  d2xdt2(x,y,z,p);
  if(d2x>0) return 1.0; //secondary derivative is positive: this cannot be a local max!
  double dx  =  dxdt(x,y,z,p);
  return dx;
}

double EVENT_x;
double xcritreached(double x,double y,double z,Parameter p) {
  return x-EVENT_x;
}

double logxcritreached(double x,double y,double z,Parameter p) {
  return log(x/EVENT_x);
}

double ydecrease(double x,double y,double z,Parameter p) {
  double dy = dydt(x,y,z,p);
  return(dy);
}

double yinflection(double x,double y,double z,Parameter p) {
  double dy = d2ydt2(x,y,z,p);
  //printf("%g %g %g %g\n",x,y,z,dy);
  return(dy);
}

double EVENT_xe,EVENT_ye,EVENT_ze,EVENT_dist;
double closetoeq(double x,double y,double z,Parameter p) {
  //derivative of the distance to a point (typically an unstable equilibrium point)
  double xst = x-EVENT_xe;
  double yst = y-EVENT_ye;
  double zst = z-EVENT_ze;  
  double d  = sqrt(xst*xst+yst*yst+zst*zst);
  double dd = (dxdt(x,y,z,p)*xst + dydt(x,y,z,p)*yst + dzdt(x,y,z,p)*zst)/d;
  return dd;
}

double EVENT_lxe,EVENT_lye,EVENT_lze,EVENT_ldist;
double distLtoeq(double x,double y,double z,Parameter p) {
  //derivative of the distance to a point (typically an unstable equilibrium point)
  double xst = log(x)-EVENT_lxe;
  double yst = log(y)-EVENT_lye;
  double zst = log(z)-EVENT_lze;  
  double d  = sqrt(xst*xst+yst*yst+zst*zst);
  return d-EVENT_ldist;
}

double** EVENT_eq;
int EVENT_nbeq;
double distToAnyEq(double x,double y,double z,Parameter p) {
  double d  = -1;
  int pos = 0;
  //distance to the nearest stable equilibrium
  for(int i=0;i<EVENT_nbeq;i++) {
    if(EVENT_eq[3][i]<0) { // the equilibrium is stable
      double dtmp = 0.0;
      dtmp += pow((x)-(EVENT_eq[0][i]),2.0);
      dtmp += pow((y)-(EVENT_eq[1][i]),2.0);
      dtmp += pow((z)-(EVENT_eq[2][i]),2.0);
      dtmp = sqrt(dtmp);
      if(d<0 || dtmp<d) {
		d = dtmp;
		pos = i;
      }
    }
  }
  //printf("%d/%d %g %g %g\n",pos+1,EVENT_nbeq,d-EVENT_ldist,x,EVENT_eq[0][pos]);
  return d-EVENT_ldist;
}

double zcritreached(double x,double y,double z,Parameter p) {
  return z-p.zcrit;
}

int dyninfWEvent(double *var,double tmax,double* eps,int maxiter,double (*event)(double,double,double,Parameter), int logx,Parameter* p){
  // integration will use runge-kutta order 8 method
  gsl_odeiv2_system sys = { logx==1?dlogxyz_dt:dxyz_dt, logx==1?jacob_dlogxyz_dt:jacob_dxyz_dt, 3, p };
  if(logx==1) var[0] = log(var[0]);
  //gsl_odeiv2_system sys = { dxyz_dt, jacob_dxyz_dt, 3, p };
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_step * s  = gsl_odeiv2_step_alloc (T, 3);
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-6, 0.0);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (3);

  int status=GSL_SUCCESS;
  int eventHasOccurred = 0;
  double tmax2 = tmax;
  double h = 1e-6;    
  double d0;
  if(logx==1) {
    d0 = (*event)(exp(var[0]),var[1],var[2],*p);
  } else {
    d0 = (*event)(var[0],var[1],var[2],*p);
  }
  double d=d0;
  double t = var[3];
  
  double* var_prev = (double*) malloc(4*sizeof(double));
  for(int i=0;i<3;i++) var_prev[i] = var[i];
  double h_prev = h;
  double t_prev = t;
  double d_prev = d;


  int    iter=0;
  while (t < tmax && iter<maxiter) {
    status = gsl_odeiv2_evolve_apply (e, c, s,
				      &sys,
				      &t, tmax2,
				      &h, var);
      if ( status!=GSL_SUCCESS ) {
#ifdef RMOD
      warning("GSL failed integrating the system!");
#else
      printf("\nGSL failed integrating the system!\n");
#endif
      eventHasOccurred = -1;
      break;
    }
    if(logx==1) {
      d = (*event)(exp(var[0]),var[1],var[2],*p);
    } else {
      d = (*event)(var[0],var[1],var[2],*p);
    }

    //printf("%d -- h=%g t=%g -- x=%g y=%g z=%g d=%g (%g)\n",iter,h,t,var[0],var[1],var[2],d,*eps);
    if(fabs(d)<=*eps & iter>5) {
      //root has been reached!
      //printf("t=%g : event! (%g)\n",t,d);
      eventHasOccurred = 1;
      break;
    }
    
    if(d0*d<0)  {
      //printf("overshoot!\n");
      //overshoot! Go one step backward and change tmax
      tmax2 = t;      
      for(int i=0;i<3;i++) var[i] = var_prev[i];
      t = t_prev;
      d = d_prev;
      //reduce h
      //printf("%g %g ",h_prev,h);
      h = h_prev/2.0;      
      //printf("%g\n",h);
    } else {
      //compute step size according to distance to root      
      double d3 = d+(d-d_prev)*h/h_prev; //prediction of next distance if step size h is applied
      if(d0*d3<0) {//predicted overshoot: step size must be reduced
        double hd = 0.9*(d*h_prev)/(d_prev-d);
	//printf("h=%g -> h=%g\n",h,hd);
	    h=hd;
      }
      //update previous var, t and h states
      for(int i=0;i<3;i++) var_prev[i] = var[i];
      t_prev = t;
      h_prev = h;
      d_prev = d;
    }
    iter++;
  }
  var[3] = t;
  if(logx==1) var[0] = exp(var[0]);

    //double d = (*event)(var[0],var[1],var[2],*p);
  if(fabs(d)>*eps) {
#ifdef RMOD
    char* msg = (char*) malloc(100*sizeof(char));
    sprintf(msg,"dyninfWEvent failed to locate event with required precision! (t=%f dist=%f)",t,fabs(d));
    if(iter>=maxiter) sprintf(msg,"\nNumber of interations exceeded %d",maxiter);
    warning(msg);
    free(msg);
#else
      printf("WARNING dyninfWEvent failed to locate event with required precision! (t=%g dist=%g)\n",t,fabs(d));
      printf("Number of interations exceeded %d\n",maxiter);
      printf("alpha=%f gamma=%f -- ",p->alpha,p->gamma);
      for(int i=0;i<3;i++) printf("%g ",var[i]);
      printf("\n");
#endif      
  }
  *eps = d;

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);

  free(var_prev);
  if(iter>=maxiter) return(-5);
  return(eventHasOccurred);
}

/* */
int tsurv(double* var, Parameter* p){
  
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_system sys = { dxyt_dz, jacob_dxyt_dz, 3, p };
  gsl_odeiv2_step * s  = gsl_odeiv2_step_alloc (T, 3);  
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (3);  
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (1e-6, 0.0);
  
  double* var2 = (double*) malloc(3*sizeof(double));
  var2[0] = var[0]; //x;
  var2[1] = var[1]; //y
  var2[2] = var[3]; //t
  double z  = var[2]; //z
  double h = 1e-6;
  if(var[2]>p->zcrit) {
    //backward integration
    h*=-1.0;
  }

  int status=GSL_SUCCESS;
  int cmpt=0,cmptmax=10000;
  //printf("%d -- (%g) -- t=%g x=%g y=%g z=%g\n",cmpt,h,var[2],var[0],var[1],z);
  while (fabs(z-p->zcrit)>1e-6 && cmpt<cmptmax) {
    
    status = gsl_odeiv2_evolve_apply (e, c, s,
				      &sys,
				      &z, p->zcrit,
    				      &h, var2);	      
    if ( status!=GSL_SUCCESS ) {
	#ifdef RMOD
		warning("GSL failed integrating the system!");
	#else
		printf("GSL failed integrating the system!\n");
	#endif
      break;
    }    
    cmpt++;
    //printf("%d -- (%g) -- t=%g x=%g y=%g z=%g\n",cmpt,h,var[2],var[0],var[1],z);
  }
  
  var[0] = var2[0];
  var[1] = var2[1];
  var[2] = z;
  var[3] = var2[2];
    
  free(var2);
  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  
  return(status);
}

double dist_to_equilibrium(double *v,double *e) {
  double d = 0.0;
  for(int i=0;i<3;i++)  d += pow(v[i]-e[i],2.0);
  d = sqrt(d);
  return(d);
}

int threshold_x(int nb,double* x,double* lo_x, double *hi_x,double* y,double* z,double tmax,double eps,double max_dist,Parameter p) {
  //computing equilibria
  int maxiter = 100;
  int nbeq = -1;
  double** eq = equilibrium(&nbeq,p);
  if(nbeq<3) return(-1); //the system is not bistable
  
  //positions of the stable equilibria
  int pos_1,pos_2;
  for(pos_1=0;pos_1<nbeq;pos_1++) {
    if(eq[3][pos_1]<0) break;
  }
  if(pos_1>=nbeq-1) return(-1); //the system is not bistable
  for(pos_2=pos_1+1;pos_2<nbeq;pos_2++) {
    if(eq[3][pos_2]<0) break;   
  }
  if(pos_2>nbeq || eq[4][pos_2]>0) return(-1); //the system is not bistable
  
  double* e1 = (double*) malloc(3*sizeof(double));
  double* e2 = (double*) malloc(3*sizeof(double));
  for(int i=0;i<3;i++) {
    e1[i] = eq[i][pos_1];
    e2[i] = eq[i][pos_2];
  }

  if(dist_to_equilibrium(e1,e2)<=2*max_dist) {
    //the system is bistable, but the equilibrium are too close
    //for threshold x to be found
    return(-1);
  }

  int cmpt = 0;
  for(int i=0;i<nb;i++) {
    if(y[i]<0) y[i] = eq[1][0]; //if y is negative, the clearance equilibrium is used as starting point
    if(z[i]<0) z[i] = eq[2][0]; //idem for z   

    //starting points
    double lx1=log(1e-24);
    double lx2=log(1-1e-24);

    double* v1 = malloc(4*sizeof(double));    
    v1[0] = exp(lx1);
    v1[1]=y[i];
    v1[2]=z[i];
    double* v2 = malloc(4*sizeof(double));
    v2[0] = exp(lx2);
    v2[1]=y[i];
    v2[2]=z[i];
    
    dyninf(v1,tmax,0,1,&p); //logx !
    dyninf(v2,tmax,0,1,&p);
    
    //printf("x1=%g -> x=%g / x=%g -- %g\n",x1,v1[0],e1[0],dist_to_equilibrium(v1,e1));
    //printf("x2=%g -> x=%g / x=%g -- %g\n",x2,v2[0],e2[0],dist_to_equilibrium(v2,e2));
    
    if(dist_to_equilibrium(v1,e1)<max_dist && dist_to_equilibrium(v2,e2)<max_dist) {
      int test = 1;
      int iter = 0;
      while(fabs(lx2-lx1)>eps && test==1 && iter<maxiter) {	
	double lx3 = (lx1+lx2)/2.0;
	double* v3 = malloc(4*sizeof(double));    
	v3[0] = exp(lx3);
	v3[1]=y[i];
	v3[2]=z[i];	
	dyninf(v3,tmax,0,1,&p); //logx
	
	double de1 = dist_to_equilibrium(v3,e1);
	double de2 = dist_to_equilibrium(v3,e2);
	//printf("%d -- %g < %g < %g -- %g %g\n",iter,lx1,lx3,lx2,de1,de2);
	
	test = 0;
	if(de1<max_dist) {
	  //x3 goes to e1
	  free(v1);
	  v1 = v3;
	  lx1 = lx3;
	  test = 1;
	}
	if(de2<max_dist) {
	  //x3 goes to e2
	  free(v2);
	  v2 = v3;
	  lx2 = lx3;
	  test = 1;
	}
	if(test==0) {
	  free(v3);
	  printf("Lowest distance to equilibria (%g) is above max_dist = %g! Try increasing tmax.\n",de1<de2?de1:de2,max_dist);
	}
	iter++;
      }//end of while
      if(test==1) {
	x[i] = exp((lx1+lx2)/2.0);
	lo_x[i] = exp(lx1);
	hi_x[i] = exp(lx2);	
	cmpt++;
      } else {
	x[i] = NAN;
	lo_x[i] = exp(lx1);
	hi_x[i] = exp(lx2);	
      }
    } else {
      if(dist_to_equilibrium(v1,e2)<max_dist) {
	x[i]=1.0;
	lo_x[i] = NAN;
	hi_x[i] = NAN;	
      } else {
	x[i] = NAN;
	lo_x[i] = NAN;
	hi_x[i] = NAN;	
      }
    }
    free(v1);
    free(v2);
  }
  
  for(int i=0;i<6;i++) free(eq[i]);
  free(eq);
  free(e1);
  free(e2);

  return(cmpt);
}

int threshold_y(int nb,double* x,double* y,double* lo, double *hi,double* z,double tmax,double eps,double max_dist,Parameter p) {
  //computing equilibria
  int maxiter = 100;
  int nbeq = -1;
  double** eq = equilibrium(&nbeq,p);
  if(nbeq<3) return(-1); //the system is not bistable
  
  //positions of the stable equilibria
  int pos_1,pos_2;
  for(pos_1=0;pos_1<nbeq;pos_1++) {
    if(eq[3][pos_1]<0) break;
  }
  if(pos_1>=nbeq-1) return(-1); //the system is not bistable
  for(pos_2=pos_1+1;pos_2<nbeq;pos_2++) {
    if(eq[3][pos_2]<0) break;   
  }
  if(pos_2>nbeq || eq[4][pos_2]>0) return(-1); //the system is not bistable
  
  double* e1 = (double*) malloc(3*sizeof(double));
  double* e2 = (double*) malloc(3*sizeof(double));
  for(int i=0;i<3;i++) {
    e1[i] = eq[i][pos_1];
    e2[i] = eq[i][pos_2];
  }

  if(dist_to_equilibrium(e1,e2)<=2*max_dist) {
    //the system is bistable, but the equilibrium are too close
    //for threshold x to be found
    return(-1);
  }

  int cmpt = 0;
  for(int i=0;i<nb;i++) {
    //starting points
    double ly1=log(1e-24);
    double ly2=log(2.0/p.delta);

    double* v1 = malloc(4*sizeof(double));    
    v1[0]=x[i];
    v1[1]=exp(ly1);
    v1[2]=z[i];
    double* v2 = malloc(4*sizeof(double));
    v2[0]=x[i];
    v2[1]=exp(ly2);
    v2[2]=z[i];
    
    dyninf(v1,tmax,0,1,&p); //logx !
    dyninf(v2,tmax,0,1,&p);
    double de1 = dist_to_equilibrium(v1,e2);
    double de2 = dist_to_equilibrium(v2,e1);
    if(de1<max_dist && de2 <max_dist) {
      //starting points must be reversed
      double tmp = ly2;
      ly2 = ly1;
      ly1 = tmp;
      tmp = de2;
      de2 = de1;
      de1 = tmp;
    } else {
      de1 = dist_to_equilibrium(v1,e1);
      de2 = dist_to_equilibrium(v2,e2);      
    }
    
    printf("y1=%g -> (%f,%f,%f) / (%f,%f,%f) %f\n",exp(ly1),v1[0],v1[1],v1[2],e1[0],e1[1],e1[2],de1);
    printf("y2=%g -> (%f,%f,%f) / (%f,%f,%f) %f\n",exp(ly2),v2[0],v2[1],v2[2],e2[0],e2[1],e2[2],de2);
    
    
    if(de1<max_dist && de2<max_dist) {
      int test = 1;
      int iter = 0;
      while(fabs(ly2-ly1)>eps && test==1 && iter<maxiter) {	
	double ly3 = (ly1+ly2)/2.0;
	double* v3 = malloc(4*sizeof(double));    
	v3[0]=x[i];
	v3[1]=exp(ly3);
	v3[2]=z[i];	
	dyninf(v3,tmax,0,1,&p); //logx
	
	de1 = dist_to_equilibrium(v3,e1);
	de2 = dist_to_equilibrium(v3,e2);
	printf("y3=%g -> (%f,%f,%f) / %f / %f\n",exp(ly3),v3[0],v3[1],v3[2],de1,de2);

	
	test = 0;
	if(de1<max_dist) {
	  //x3 goes to e1
	  free(v1);
	  v1 = v3;
	  ly1 = ly3;
	  test = 1;
	}
	if(de2<max_dist) {
	  //x3 goes to e2
	  free(v2);
	  v2 = v3;
	  ly2 = ly3;
	  test = 1;
	}
	if(test==0) {
	  free(v3);
	  printf("Lowest distance to equilibria (%g) is above max_dist = %g! Try increasing tmax.\n",de1<de2?de1:de2,max_dist);
	}
	iter++;
      }//end of while
      if(test==1) {
	y[i] = exp((ly1+ly2)/2.0);
	lo[i] = exp(ly1);
	hi[i] = exp(ly2);	
	cmpt++;
      } else {
	y[i] = NAN;
	lo[i] = exp(ly1);
	hi[i] = exp(ly2);	
      }
    } else {
      if(dist_to_equilibrium(v1,e2)<max_dist) {
	y[i]= NAN;
	lo[i] = NAN;
	hi[i] = NAN;	
      } else {
	y[i] = NAN;
	lo[i] = NAN;
	hi[i] = NAN;	
      }
    }
    free(v1);
    free(v2);
  }
  
  for(int i=0;i<6;i++) free(eq[i]);
  free(eq);
  free(e1);
  free(e2);

  return(cmpt);
}

int threshold_z(int nb,double* x,double* y,double* z,double* lo, double *hi,double tmax,double eps,double max_dist,Parameter p) {
  //computing equilibria
  int maxiter = 100;
  int nbeq = -1;
  double** eq = equilibrium(&nbeq,p);
  if(nbeq<3) return(-1); //the system is not bistable
  
  //positions of the stable equilibria
  int pos_1,pos_2;
  for(pos_1=0;pos_1<nbeq;pos_1++) {
    if(eq[3][pos_1]<0) break;
  }
  if(pos_1>=nbeq-1) return(-1); //the system is not bistable
  for(pos_2=pos_1+1;pos_2<nbeq;pos_2++) {
    if(eq[3][pos_2]<0) break;   
  }
  if(pos_2>nbeq || eq[4][pos_2]>0) return(-1); //the system is not bistable
  
  double* e1 = (double*) malloc(3*sizeof(double));
  double* e2 = (double*) malloc(3*sizeof(double));  
  for(int i=0;i<3;i++) {
    e1[i] = eq[i][pos_1];
    e2[i] = eq[i][pos_2];
  }
  if(dist_to_equilibrium(e1,e2)<=2*max_dist) {
    //the system is bistable, but the equilibrium are too close
    //for threshold x to be found
    return(-1);
  }

  int cmpt = 0;
  for(int i=0;i<nb;i++) {
    //starting points
    double lz1=log(1e-24);
    double lz2=log(2.0*p.zcrit);

    double* v1 = malloc(4*sizeof(double));    
    v1[0]=x[i];
    v1[1]=y[i];
    v1[2]=exp(lz1);
    double* v2 = malloc(4*sizeof(double));
    v2[0]=x[i];
    v2[1]=y[i];
    v2[2]=exp(lz2);
    
    dyninf(v1,tmax,0,1,&p); //logx !
    dyninf(v2,tmax,0,1,&p);
    //test if conditions are in the right order
    double de1 = dist_to_equilibrium(v1,e2);
    double de2 = dist_to_equilibrium(v2,e1);
    if(de1<max_dist && de2 <max_dist) {
      //starting points must be reversed
      double tmp = lz2;
      lz2 = lz1;
      lz1 = tmp;
      tmp = de2;
      de2 = de1;
      de1 = tmp;
    } else {
      de1 = dist_to_equilibrium(v1,e1);
      de2 = dist_to_equilibrium(v2,e2);      
    }
    
    /* printf("z1=%g -> (%f,%f,%f) / (%f,%f,%f) %f\n",exp(lz1),v1[0],v1[1],v1[2],e1[0],e1[1],e1[2],de1); */
    /* printf("z2=%g -> (%f,%f,%f) / (%f,%f,%f) %f\n",exp(lz2),v2[0],v2[1],v2[2],e2[0],e2[1],e2[2],de2); */
    
    
    if(de1<max_dist && de2<max_dist) {
      int test = 1;
      int iter = 0;
      while(fabs(lz2-lz1)>eps && test==1 && iter<maxiter) {	
	double lz3 = (lz1+lz2)/2.0;
	double* v3 = malloc(4*sizeof(double));    
	v3[0]=x[i];
	v3[1]=y[i];
	v3[2]= exp(lz3);	
	dyninf(v3,tmax,0,1,&p); //logx
	
	de1 = dist_to_equilibrium(v3,e1);
	de2 = dist_to_equilibrium(v3,e2);
	/* printf("z3=%g -> (%f,%f,%f) / %f / %f\n",exp(lz3),v3[0],v3[1],v3[2],de1,de2); */
	
	test = 0;
	if(de1<max_dist) {
	  //x3 goes to e1
	  free(v1);
	  v1 = v3;
	  lz1 = lz3;
	  test = 1;
	}
	if(de2<max_dist) {
	  //x3 goes to e2
	  free(v2);
	  v2 = v3;
	  lz2 = lz3;
	  test = 1;
	}
	if(test==0) {
	  free(v3);
	  printf("Lowest distance to equilibria (%g) is above max_dist = %g! Try increasing tmax.\n",de1<de2?de1:de2,max_dist);
	}
	iter++;
      }//end of while
      if(test==1) {
	z[i] = exp((lz1+lz2)/2.0);
	lo[i] = exp(lz1);
	hi[i] = exp(lz2);	
	cmpt++;
      } else {
	z[i] = NAN;
	lo[i] = exp(lz1);
	hi[i] = exp(lz2);	
      }
    } else {
      if(de2<max_dist) {
	z[i]=p.zcrit;
	lo[i] = NAN;
	hi[i] = NAN;	
      } else {
	z[i] = NAN;
	lo[i] = NAN;
	hi[i] = NAN;	
      }
    }
    free(v1);
    free(v2);
  }
  
  for(int i=0;i<6;i++) free(eq[i]);
  free(eq);
  free(e1);
  free(e2);

  return(cmpt);
}


//gamma value above which the host is protected against infection by a dose x0
//meaningless if the system is not bistable
int f_infection(double x0,double dlogy0,double dz0,double tmax,double maxdist,Parameter p) {
  int nbeq;
  double** eq = equilibrium(&nbeq,p);

  //printf("gamma = %g, %d equilibria\n",p.gamma,nbeq);
  if(nbeq<3) printf("WARNING f_infection should not be called on systems with only %d equilibria!\n",nbeq);
  
  //set starting conditions
  double* var = (double*) malloc(4*sizeof(double));
  var[0] = x0;
  var[1] = eq[1][0]*exp(dlogy0);
  var[2] = eq[2][0]+dz0*p.zcrit;
  var[3] = 0.0;

  //dynamics
  //dyninf(var,tmax,0,1,&p);
  EVENT_eq   = eq;
  EVENT_nbeq = nbeq;
  EVENT_ldist = maxdist;
  int maxiter = 10000;
  double eps = 1e-6; //precision of the estimation of the time at which event has occurred (only rough estimates are needed here!)
  double epsold = eps;
  //printf("alpha=%f gamma=%f (%f %f %f) (%f %f %f)\n",p.alpha,p.gamma,var[0],var[1],var[2],eq[0][nbeq-1],eq[1][nbeq-1],eq[2][nbeq-1]);
  int test = dyninfWEvent(var,tmax,&eps,maxiter,&distToAnyEq,1,&p);

  //computed distance to equilibria
  double dist = 0.0;
  for(int i=0;i<3;i++) dist += pow((var[i])-(eq[i][nbeq-1]),2.0);
  dist = sqrt(dist);

  //printf("alpha=%f gamma=%f (%f %f %f) (%f %f %f) %f\n",p.alpha,p.gamma,var[0],var[1],var[2],eq[0][nbeq-1],eq[1][nbeq-1],eq[2][nbeq-1],dist);

  for(int i=0;i<6;i++) free(eq[i]);  
  free(eq);
  free(var);
  
  if(fabs(dist-maxdist)<epsold) return 1;
  return 0;
}

int gamma_infection(double x0,double dy0,double dz0,double* g_crit, double tmax,double maxdist,double eps,int maxiter,double* estim,Parameter p) {

  Parameter p2;
  copyParameter(p,&p2); 
    
  int g_crit_computed = 0;
  if(g_crit == NULL) {
    //gamma values that bound the bistable region need to be computed
    g_crit = (double*) malloc(10*sizeof(double));
    int    *curv   = (int*) malloc(10*sizeof(int));
    int nb = gamma_crit (10,g_crit,curv,p);

    if(nb>2) {
      //should never happen!
      printf("gamma_infection: More than two critical values of gamma have been found!\n");
      freeParameter(&p2);
      return -2;
    }
    if(nb < 1){
     freeParameter(&p2);
     return -1;
     } //the system is not bistable
    if(nb<2) {
      double g_hi=-1.0,g_lo=-1.0;
      for(int i=0;i<nb;i++) {
	if(curv[i]==1)  g_hi=g_crit[i];
	if(curv[i]==-1) g_lo=g_crit[i];
	if(curv[i]==0) {
	  printf("WARNING! Cannot tell whether gamma_crit is a maximum or a minimum!\n");
	}
      }
      if(g_lo>=0) {
	g_crit[0] = g_lo;
      } else {
	if(nb>1) {
	  g_crit[0] = 0.0;
	} else {
	  g_crit[0] = gamma_clearance(p);
	}
      }
      if(g_hi>=0) {
	g_crit[1] = g_hi;
      }
      else {
	g_crit[1] = 0.0;
      }
    }
    free(curv);
    g_crit_computed = 1;
  }
  
  //move starting points away from gamma_crit!
  //1e-4 worked fine in the few cases I tested...
  double lo_g = g_crit[0]+1e-4;
  double hi_g = g_crit[1]-1e-4;
  
  
  if(hi_g<=0){
   freeParameter(&p2);
   return -1; //upper bound of gamma is zero: the system is not bistable
  } else {
    estim[2] = hi_g;
  }
  if(lo_g<=0) {
    lo_g = 1e-4;
    estim[1]=1e-4;
  } else {
    estim[1] = lo_g;
  }
  //printf("lo = %f,hi = %f\n",lo_g,hi_g);  
  p2.gamma = lo_g;
  if(f_infection(x0,dy0,dz0,tmax,maxdist,p2)==0) {
    if(g_crit_computed==1) free(g_crit);    
    if(g_crit[0]>0) {
      for(int i=0;i<3;i++) estim[i] = lo_g;
    } else {
      estim[0] = -INFINITY;
    }
    freeParameter(&p2);
    //printf("LOW\n");
    return 1; //the infection does not start even with the lowest possible gamma
  }
  
  p2.gamma = hi_g;
  if(f_infection(x0,dy0,dz0,tmax,maxdist,p2)==1) {
    if(g_crit_computed==1) free(g_crit);
    for(int i=0;i<3;i++) estim[i] = hi_g;
    freeParameter(&p2);
    //printf("HIGH\n");
    return 1; //the infection starts even with the highest possible gamma
  }


  //start bisection
  int iter = 0;  
  //printf("%d %g %g\n",iter,lo_g,hi_g);
  while(hi_g-lo_g>eps & iter<maxiter) {
    p2.gamma = (lo_g+hi_g)/2.0;
    if(f_infection(x0,dy0,dz0,tmax,maxdist,p2)==1) {
      lo_g = p2.gamma;
    } else {
      hi_g = p2.gamma;
    }
    iter++;
    //printf("%d %g %g\n",iter,lo_g,hi_g);
  }
  
  estim[0] = (hi_g+lo_g)/2.0;
  estim[1] = lo_g;
  estim[2] = hi_g;

  if(g_crit_computed==1) free(g_crit);
  freeParameter(&p2);
  if(iter>=maxiter) return 0;
  return 1;
}

//alpha infection
int alpha_infection(double x0,double dy0,double dz0,double* a_crit, double tmax,double maxdist,double eps,int maxiter,double* estim,Parameter p) {

  Parameter p2;
  copyParameter(p,&p2); 
    
  int a_crit_computed = 0;
  if(a_crit == NULL) {
    //alpha values that bound the bistable region need to be computed
    a_crit = (double*) malloc(10*sizeof(double));
    int    *curv   = (int*) malloc(10*sizeof(int));
    int nb = alpha_crit (10,a_crit,curv,p);

    if(nb>2) {
      //should never happen!
      printf("alpha_infection: More than two critical values of alpha have been found!\n");
      freeParameter(&p2);
      return -2;
    }
    if(nb < 1){
     freeParameter(&p2);
     return -1;
     } //the system is not bistable
    if(nb<2) {
      double a_hi=-1.0,a_lo=-1.0;
      for(int i=0;i<nb;i++) {
	if(curv[i]==1)  a_hi=a_crit[i];
	if(curv[i]==-1) a_lo=a_crit[i];
	if(curv[i]==0) {
	  printf("WARNING! Cannot tell whether alpha_crit is a maximum or a minimum!\n");
	}
      }
      if(a_lo>=0) {
	a_crit[0] = a_lo;
      } else {
	  a_crit[0] = 0.0;
      }
      if(a_hi>=0) {
	a_crit[1] = a_hi;
      }
      else {
	a_crit[1] = 0.0;
      }
    }
    free(curv);
    a_crit_computed = 1;
  }
  
  //move starting points away from alpha_crit!
  //1e-4 worked fine in the few cases I tested...
  double lo_a = a_crit[0]+1e-4;
  double hi_a = a_crit[1]-1e-4;
 
  if(hi_a<=0){
   freeParameter(&p2);
   return -1; //upper bound of alpha is zero: the system is not bistable
  } else {
    estim[2] = hi_a;
  }
  if(lo_a<=0) {
    lo_a = 0.0;
    estim[1]=0.0;
  } else {
    estim[1] = lo_a;
  }
    
  p2.alpha = lo_a;
  if(f_infection(x0,dy0,dz0,tmax,maxdist,p2)==0) {
    if(a_crit[0]>0) {
      for(int i=0;i<3;i++) estim[i] = lo_a;
    } else {
      estim[0] = -INFINITY;
    }
    if(a_crit_computed==1) free(a_crit);    
    freeParameter(&p2);
    return 1; //the infection does not start even with the lowest possible alpha
  }
  
  p2.alpha = hi_a;
  if(f_infection(x0,dy0,dz0,tmax,maxdist,p2)==1) {
    for(int i=0;i<3;i++) estim[i] = hi_a;
    if(a_crit_computed==1) free(a_crit);
    freeParameter(&p2);
    return 1; //the infection starts even with the highest possible alpha
  }


  //start bisection
  int iter = 0;  
  //printf("%d %g %g\n",iter,lo_a,hi_a);
  while(hi_a-lo_a>eps & iter<maxiter) {
    p2.alpha = (lo_a+hi_a)/2.0;
    if(f_infection(x0,dy0,dz0,tmax,maxdist,p2)==1) {
      lo_a = p2.alpha;
    } else {
      hi_a = p2.alpha;
    }
    iter++;
    //printf("%d %g %g\n",iter,lo_a,hi_a);
  }
  
  estim[0] = (hi_a+lo_a)/2.0;
  estim[1] = lo_a;
  estim[2] = hi_a;

  if(a_crit_computed==1) free(a_crit);
  freeParameter(&p2);
  if(iter>=maxiter) return 0;
  return 1;
}

int f_spbl(double x0,double dlogy0,double dz0,double tmax,double* eps,int maxiter,Parameter p) {
  int nbeq = 0;
  double** eq = equilibrium(&nbeq,p);
    
/*  printf("alpha = %g, %d equilibria\n",p.alpha,nbeq);*/
  if(nbeq<3) printf("WARNING f_spbl should not be called on systems with only %d equilibria!\n",nbeq);
  
  //set starting conditions
  double* var = (double*) malloc(7*sizeof(double));
  var[0] = x0;
  var[1] = eq[1][0]*exp(dlogy0);
  var[2] = eq[2][0]+dz0*p.zcrit;
  for(int i=3;i<7;i++) var[i] = 0.0;

  //spbl
  int test = spbl(var,eq,nbeq,tmax,eps,maxiter,p);
  
/*  for(int i = 0; i<7;i++)printf("%f ",var[i]);*/
/*  printf("\n%d\n",test);*/
  for(int i=0;i<6;i++) free(eq[i]);  
  free(eq);
  free(var);

  if(test>0) return 1;
  return 0;
}

int alpha_spbl(double* estim,Parameter p,double x0,double dlogy0,double dz0,double lo_a,double hi_a,double tmax,double eps,int maxiter) {
  if(lo_a>hi_a || lo_a<0 || hi_a<0) return -1;
    
  Parameter p2;
  copyParameter(p,&p2);
  
  double eps_tmp = eps;
  int test;
  
  p2.alpha = hi_a;
  test = f_spbl(x0,dlogy0,dz0,tmax,&eps_tmp,maxiter,p2);
  if(test==0) {
    estim[0] = NAN;
    estim[1] = lo_a;
    estim[2] = hi_a;
    freeParameter(&p2);
    return 0;
  }

  p2.alpha = lo_a;
  eps_tmp = eps;
  test = f_spbl(x0,dlogy0,dz0,tmax,&eps_tmp,maxiter,p2);  
  if(test==1) {
    estim[0] = lo_a;
    estim[1] = lo_a;
    estim[2] = lo_a;
    freeParameter(&p2);
    return 1;
  }

  //start bisection
  int iter = 0;  
  //printf("%d %g %g\n",iter,lo_a,hi_a);
  while(hi_a-lo_a>eps & iter<maxiter) {
    p2.alpha = (lo_a+hi_a)/2.0;
    eps_tmp = eps;
    if(f_spbl(x0,dlogy0,dz0,tmax,&eps_tmp,maxiter,p2)==0) {
      lo_a = p2.alpha;
    } else {
      hi_a = p2.alpha;
    }
    iter++;
    //printf("%d %g %g\n",iter,lo_a,hi_a);
  }
  
  estim[0] = (hi_a+lo_a)/2.0;
  estim[1] = lo_a;
  estim[2] = hi_a;
  freeParameter(&p2);
  return 1;
}


int alpha_blud(double* estim,double target_blud,Parameter p,double tmax,double eps,int maxiter) {
  double ak = alpha_kill(p);
  if(ak<2*eps || ak==NAN) return -1;

  double min_blud = (p.ksi*p.zcrit-p.eta/p.delta)/(p.omega-p.eta/p.delta);
  
  if(target_blud<min_blud) {
    printf("WARNING alpha_blud : target_blud is below minimum blud value %f!\n",min_blud);
    estim[0] = ak;
    estim[1] = NAN;
    estim[2] = NAN;
    return -2;    
  }
    
  Parameter p2;
  copyParameter(p,&p2);

  double f_alpha_blud(double a, void *vp){
    Parameter* p = (Parameter*) vp;
    p->alpha = a;
    double* var = (double*) malloc(4*sizeof(double));
    for(int i=0;i<4;i++) var[i] = -1.0;
    double eps_tmp = eps;

    int code = blud(var,NULL,-1,tmax,&eps_tmp,maxiter,*p);
    double b = var[0];
    free(var);

    //printf("%f -- %d %f\n",a,code,b);
    
    if(code<0) return NAN;
    return b-target_blud;
  }

  void gsl_permisive_error_handler(const char * reason,
				      const char * file,
				      int line,
				      int gsl_errno) {
    printf("\nError while computing alpha_blud: %s\n", gsl_strerror (gsl_errno));
    printf("%s\n\n",reason);
  }
  gsl_error_handler_t*  old_handler = gsl_set_error_handler (&gsl_permisive_error_handler);

  double lo = eps;
  double hi = ak-eps;
  
  double b_lo = f_alpha_blud(lo,&p2);
  if(b_lo==NAN)  {
    freeParameter(&p2);    
    return -3;
  }
  
  if(fabs(b_lo)<eps) {
    estim[0]=lo;
    estim[1]=0;
    estim[2]=lo;
    freeParameter(&p2);    
    return 2;    
  }

  double b_hi = f_alpha_blud(hi,&p2);
  if(b_hi==NAN)  {
    freeParameter(&p2);    
    return -3;
  }
  
  if(fabs(b_hi)<eps) {
    estim[0]=hi;
    estim[1]=hi;
    estim[2]=ak;
    freeParameter(&p2);    
    return 2;    
  }

  if(b_lo<0) {
    printf("WARNING alpha_blud : target_blud seems to exceed maximum possible blud value %f!\n",b_lo+target_blud);
    estim[0]=NAN;
    estim[1]=NAN;
    estim[2]=NAN;
    freeParameter(&p2);    
    return -4;        
  }
  
  gsl_function F_ab;
  F_ab.function = &f_alpha_blud;
  F_ab.params = &p2;

  const gsl_root_fsolver_type *Tr = gsl_root_fsolver_brent;
  gsl_root_fsolver *sr = gsl_root_fsolver_alloc (Tr);
  gsl_root_fsolver_set (sr, &F_ab, lo, hi);

  double atmp=hi;  
  int status,iter=0;
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (sr);
      atmp = gsl_root_fsolver_root (sr);
      lo = gsl_root_fsolver_x_lower (sr);
      hi = gsl_root_fsolver_x_upper (sr);
      //printf("%d/%d -- %f %f %f -- %d\n",iter,maxiter,atmp,lo,hi,status);
      status = gsl_root_test_interval (lo,hi,0,1e-3);
    }
  while (status == GSL_CONTINUE && iter < maxiter);
  gsl_set_error_handler (old_handler);

  estim[0] = atmp;
  estim[1] = lo;
  estim[2] = hi;
  
  freeParameter(&p2);
  gsl_root_fsolver_free (sr);

  if(iter >= maxiter) return -3;
  
  if (status == GSL_SUCCESS) return 1;
  return 0;
}

//compute the derivative of x(t) according to x(0) and derive
//direct Lyapunov Exponent (DLE) from it
//
double DLEunivar(double x,double y,double z,double tmax,double eps,char who,Parameter* p) {
  
  double lvar;
  switch(who) {
  case 'x':
    lvar = log(x);
    break;
  case 'y':
    lvar = log(y);
    break;
  case 'z':
    lvar = log(z);
    break;
  default:
    printf("Error in DLEunivar: variable %c does not exist!\n",who);
  }
  
  double h = eps/2.0;

  double** var = (double**) malloc(4*sizeof(double*)); 
  for(int i=0;i<4;i++) var[i] = (double*) malloc(4*sizeof(double));
  double* deriv = (double*) malloc(3*sizeof(double));

  //compute the four neighboring trajectories
  for(int j=0;j<4;j++) {
    var[j][0] = x;
    var[j][1] = y;
    var[j][2] = z;
    var[j][3] = 0.0;
    int k = j<2?(j-2):(j-1);
    double lvar2 = lvar+ ((double) k)*h;
    switch(who) {
    case 'x':
      var[j][0] = exp(lvar2);
      break;
    case 'y':
      var[j][1] = exp(lvar2);
      break;
    case 'z':
      var[j][2] = exp(lvar2);
      break;
    }
    dyninf(var[j],tmax,0,1,p);
    //derivatives will be computed on log-transformed variables?
    //probably necessary if dle needs to be computed with one variable close to zero
    for(int i=0;i<3;i++) var[j][i] = log(var[j][i]);
  }

  //compute derivative approximations for each variable, and Lyapunov exponent
  //Estimating the error would require to compute a 3 point approx, and therefore an additionnal
  //fifth trajectory, starting from x.
  double dle = 0.0;
  for(int i=0;i<3;i++) {
    deriv[i]  = (var[0][i]-8.0*var[1][i]+8.0*var[2][i]-var[3][i])/(12.0*h); //5 points approx
    dle += deriv[i]*deriv[i];
  }

  free(deriv);
  for(int i=0;i<4;i++) free(var[i]);
  free(var);
  
  return log(dle);
}

double DLE(double x,double y,double z,double tmax,double eps,Parameter* p) {
  double h = eps/2.0;

  double** J =  (double**) malloc(3*sizeof(double*));
  for(int i=0;i<3;i++) J[i] = (double*) malloc(3*sizeof(double));
  
  double** var = (double**) malloc(4*sizeof(double*)); 
  for(int i=0;i<4;i++) var[i] = (double*) malloc(4*sizeof(double));

  double lvar;
  for(int who=0;who<3;who++) {
    switch(who) {
    case 0:
      lvar = log(x);
      break;
    case 1:
      lvar = log(y);
      break;
    case 2:
      lvar = log(z);
      break;
    }  
    //compute the four neighboring trajectories
    for(int j=0;j<4;j++) {
      var[j][0] = x;
      var[j][1] = y;
      var[j][2] = z;
      var[j][3] = 0.0;
      int k = j<2?(j-2):(j-1);
      double lvar2 = lvar+ ((double) k)*h;
      var[j][who] = exp(lvar2);      
      dyninf(var[j],tmax,0,1,p);      
      //derivatives will be computed on log-transformed variables?
      //probably necessary if dle needs to be computed with one variable close to zero
      for(int i=0;i<3;i++) var[j][i] = log(var[j][i]);
    }

  //compute derivative approximations for each variable, and Lyapunov exponent
  //Estimating the error would require to compute a 3 point approx, and therefore an additionnal
  //fifth trajectory, starting from x.
    for(int i=0;i<3;i++) {
      J[i][who]  = (var[0][i]-8.0*var[1][i]+8.0*var[2][i]-var[3][i])/(12.0*h); //5 points approx
    }
  }

  //compute tr(J) * J and store it in L
  double* L =  (double*) malloc(3*3*sizeof(double));
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      L[i*3+j]=0.0;
      for(int k=0;k<3;k++) {
	L[i*3+j] += J[k][i]*J[k][j];
      }
    }
  }
  
  for(int i=0;i<4;i++) free(var[i]);
  free(var);
  for(int i=0;i<3;i++) free(J[i]);
  free(J);

  gsl_vector_complex *eval = gsl_vector_complex_alloc (3);		 
  gsl_eigen_nonsymm_workspace * w = gsl_eigen_nonsymm_alloc (3);
  gsl_matrix_view m = gsl_matrix_view_array (L,3,3);

  void gsl_permisive_error_handler(const char * reason,
				      const char * file,
				      int line,
				      int gsl_errno) {
    printf("\nError while DLE: %s\n", gsl_strerror (gsl_errno));
    printf("%s\n\n",reason);
  }
  gsl_error_handler_t*  old_handler = gsl_set_error_handler (&gsl_permisive_error_handler);
  int test  = gsl_eigen_nonsymm (&m.matrix, eval, w);                //eigenvalues are stored in eval
  gsl_set_error_handler (old_handler);

  free(L);
  if(test != GSL_SUCCESS) {
      return NAN;
  }
    
  gsl_vector_view re_ev = gsl_vector_complex_real(eval); //retrieve real part
  double dle = gsl_vector_get(&(re_ev.vector),0);
  for(int i=1;i<3;i++) {
    double dle_tmp = gsl_vector_get(&(re_ev.vector),i);
    if(dle_tmp>dle) dle=dle_tmp;
  }
  return log(dle)/(2.0*tmax);
}


void DLEmatrix(double* L,double x,double y,double z,double tmax,double eps,Parameter* p) {
  double h = eps/2.0;

  double** J =  (double**) malloc(3*sizeof(double*));
  for(int i=0;i<3;i++) J[i] = (double*) malloc(3*sizeof(double));
  
  double** var = (double**) malloc(4*sizeof(double*)); 
  for(int i=0;i<4;i++) var[i] = (double*) malloc(4*sizeof(double));

  double lvar;
  for(int who=0;who<3;who++) {
    switch(who) {
    case 0:
      lvar = log(x);
      break;
    case 1:
      lvar = log(y);
      break;
    case 2:
      lvar = log(z);
      break;
    }  
    //compute the four neighboring trajectories
    for(int j=0;j<4;j++) {
      var[j][0] = x;
      var[j][1] = y;
      var[j][2] = z;
      var[j][3] = 0.0;
      int k = j<2?(j-2):(j-1);
      double lvar2 = lvar+ ((double) k)*h;
      var[j][who] = exp(lvar2);      
      dyninf(var[j],tmax,0,1,p);      
      //derivatives will be computed on log-transformed variables?
      //probably necessary if dle needs to be computed with one variable close to zero
      for(int i=0;i<3;i++) var[j][i] = log(var[j][i]);
    }

  //compute derivative approximations for each variable, and Lyapunov exponent
  //Estimating the error would require to compute a 3 point approx, and therefore an additionnal
  //fifth trajectory, starting from x.
    for(int i=0;i<3;i++) {
      J[i][who]  = (var[0][i]-8.0*var[1][i]+8.0*var[2][i]-var[3][i])/(12.0*h); //5 points approx
    }
  }

  //compute tr(J) * J and store it in L
  for(int i=0;i<3;i++) {
    for(int j=0;j<3;j++) {
      L[i*3+j]=0.0;
      for(int k=0;k<3;k++) {
	L[i*3+j] += J[k][i]*J[k][j];
      }
    }
  }
  
  for(int i=0;i<4;i++) free(var[i]);
  free(var);
  for(int i=0;i<3;i++) free(J[i]);
  free(J);
}

double threshold_spbl_x(double dy,double dz,double tmax,double* eps,int maxiter,Parameter p) {    
  int nbeq = -1;
  double** eq = equilibrium(&nbeq,p);
  if(nbeq<3) return NAN;
  
  double eps_old = *eps;

  double x,lo_x,hi_x;
  double y = eq[1][0]+dy;
  double z = eq[2][0]+dz;
  /* printf("%g %g\n",y,z); */
  if(threshold_x(1,&x,&lo_x,&hi_x,&y,&z,tmax,*eps,0.001,p)!=1) return NAN;

  double* var = (double*) malloc(7*sizeof(double));
  double lo_lx = log(x+1e-3);
  double hi_lx = log(eq[0][nbeq-1]-1e-6);

  var[0] = exp(lo_lx);
  var[1] = eq[1][0]+dy;
  var[2] = eq[2][0]+dz;
  var[3] = 0.0;
  *eps = eps_old;
  int code_lo = spbl(var,eq,nbeq,tmax,eps,maxiter,p);
  if(code_lo!=1) {//spbl is not reached even for the lowest value of x
                  //return this lowest value, assuming that spbl is necessarily 
                  //reached when the initial point lies on the separatrix.
    free(var);
    for(int i=0;i<6;i++) free(eq[i]); free(eq);
    return exp(lo_lx);
  }
    
  var[0] = exp(hi_lx);
  var[1] = eq[1][0]+dy;
  var[2] = eq[2][0]+dz;
  var[3] = 0.0;
  *eps = eps_old;
  int code_hi = spbl(var,eq,nbeq,tmax,eps,maxiter,p);
  if(code_hi==1) {//spbl is  reached even for the highest value of x
                  //return this higest value
    free(var);
    for(int i=0;i<6;i++) free(eq[i]); free(eq);
    return exp(hi_lx);
  }

  int iter = 0;
  while(hi_lx-lo_lx>eps_old && iter<maxiter) {
    //printf("%d -- %g %g\n",iter,exp(lo_lx),exp(hi_lx));
    double lx = (lo_lx+hi_lx)/2.0;
    var[0] = exp(lx);
    var[1] = eq[1][0]+dy;
    var[2] = eq[2][0]+dz;
    var[3] = 0.0;
    *eps = eps_old;
    int code = spbl(var,eq,nbeq,tmax,eps,maxiter,p);
    //printf("%g -- %g %g %g -- %d\n\n",exp(lx),var[0],var[1],var[2],code);
    if(code==1 & var[2] < p.zcrit) {
      lo_lx = lx;
    } else {
      hi_lx = lx;
    }
    iter++;
  }
  
    free(var);
    for(int i=0;i<6;i++) free(eq[i]); free(eq);
    return exp((lo_lx+hi_lx)/2.0);
}

/* returns code that tells why BLUD could not be computed */
/* 2 -- SPBL has been found but closest eq is stable */
/* 1 -- SPBL has been computed */
/* 0 -- x does not exceeed local peak before tmax*/
/* -1 -- no local minimum has been found before tmax*/
/* -2 -- local maximum has not been found before tmax*/
/* -3 -- the system is not bistable*/
int spbl(double* var,double** eq,int nbeq,double tmax,double* eps,int maxiter,Parameter p) {
  int compute_eq = 0;
  if(!eq) {
    //equilibria are computed when not provided by calling function
    nbeq = -1;
    eq = equilibrium(&nbeq,p);
    compute_eq = 1;
  }

  int nbsteq = 0;
  EVENT_xe = -1.0; //initialization
  int poseq = -1;
  for(int i=0;i<nbeq;i++) {    
    if(eq[3][i]<0) {      
      nbsteq++;
    } else {
      //position of the unstable equilibrium with x>0
      if(EVENT_xe<eq[0][i]) {
	poseq = i;
	EVENT_xe = eq[0][poseq];
	EVENT_ye = eq[1][poseq];
	EVENT_ze = eq[2][poseq];
      }
    }
  }
  
  if(nbsteq!=2) {
    //the system is not bistable
    if(compute_eq==1) {
      for(int i=0;i<6;i++) free(eq[i]);
      free(eq);
    }
    return(-3);
  }

  // 1. compute the first local maximum
  var[3]= 0.0; //var[3] is time!
  double eps_old = *eps;
  EVENT_xstart = var[0];
  EVENT_ystart = var[1];
  EVENT_zstart = var[2];
  EVENT_mindist = 1e-6;
  int test = dyninfWEvent(var,tmax,eps,maxiter,&xlocalmax,1,&p);
/*   for(int i=0;i<4;i++) printf("%g ",var[i]); */
/*   printf(" -- %d\n",test); */

  if(test!=1) return -2; // local maximum has not been found before tmax
  var[4] = var[0]; //time at peak
  var[5] = var[3]; //load at peak

  // 2. compute the following local minimum
  *eps = eps_old; 
  EVENT_xstart = var[0];
  EVENT_ystart = var[1];
  EVENT_zstart = var[2];
  test = dyninfWEvent(var,tmax,eps,maxiter,&xlocalmin,1,&p);
/*   for(int i=0;i<4;i++) printf("%g ",var[i]); */
/*   printf(" -- %d\n",test); */
  if(test!=1) return -1; // no local mininum found before tmax 

  //check that point 2 is closest to the unstable equilibrium  
  double dist,mindist;
  int posmindist = -1;
  for(int i=0;i<nbeq;i++) {
    dist = 0.0;
    for(int j=0;j<3;j++) dist += pow(var[j]-eq[j][i],2.0);
    dist = sqrt(dist);
    if(posmindist==-1 || mindist>dist) {
      mindist = dist;
      posmindist = i;
    }
  }
  
  //if(posmindist != poseq) return -1; //the closest equilibrium clearance or chronic infection
  if(posmindist < poseq) return -1; //the closest equilibrium clearance or chronic infection

  double* var_old = (double*) malloc(4*sizeof(double));
  for(int i=0;i<4;i++) var_old[i] = var[i];


  // 3. compute the point at which loads exceeds the first peak
  *eps = eps_old; 
  EVENT_x = var[4];  
  test = dyninfWEvent(var,tmax,eps,maxiter,&logxcritreached,1,&p);
/*   for(int i=0;i<4;i++) printf("%g ",var[i]); */
/*   printf(" -- %d\n",test); */
/*   */
   
  var[6] = var[3]-var[5];
  for(int i=0;i<4;i++) var[i] = var_old[i];
  free(var_old);
  if(test!=1) return 0;//x does not reach load at peak before tmax
  
  if(posmindist > poseq) return 2; //SPBL found, but closest eq is stable
  return 1; //SPBL is OK!

}

/* returns code that tells why BLUD could not be computed */
/* 0 -- no BLUD computed, because the host did not die */
/* 1 -- BLUD has been computed and lies below equilibrium x */
/* 2 -- BLUD has been computed but is above equilibrium x */
/* -1 -- equilibria could not be computed (check ksi value ??)*/
/* -2 -- BLUD could not be computed, because the host did not die. Equilibrium z is below zcrit */
/* -3 -- no unstable equilibrium could be used as starting condition */
/* -4 -- the starting z value is above zcrit */
int blud(double* var,double** eq,int nbeq,double tmax,double* eps,int maxiter,Parameter p) {
  int compute_eq = 0;
  if(!eq) {
    //equilibria are computed when not provided by calling function
    nbeq = -1;
    eq = equilibrium(&nbeq,p);
    compute_eq = 1;
  }
  
  if(nbeq<1) {
    //problem with equilibria computation??
    //check value of ksi!!
    if(compute_eq==1) {
      for(int i=0;i<6;i++) free(eq[i]);
      free(eq);
    }
    return(-1);
  }
  
  //starting point must be one of the unstable equilibria
  // * clearance if the system is not bistable
  // * the unstable chronic infection otherwise
  int starting_eq = nbeq-1;
  while(eq[3][starting_eq]<0) starting_eq--;
  if(starting_eq<0) {
    if(compute_eq==1) {
      for(int i=0;i<6;i++) free(eq[i]);
      free(eq);
    }
    return(-3);
  }
  if(eq[2][starting_eq]>=p.zcrit) {
    //the starting equilibrium is above zcrit
    //BLUD is no defined, which does not mean that the host cannot die from the infection
    if(compute_eq==1) {
      for(int i=0;i<6;i++) free(eq[i]);
      free(eq);
    }
    return(-4);
  }
  
  //compute perturbation
  //If starting equilibrium is clearance, the perturbation can easily be computed
  //by hand as the dominant eigenvalue is 1-\delta y and the associated eigenvector is (1,0,0).
  //Trying to compute it from the jacobian will anyway fail if u<1 or v<1.

  double* y = malloc(3*sizeof(double)); 
  for(int j=0;j<3;j++) y[j] = eq[j][starting_eq];
  
  if(starting_eq==0) {
    y[0] = 1.0;
    y[1] = 0.0;
    y[2] = 0.0;
  } else {    
    gsl_vector_complex *eval = gsl_vector_complex_alloc (3);
    gsl_matrix_complex *evec = gsl_matrix_complex_alloc (3, 3);
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (3);
    gsl_eigen_nonsymmv_params(1,w);

    double *dfdy = malloc(9*sizeof(double));    
    if(compute_jacob_dxyz_dt(y,dfdy,1,p)==0) {
      gsl_matrix_view m = gsl_matrix_view_array (dfdy,3,3);    
      gsl_eigen_nonsymmv (&m.matrix, eval, evec, w); //eigenvalues are stored in eval and eigenvectors in evec
      gsl_eigen_nonsymmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC);
      //check that last vector goes with the only eigen value with positive real part
      int pos_ev = 2;    
      gsl_complex eval_i = gsl_vector_complex_get (eval, pos_ev);
      while(GSL_REAL(eval_i)<0 && pos_ev>0) {
	pos_ev--;
	eval_i = gsl_vector_complex_get (eval, pos_ev);
      }
      if(GSL_REAL(eval_i)<0) {
#ifdef RMOD
	warning("ERROR in BLUD computation: eigenvalue should have positive real part!");
#else
	printf("ERROR in BLUD computation: eigenvalue should have positive real part!\n");
#endif
      }
      gsl_vector_complex_view evec_i = gsl_matrix_complex_column (evec,pos_ev);
      for (int j=0;j<3;j++) {
	gsl_complex z = gsl_vector_complex_get(&evec_i.vector, j);
	y[j]=GSL_REAL(z);
      }
    }
    free(dfdy);
    gsl_eigen_nonsymmv_free(w);
    gsl_vector_complex_free(eval);
    gsl_matrix_complex_free(evec);
  }
  
  if(y[2]<0) {
    //move in the direction of increasing z
    for(int j=0;j<3;j++) y[j]*=-1;
  }
    
  for(int i=0;i<3;i++) var[i] = eq[i][starting_eq]+y[i]*(1e-4); //!!!!!!!!!!!
  var[3] = 0.0;
  
  free(y);

  //method 1: one call to dyninf and another to tsurv
  //the existence of BLUD being unguaranteed, direct integration over z may fail.
  int test = dyninf(var,tmax,1,1,&p);
  //printf("x=%g y=%g z=%g t=%g\n",var[0],var[1],var[2],var[3]);
  if(test==1) { //the host is dead
    tsurv(var,&p);
    if(eq[0][nbeq-1]<var[0]) {
      test = 2; //BLUD is above equilibrium x
    } 
  } else { //the host survived
    if(eq[2][nbeq-1]<p.zcrit) test = -2; //no BLUD and equilibrium z is below zcrit
  }

  if(compute_eq==1) {
    for(int i=0;i<6;i++) free(eq[i]);
    free(eq);
  }
  
  return(test);
}


double getLT(int* dead,double x,double tmax,double** eq,Parameter* p) {
    double* var = (double*) malloc(4*sizeof(double));
    var[0] = x;
    var[1] = eq[1][0]; //y homeostatic state
    var[2] = eq[2][0]; //z homeostatic state
    var[3] = 0.0;

    double t = tmax;
    *dead = dyninf(var,tmax,1,0,p); //stops the dynamics as soon as z exceed zcrit
    if(*dead==1){ //the host has died: compute time at death
      if(tsurv(var,p)!= GSL_SUCCESS) *dead = -1;
      t = var[3];
    }
    free(var);
    if(*dead==0) return(NAN);
    return(t);
}

double logLD(double logxmin,double logxmax,double tmax,double eps,double** eq,Parameter* p) {
  double logx_lo  = logxmin;
  double logx_hi  = logxmax;
  int     dead_lo,dead_hi;
  
  double t_hi = getLT(&dead_hi,exp(logx_hi),tmax,eq,p);
  //printf("%g %d\n",t_hi,dead_hi);
  if(dead_hi==0) return(NAN);
  
  double t_lo = getLT(&dead_lo,exp(logx_lo),tmax,eq,p);
  //printf("%g %d\n",t_lo,dead_lo);
  if(dead_lo==1) return(logxmin);

  while(fabs(logx_lo-logx_hi)>eps) {
    int    dead_mid;
    double logx_mid = (logx_lo+logx_hi)/2.0;
    double t_mid = getLT(&dead_mid,exp(logx_mid),tmax,eq,p);
    //printf("%g %g %d -- ",logx_mid,t_mid,dead_mid);
    if(dead_mid==1) {
      logx_hi = logx_mid;
    } else {
      logx_lo = logx_mid;
    }
    //printf("[%g;%g] %g %g\n",logx_lo,logx_hi,fabs(logx_lo-logx_hi),eps);
  }
  return((logx_lo+logx_hi)/2.0);
}

double avgLT(double tmax,double** eq,Parameter *p,StartingCondition* sc) {
  double gg(double logx,void* param) {
    (void*) param;
    int s;
    double lt = getLT(&s,exp(logx),tmax,eq,p);
    //printf("log(x)=%g LT=%g",logx,lt);
    
    if(lt!=lt) return(NAN);
    double fx = gsl_ran_gaussian_pdf(logx - sc->mu_logx,sc->sigma_logx);
    //printf(" f=%g\n",fx);
    return(lt*fx);
  }

  
  double logx_hi = 0.0;
  double logx_lo = sc->mu_logx + gsl_cdf_gaussian_Pinv(1e-12,sc->sigma_logx);  
  logx_lo = logLD(logx_lo,0.0,tmax,1e-6,eq,p);
  
  if(logx_lo!=logx_lo) {
    //the host does not die, even at highest dose
    return(NAN);
  }
    
  //probability to die
  double pr = 1.0 - gsl_cdf_gaussian_P(logx_lo-sc->mu_logx,sc->sigma_logx);
  if(pr<=1e-48) {
    //the host could die, but at doses which are highly unlikely
    //time to death should then be infinite
    return(INFINITY);
  }

  void gsl_permisive_error_handler(const char * reason,
				 const char * file,
				 int line,
				 int gsl_errno) {
    printf("\nError while computing average LT: %s\n", gsl_strerror (gsl_errno));
    printf("%s\n\n",reason);    
  }
  gsl_error_handler_t*  old_handler = gsl_set_error_handler (&gsl_permisive_error_handler);

  //numerical integration of ff on the distribution of xa
  gsl_integration_workspace * w  = gsl_integration_workspace_alloc (1000);
  double result = 0.0,error = 0.0;

  gsl_function FF;
  FF.function = &gg;
  FF.params = &p;
  
  gsl_integration_qag (&FF,logx_lo,logx_hi, 0, 1e-6, 1000, 6, w, &result, &error);  
  gsl_integration_workspace_free (w);
  gsl_set_error_handler(old_handler);
 
  result /= pr;
  
  return(result);
}

