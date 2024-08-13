#include <stdlib.h>
#include <stdio.h>

#include <unistd.h>

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>

#include <R.h>
#include <Rmath.h>

#include "WHD_model.h"

#define RMOD = 1

/*function taking as input 'par' vector from R and make struct type object of named parameters*/
void initParameter(Parameter* p, double *par) {
	p->alpha = par[0];
	p->beta = par[1];
	p->gamma = par[2];
	p->delta = par[3];
	p->eta = par[4];
	p->ksi = par[5];
	p->phi = par[6];
	p->psi = par[7];
	p->omega = par[8];
	p->u = par[9];
	p->v = par[10];
	p->w = par[11];
	p->l = par[12];
	p->zcrit = par[13];
	
	initializeParameter(p);
}

/*same function but for the vector of initial conditions*/
void initStartingCondition(StartingCondition* sc, double *par_ini,int logscale){
  if(logscale==1) {
    sc->mu_logx = par_ini[0];
    sc->mu_logy = par_ini[1];
    sc->mu_logz = par_ini[2];
    sc->mu_x = exp(par_ini[0]);
    sc->mu_y = exp(par_ini[1]);
    sc->mu_z = exp(par_ini[2]);
  } else {    
    sc->mu_x = par_ini[0];
    sc->mu_y = par_ini[1];
    sc->mu_z = par_ini[2];
    sc->mu_logx = log(par_ini[0]);
    sc->mu_logy = log(par_ini[1]);
    sc->mu_logz = log(par_ini[2]);
  }
  sc->sigma_logx = par_ini[3];	
  sc->sigma_logy = par_ini[4];
  sc->sigma_logz = par_ini[5];
}


/*function to sample initial conditions*/

double getXini(StartingCondition sc){
  GetRNGstate();
  double res = sc.mu_logx + norm_rand()*sc.sigma_logx;  
  PutRNGstate();
  return (exp(res));
}

double getYini(StartingCondition sc){
  GetRNGstate();
  double res = sc.mu_logy + norm_rand()*sc.sigma_logy;
  PutRNGstate();
  return (exp(res));
}

double getZini(StartingCondition sc){
  GetRNGstate();
  double res = sc.mu_logz + norm_rand()*sc.sigma_logz;
  PutRNGstate();
  return (exp(res));
}

void _F(int* n,double* x,double* f, double *par) {
  Parameter p;
  initParameter(&p,par);
  for(int i=0;i<*n;i++) f[i] = F(x[i],p);
  freeParameter(&p);
}

void _G(int* n,double* y,double* z,double* g, double *par) {
  Parameter p;
  initParameter(&p,par);
  for(int i=0;i<*n;i++) g[i] =  G(y[i],z[i],p);
  freeParameter(&p);
}

void _dFdx(int* n,double* x,double* f, double *par) {
  Parameter p;
  initParameter(&p,par);
  for(int i=0;i<*n;i++) {
    f[i] = dFdx(x[i],p);
    //printf("%d %g %g\n",i,x[i],f[i]);
  }
  freeParameter(&p);
}

void _dGdy(int* n,double* y,double* z,double* g, double *par) {
  Parameter p;
  initParameter(&p,par);
  for(int i=0;i<*n;i++) g[i] =  dGdy(y[i],z[i],p);
  freeParameter(&p);
}

void _dGdz(int* n,double* y,double* z,double* g, double *par) {
  Parameter p;
  initParameter(&p,par);
  for(int i=0;i<*n;i++) g[i] =  dGdz(y[i],z[i],p);
  freeParameter(&p);
}

void _P_gamma_eq(double* x,double* g,double* par) {
  Parameter p;
  initParameter(&p,par);
  *g = P_gamma_eq(*x,&p);
  freeParameter(&p);
}

void _dP_gamma_eq(double* x,double* dp,double* par) {
  Parameter p;
  initParameter(&p,par);
  *dp = dP_gamma_eq(*x,&p);
  freeParameter(&p);
}

void _P_alpha_eq(double* x,double* a,double* par) {
  Parameter p;
  initParameter(&p,par);
  *a = P_alpha_eq(*x,&p);
  freeParameter(&p);
}

void _dP_alpha_eq(double* x,double* dp,double* par) {
  Parameter p;
  initParameter(&p,par);
  *dp = dP_alpha_eq(*x,&p);
  freeParameter(&p);
}

void _gamma_clearance(double* gc,double* par) {
  Parameter p;
  initParameter(&p,par);
  *gc = gamma_clearance(p);
  freeParameter(&p);
}

void _gamma_survival(double* gc,double* par) {
  Parameter p;
  initParameter(&p,par);
  *gc = gamma_survival(p);
  freeParameter(&p);
}

void _gamma_kill(double* gk,double* par) {
  Parameter p;
  initParameter(&p,par);
  *gk = gamma_kill(p);
  freeParameter(&p);
}

void _alpha_kill(double* ak,double* par) {
  Parameter p;
  initParameter(&p,par);
  *ak = alpha_kill(p);
  freeParameter(&p);
}

void _gamma_crit (int* maxg,double* gcrit,int* curv,double* par) {
  Parameter p;
  initParameter(&p,par); 
  *maxg = gamma_crit(*maxg,gcrit,curv,p);
  freeParameter(&p);
}

void _alpha_crit (int* maxa,double* acrit,int* curv,double* par) {
  Parameter p;
  initParameter(&p,par); 
  *maxa = alpha_crit(*maxa,acrit,curv,p);
  freeParameter(&p);
}

void _gamma_infection (double* estim,int* test,double* x0,double* dy0,double* dz0,double* tmax,double* maxdist,double* eps,int* maxiter,double* par) {
  Parameter p;
  initParameter(&p,par); 
  *test = gamma_infection(*x0,*dy0,*dz0,NULL,*tmax,*maxdist,*eps,*maxiter,estim,p);
  freeParameter(&p);
}

void _alpha_infection (double* estim,int* test,double* x0,double* dy0,double* dz0,double* tmax,double* maxdist,double* eps,int* maxiter, double* par) {
  Parameter p;
  initParameter(&p,par); 
  *test = alpha_infection(*x0,*dy0,*dz0,NULL,*tmax,*maxdist,*eps,*maxiter,estim,p);
  freeParameter(&p);
}

void _alpha_spbl(double* estim,int* test,double* x0,double* dlogy0,double* dz0,double* lo_a,double* hi_a,double* tmax,double* eps,int* maxiter,double* par) {
  Parameter p;
  initParameter(&p,par); 
  *test = alpha_spbl(estim,p,*x0,*dlogy0,*dz0,*lo_a,*hi_a,*tmax,*eps,*maxiter);
  freeParameter(&p);
  }

void _alpha_blud(double* estim,int* test,double* target_blud,double* tmax,double* eps,int* maxiter,double*  par) {
  Parameter p;
  initParameter(&p,par); 
  *test = alpha_blud(estim,*target_blud,p,*tmax,*eps,*maxiter);
  freeParameter(&p);
}

void _H(int* n,double *x,double* h,double* par) {
  Parameter p;
  initParameter(&p,par);
  for(int i=0;i<*n;i++) {
    h[i] = H(x[i],&p);
  }
  freeParameter(&p);
}

void _dH(int* n,double *x,double* dh,double* par) {
  Parameter p;
  initParameter(&p,par);
  for(int i=0;i<*n;i++) {
    dh[i] = dH(x[i],&p);
    //printf("%g\t%g\n",x[i],dh[i]);
  }
  freeParameter(&p);
}

void _jacobian(double* eq,double* dfdy,double* pars) {
  Parameter p;
  initParameter(&p,pars);
  compute_jacob_dxyz_dt(eq,dfdy,1,p);
  freeParameter(&p);
}

void _jacobianLogx(double* eq,double* dfdy,double* pars) {
  Parameter p;
  initParameter(&p,pars);
  double t;
  double* dfdt = (double*) malloc(3*sizeof(double));
  jacob_dlogxyz_dt (t,eq,dfdy,dfdt, (void*)  &p);
  free(dfdt);
  freeParameter(&p);
}

void _dyninf(int* nb,double* x,double* y,double* z,double* t,double* par) {
  Parameter p;
  initParameter(&p,par);  
  *nb = dyninf_complete(*nb,x,y,z,t,&p);
  freeParameter(&p);
}

void _tsurv(double* x,double* y,double* z,double* t,double* tmax,double* eps,int* maxiter,double* par) {
  Parameter p;
  initParameter(&p,par);  
  
  double* var = (double*) malloc(4*sizeof(double));
  var[0] = *x;
  var[1] = *y;
  var[2] = *z;
  var[3] = 0.0; 
  int dead = dyninfWEvent(var,*tmax,eps,*maxiter,&zcritreached,0,&p);
  if(dead==1) {
    *x = var[0];
    *y = var[1];
    *z = var[2];
    *t = var[3];
  } else {
    *t = NAN;
  }
  free(var);
  freeParameter(&p);
}

void _tsurvV2(double* x,double* y,double* z,double* t,double* tmax,double* eps,int* maxiter,double* par) {
  Parameter p;
  initParameter(&p,par);  
  
  double* var = (double*) malloc(4*sizeof(double));
  var[0] = *x;
  var[1] = *y;
  var[2] = *z;
  var[3] = 0.0;
  int dead = dyninf(var,*tmax,1,0,&p); //stops integration when z>=zcrit
  if(dead==1) {
    tsurv(var,&p); //integrates over z, starting from last point of previous integration and going backward
    *x = var[0];
    *y = var[1];
    *z = var[2];
    *t = var[3];
  } else {
    *t = NAN;
  }
  free(var);
  freeParameter(&p);
}

void _dxyt_dz(double* z,double* y, double* f,double* param) {
  Parameter p;
  initParameter(&p,param);    
  dxyt_dz (*z,y,f,&p);
  freeParameter(&p);
}

void _dxyz_dt(double* t,double* y, double* f,double* param) {
  Parameter p;
  initParameter(&p,param);    
  dxyz_dt (*t,y,f,&p);
  freeParameter(&p);
}

void _survsim(int* n,double* par, double* par_ini, int* logscale,double* tmax,double* eps,int* maxiter,
	      double* x0, double* y0, double* z0, double* x, double* y, double* z, double* time, int* dead){
  Parameter p;
  StartingCondition sc;
	
  initParameter(&p,par);
  //printParameter(&p);
	
  initStartingCondition(&sc,par_ini,*logscale);
  double* var = (double*) malloc(4*sizeof(double));
  double precision;

  int nbeq = 0;
  double** eq = equilibrium(&nbeq,p);

  
  for(int i=0;i<*n;i++){  
    double x_ini= getXini(sc);
    double y_ini= eq[1][0] - getYini(sc);
    y_ini= y_ini<0.0?0.0:y_ini;
    double z_ini= eq[2][0] + getZini(sc);

		
    x0[i] = x_ini;
    y0[i] = y_ini;
    z0[i] = z_ini;

    var[0] = x_ini;
    var[1] = y_ini;
    var[2] = z_ini;
    var[3] = 0.0;  
    precision = *eps;

    dead[i] = dyninfWEvent(var,*tmax,&precision,*maxiter,&zcritreached,1,&p);
	
    /*printf("%d %g %g %d\n",i,*eps,precision,dead[i]);*/	

    x[i] = var[0];
    y[i] = var[1];
    z[i] = var[2];
    time[i] = var[3];
  }
  free(var);
  free(eq);
  freeParameter(&p);
}

void _equilibrium(int* nbeq,double *xeq,double* yeq,double* zeq,double* re_ev,double* im_ev,double* im_ev_glob,double* par) {
  Parameter p;
  initParameter(&p,par);

  double** eq = equilibrium(nbeq,p);
  for(int i=0;i<*nbeq;i++) {
    xeq[i] = eq[0][i];
    yeq[i] = eq[1][i];
    zeq[i] = eq[2][i];
    re_ev[i] = eq[3][i];
    im_ev[i] = eq[4][i];
    im_ev_glob[i] = eq[5][i];
  }
  for(int i=0;i<6;i++) free(eq[i]);
  free(eq);
  freeParameter(&p);
}

void _yhomeo(double* yeq,int* status,double* par) {
  Parameter p;
  initParameter(&p,par);
  *status =  NoInfeq(yeq,p);
  freeParameter(&p);
}

void _threshold_x(int* nb,double* x,double* lo_x, double *hi_x,double* y,double* z,double* tmax,double* eps,double* max_dist,double* par) {
  Parameter p;
  initParameter(&p,par); 
  *nb = threshold_x(*nb,x,lo_x,hi_x,y,z,*tmax,*eps,*max_dist,p);
  freeParameter(&p);
}

void _threshold_y(int* nb,double* x,double* y,double* lo_y, double *hi_y,double* z,double* tmax,double* eps,double* max_dist,double* par) {
  Parameter p;
  initParameter(&p,par); 
  *nb = threshold_y(*nb,x,y,lo_y,hi_y,z,*tmax,*eps,*max_dist,p);
  freeParameter(&p);
}

void _threshold_z(int* nb,double* x,double* y,double* z,double* lo_z, double *hi_z,double* tmax,double* eps,double* max_dist,double* par) {
  Parameter p;
  initParameter(&p,par); 
  *nb = threshold_z(*nb,x,y,z,lo_z,hi_z,*tmax,*eps,*max_dist,p);
  freeParameter(&p);
}

void _spbl(double* v,int* test,double* tmax,double* eps,int* maxiter,double* par) {
  Parameter p;
  initParameter(&p,par);  
  *test = spbl(v,NULL,-1,*tmax,eps,*maxiter,p);
  freeParameter(&p);
}

double   _threshold_spbl_x(double *lx,double* dy,double* dz,double* tmax,double* eps,int* maxiter,double* par) {
  Parameter p;
  initParameter(&p,par);  
  *lx = threshold_spbl_x(*dy,*dz,*tmax,eps,*maxiter,p);
  freeParameter(&p);
  }

void _blud(double* v,int* test,double* tmax,double* eps,int* maxiter,double* par) {
  Parameter p;
  initParameter(&p,par);
  double* var = (double*) malloc(4*sizeof(double));
  *test = blud(var,NULL,-1,*tmax,eps,*maxiter,p);
  for(int i=0;i<4;i++)  v[i] = var[i];
  free(var);
}

void _dyninf_toControl(double* v,double* x0,double* tmax,double* eps,int* maxiter,double* par) {
    Parameter p;
    initParameter(&p,par);

    //homeostatic state
    //starting condition : y is reduced from homeostatic state so that dx/dt is maximized
    double ye;
    NoInfeq(&ye,p);
    double z0 = p.eta*ye/p.ksi;
    double y0 = (1.0-2.0*(*x0))/p.delta;
    /* double ye; */
    /* NoInfeq(&ye,p); */
    /* double ze = p.eta*ye/p.ksi; */
    /* double x0 = (1.0-p.delta*ye)/2.0; */
    //if x0 is negative, control cannot occur
    if(*x0<=0.0 || y0<0.0) {
      for(int i=0;i<4;i++) v[i]=NAN;
      return;
    }

    v[0] = *x0;
    v[1] = y0;
    v[2] = z0;
    v[3] = 0.0;
    int controlled = dyninfWEvent(v,*tmax,eps,*maxiter,&xdecrease,1,&p);
    double dx = dxdt(v[0],v[1],v[2],p);
    double dy = dxdt(v[0],v[1],v[2],p);
    double dz = dxdt(v[0],v[1],v[2],p);
    double dv  = sqrt(dx*dx+dy*dy+dz*dz);
    if(dv<=*eps) controlled=0; //dxdt is zero, but other derivatives are also close to zero: an equilibrium point has been reached!

    if(controlled<1) {
      for(int i=0;i<4;i++) v[i]=NAN;
    }
    freeParameter(&p);
}

void _dyninf_toMaxResponse(double* v,double* x0,double* tmax,double* eps,int* maxiter,double* par) {
    Parameter p;
    initParameter(&p,par);

    //starting condition : y is reduced from homeostatic state so that dx/dt is maximized
    double ye;
    NoInfeq(&ye,p);
    double z0 = p.eta*ye/p.ksi;
    double y0 = (1.0-2.0*(*x0))/p.delta;
    //double x0 = (1.0-p.delta*ye)/2.0;
    //if x0 is negative, control cannot occur
    if(*x0<=0.0 || y0<0.0) {
      for(int i=0;i<4;i++) v[i]=NAN;
      return;
    }

    v[0] = *x0;
    v[1] = y0;
    v[2] = z0;
    v[3] = 0.0;
    int controlled = dyninfWEvent(v,*tmax,eps,*maxiter,&ydecrease,1,&p);
    double dx = dxdt(v[0],v[1],v[2],p);
    double dy = dxdt(v[0],v[1],v[2],p);
    double dz = dxdt(v[0],v[1],v[2],p);
    double dv  = sqrt(dx*dx+dy*dy+dz*dz);
    if(dv<=*eps) controlled=0; //dxdt is zero, but other derivatives are also close to zero: an equilibrium point has been reached!
    if(controlled<1) {
      for(int i=0;i<4;i++) v[i]=NAN;
    }
    freeParameter(&p);
}

void _DLEx(int* nb,double* dle,double* x,double* y,double* z,double* tmax,double* eps,double* par) {
   Parameter p;
   initParameter(&p,par);
   for(int i=0;i<*nb;i++) {
     dle[i] = DLEunivar(x[i],y[i],z[i],*tmax,*eps,'x',&p);
   }
   freeParameter(&p);
}

void _DLE(int* nb, double* dle,double* x,double* y,double* z,double* tmax,double* eps,double* par) {
   Parameter p;
   initParameter(&p,par);
   for(int i=0;i<*nb;i++) {
     dle[i] = DLE(x[i],y[i],z[i],*tmax,*eps,&p);
   }
   freeParameter(&p);
}

void _DLEmatrix(double* L,double* x,double* y,double* z,double* tmax,double* eps,double* par) {
   Parameter p;
   initParameter(&p,par);
   DLEmatrix(L,*x,*y,*z,*tmax,*eps,&p);
   freeParameter(&p);
}
