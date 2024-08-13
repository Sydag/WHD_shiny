#include "WHD_startingCondition.h"

typedef struct Parameter Parameter;
struct Parameter
{
  double alpha;
  double beta;
  double gamma;
  double delta;
  double eta;
  double ksi;
  double phi;
  double psi;
  double omega;
  double u;
  double v;
  double w;
  double l;
  double zcrit;
  
  double alpha_min;
  double beta_min;
  double gamma_min;
  double delta_min;
  double eta_min;
  double ksi_min;
  double phi_min;
  double psi_min;
  double omega_min;
  double u_min;
  double v_min;
  double w_min;
  double l_min;
  double zcrit_min;
  
  double alpha_max;
  double beta_max;
  double gamma_max;
  double delta_max;
  double eta_max;
  double ksi_max;
  double phi_max;
  double psi_max;
  double omega_max;
  double u_max;
  double v_max;
  double w_max;
  double l_max;
  double zcrit_max;
  
  int alpha_nb;
  int beta_nb;
  int gamma_nb;
  int delta_nb;
  int eta_nb;
  int ksi_nb;
  int phi_nb;
  int psi_nb;
  int omega_nb;
  int u_nb;
  int v_nb;
  int w_nb;
  int l_nb;
  int zcrit_nb;
  
  double alpha_incr;
  double beta_incr;
  double gamma_incr;
  double delta_incr;
  double eta_incr;
  double ksi_incr;
  double phi_incr;
  double psi_incr;
  double omega_incr;
  double u_incr;
  double v_incr;
  double w_incr;
  double l_incr;
  double zcrit_incr;

  //const char par_names[14][20];
  char**    par_names;
  double**  par_values;
  double**  par_min_values;
  double**  par_max_values;
  int**     par_nb_values;
  double**  par_incr_values;
  int*      count;

  int initialized; //will be zero by default (hopefully)
};


//reading parameters
void initializeParameter(Parameter* p);
void freeParameter(Parameter* p);
int  readParameter(Parameter *p,char* filename);
int  copyParameter(Parameter from,Parameter* to);
int  readParameterRange(Parameter* p,char* filename);
int  nextParameter(Parameter* p);
void printVariableParameter(FILE* f,Parameter p);
void printParameterState(FILE* f,Parameter p);
void printParameter(FILE* f,Parameter p);

//functionnal responses (and related functions)
double F(double x, Parameter p);
double G(double y, double z, Parameter p);
double dFdx(double x, Parameter p);
double dGdy(double y, double z, Parameter p);
double dGdz(double y, double z, Parameter p);
double H(double x, void* vp);
double dH(double x, void *vp);

//derivatives of x,y,z
double dxdt(double x, double y,double z, Parameter p);
double dydt(double x, double y,double z, Parameter p);
double dzdt(double x, double y,double z, Parameter p);
double d2xdt2(double x, double y,double z, Parameter p);
double d2ydt2(double x, double y,double z, Parameter p);
double d2zdt2(double x, double y,double z, Parameter p);
int dxyt_dz (double z, const double y[], double f[],void *param);
int dxyz_dt (double t, const double y[], double f[],void *param);

//equilibrium
double   P_gamma_eq(double x, void* p);
double   dP_gamma_eq(double x, void* vp);
double   P_alpha_eq(double x, void* p);
double   dP_alpha_eq(double x, void* vp);
int      NoInfeq(double *yeq,Parameter p);
double** equilibrium(int* nbeq,Parameter p);

int      compute_jacob_dxyz_dt(const double y[], double* dfdy, int complete_system,Parameter p);
int      jacob_dlogxyz_dt (double t, const double y[], double *dfdy,double dfdt[], void *param);

//threshold gamma values
double gamma_clearance(Parameter p);
double gamma_survival(Parameter p);
double gamma_kill(Parameter p);
int    gamma_crit (int maxg,double g_crit[],int curv[],Parameter p);
int    gamma_infection(double x0,double dlogy0,double dlogz0,double* g_crit, double tmax,double maxdist,double eps,int maxiter,double* estim,Parameter p);

//threshold alpha values
double alpha_kill(Parameter p);
int alpha_crit (int maxa,double* a_crit,int* curv,Parameter p);
int alpha_infection(double x0,double dy0,double dz0,double* a_crit, double tmax,double maxdist,double eps,int maxiter,double* estim,Parameter p);
int alpha_spbl(double* estim,Parameter p,double x0,double dlogy0,double dlogz0,double lo_a,double hi_a,double tmax,double eps,int maxiter);
int alpha_blud(double* estim,double target_blud,Parameter p,double tmax,double eps,int maxiter);

//event functions (intended to be called by dyninfWEvent)
double   xdecrease(double x,double y,double z,Parameter p);    //dx/dt
double   ydecrease(double x,double y,double z,Parameter p);    //dy/dt
double   yinflection(double x,double y,double z,Parameter p);  //d2y/dt2
double   zcritreached(double x,double y,double z,Parameter p); //z-zcrit

//Lyapunov exponents / separatrix
int      threshold_x(int nb,double* x,double* lo_x, double *hi_x,double* y,double* z,double tmax,double eps,double max_dist,Parameter p);
int      threshold_y(int nb,double* x,double* y,double* lo_y, double *hi_y,double* z,double tmax,double eps,double max_dist,Parameter p);
int      threshold_z(int nb,double* x,double* y,double* z,double* lo_z, double *hi_z,double tmax,double eps,double max_dist,Parameter p);
double   DLEunivar(double x,double y,double z,double tmax,double eps,char who,Parameter* p);
double   DLE(double x,double y,double z,double tmax,double eps,Parameter* p);
void     DLEmatrix(double* L,double x,double y,double z,double tmax,double eps,Parameter* p);

//various integration functions
int      dyninf_complete(int nb,double* x, double* y, double* z,double* t,Parameter* p);
int      dyninf(double* var,double tmax, int stopWhenDead, int logx, Parameter* p);
int      dyninfWEvent(double* var,double tmax,double* eps,int maxiter,double (*event)(double,double,double,Parameter), int logx, Parameter* p);
int      tsurv(double* var, Parameter* p);

int      spbl(double* var,double** eq,int nbeq,double tmax,double* eps,int maxiter,Parameter p);
double   threshold_spbl_x(double dy,double dz,double tmax,double* eps,int maxiter,Parameter p);

int      blud(double* var,double** eq,int nbeq,double tmax,double* eps,int maxiter,Parameter p);
double   getLT(int* dead,double x,double tmax,double** eq,Parameter* p);
double   avgLT(double tmax,double** eqa,Parameter *pa,StartingCondition* sc);
