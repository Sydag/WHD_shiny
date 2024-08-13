#ifndef STARTINGCONDITION
#define STARTINGCONDITION

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

typedef struct StartingCondition StartingCondition ;
struct StartingCondition 
{
  double mu_x;
  double mu_y;
  double mu_z;
  
  double mu_logx;
  double mu_logy;
  double mu_logz;
  
  double sigma_logx;
  double sigma_logy;
  double sigma_logz;
   
  double loop_x_min;
  double loop_y_min;
  double loop_z_min;

  double loop_x_max;
  double loop_y_max;
  double loop_z_max;

  int loop_x_nb;
  int loop_y_nb;
  int loop_z_nb;

  double loop_x_incr;
  double loop_y_incr;
  double loop_z_incr;

  double sigma_logx_min;
  double sigma_logy_min;
  double sigma_logz_min;
  
  double sigma_logx_max;
  double sigma_logy_max;
  double sigma_logz_max;
  
  int sigma_logx_nb;
  int sigma_logy_nb;
  int sigma_logz_nb;
    
  double sigma_logx_incr;
  double sigma_logy_incr;
  double sigma_logz_incr;

  int logscale[3]; //boolean : 1 if log scale, zero otherwise
  
  int count[6];
};

int readStartingCondition(StartingCondition *sc,char* filename);
void printStartingX(StartingCondition sc);
void printStartingY(StartingCondition sc);
void printStartingZ(StartingCondition sc);
void printStarting(StartingCondition sc);
int nextStartingCondition(StartingCondition *sc);
int setLoopXYZ(char who,double min,double max,int nb,int logscale,StartingCondition *sc);
int setLoopSigmaXYZ(int who_s,double min_s,double max_s,int nb_s,StartingCondition *sc);
double getXini_gsl(StartingCondition sc,gsl_rng * r);
double getYini_gsl(StartingCondition sc,gsl_rng * r);
double getZini_gsl(StartingCondition sc,gsl_rng * r);

#endif
