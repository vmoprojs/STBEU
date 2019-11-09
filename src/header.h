#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <assert.h>
#include <stdint.h>

//#include "err_code.h"
//#include "device_picker.h"
#define LOW -1.0e15
#define MAXERR 1e-6
#define REARTH 6378.388
#define MAX_BINARY_SIZE (0x1000000)

#define MAX_SOURCE_SIZE (0x100000)
//#define BIN_PATH "Kernel.clbin"

#define EPS1 1.0e-5
#define SQE 3.162278e-30
//#define SQE 1e-15

#define SEP Rprintf("-----------------------------------------------------------\n")



//---------START GLOBAL VARIABLES-----------


//int *first;//vector of index in the bivariate case

//---------END GLOBAL VARIABLES-------------


//---------START DECLARING FUNCTIONS-----------

void Grad_Pair_Gauss(double rho, int *flag, double *gradcor, double *grad,
                     int *npar, double *par, double u, double v);
double DGneiting_sep(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_pw_s(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_pw_t(double h,double u, double power_s,double power_t,double scale_s,double scale_t,double sep);
double DGneiting_sc_t(double h,double u, double power_s,double power_t,
                      double scale_s,double scale_t,double sep);
double DGneiting_sc_s(double h,double u,double power_s,double power_t,
                      double scale_s,double scale_t,double sep);
double DStabSc(double lag, double power, double scale, double rho);
double CorFunBohman(double lag,double scale);
double CorFunStable(double lag, double power, double scale);
double CorFct(int *cormod, double h, double u, double *par);
void GradCorrFct(double rho, int *cormod, int *flag,
                 double *grad, double h, double u,double *par);
double Dist_geodesic(double loni, double lati, double lonj, double latj);
void Range(double *x, double *ran, int size);
int is_equal(double val1, double val2);
void SeqStep(double *x, int len, double step, double *res);
void SetSampling_s(int ncoord,int ntime,double *coordx, double *coordy, double *data, int *npts,double *scoordx, double *scoordy, double *sdata, double xmax,double xmin, double ymax,double ymin);
void SetSampling_t(double *data,double *sdata,int ncoord,int ntime,int wint,int k);
void scalar_space(int *npts, int *ntime,double *coordt, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *sdata,int *weigthed, double *mom_cond, int *dist, double *scoordx,double *scoordy,double *gradcor,double  *grad, double *ww );
void SubSamp_space(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *a, double *b,double *block_mean,int *weigthed, int *local_wi, int *dev);

void scalar_time(int *ncoord,int *nstime,double *sublagt,int *cormod,double *parcor,int *flagcor, double *gradcor,int *flagnuis, double *grad,int *npar,double *nuis, double *sdata,int *weigthed,double *maxtime, double *ww, double *mom_cond,int *dist, double *coordx, double *coordy,double *maxdist);
void SubSamp_time(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *a, double *b,double *winc, double *winstp,double *block_mean,int *weigthed, int *local_wi, int *dev);

void scalar_spacetime(int *npts,int *nstime,double *sublagt,double *maxtime,int *cormod,double *parcor,int *flagcor, double *gradcor,int *flagnuis, double *grad,int *npar,double *nuis,double *s2data,int *weigthed, double *ww, double *mom_cond,double *maxdist,int *dist, double *scoordx, double *scoordy);
void SubSamp_spacetime(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *winc_t, double *winstp_t,double *block_mean,int *weigthed, int *local_wi, int *dev);







//---------END DECLARING FUNCTIONS-----------





//---------START WENDLAND FUNCTIONS-----------

/* START Wendland covariance */

/* integrand  in  generalized wendland function*/
double int_gen(double x,double mu, double alpha,double lag,double supp);
void integr_gen(double *x, int n, void *ex);
// function computing generalized wendland
double wendintegral(double x, double *param);
/* generalized wendland function*/
double CorFunW_gen(double lag,double R_power1,double smooth,double scale);  // mu alpha beta
double wen_time(double *par, double h,double u);

double RES_CorFunW_gen(double *lag,double *R_power1,double *smooth,double *scale, double *res);
double RES_wen_time(double *par, double *h,double *u, double *res);
/* END Wendland covariance */


/* START DERIVATIVES Wendland covariance */

// SCALE_S:
double deri_scale_s_wen_time(double *par, double h,double u);
double RES_deri_scale_s_wen_time(double *par, double *h,double *u, double *res);
// SCALE_T:
double deri_scale_t_wen_time(double *par, double h,double u);
double RES_deri_scale_t_wen_time(double *par, double *h,double *u, double *res);
// SMOOTH:
double deri_smooth_wen_time(double *par, double h,double u);
double RES_deri_smooth_wen_time(double *par, double *h,double *u, double *res);
// SILL:
double deri_sill_wen_time(double *par, double h,double u);
double RES_deri_sill_wen_time(double *par, double *h,double *u, double *res);
// SEP:
double deri_sep_wen_time(double *par, double h,double u);
double RES_deri_sep_wen_time(double *par, double *h,double *u, double *res);
// R_power:
double deri_R_power_wen_time(double *par, double h,double u);
double RES_deri_R_power_wen_time(double *par, double *h,double *u, double *res);
// R_power_t:
double deri_R_power_t_wen_time(double *par, double h,double u);
double RES_deri_R_power_t_wen_time(double *par, double *h,double *u, double *res);


/* END DERIVATIVES Wendland covariance */


//---------END WENDLAND FUNCTIONS-----------












/*----------------------------------------------------------------
 File name: Host.c
 Description: procedures for OCL computation.
 Start
 ---------------------------------------------------------------*/
#define CL_SILENCE_DEPRECATION
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <unistd.h>
#else
#include <CL/cl.h>
#endif


#pragma once

#include <string.h>
//#include "err_code.h"

#define MAX_PLATFORMS     8
#define MAX_DEVICES      16
#define MAX_INFO_STRING 256

#ifdef __cplusplus
#include <cstdio>
#endif

const char *err_code (cl_int err_in);
void check_error(cl_int err, const char *operation, char *filename, int line);
#define checkError(E, S) check_error(E,S,__FILE__,__LINE__)
unsigned getDeviceList(cl_device_id devices[MAX_DEVICES]);
void getDeviceName(cl_device_id device, char name[MAX_INFO_STRING]);
int parseUInt(const char *str, cl_uint *output);
void parseArguments(int argc, char *argv[], cl_uint *deviceIndex);
void param_OCL(int *npts,int *ntime,int *flagcor,int *flagnuis,int *npar,int *weigthed, int *dist,double *maxtime,double *maxdist,double *parcor,double *nuis,int *nparc,int *cormod,int *int_par, double *dou_par,int *nparnuis);

void exec_kernel(int *npts, int *ntime,double *coordt, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *data,int *weigthed, double *mom_cond, int *dist, double *coordx,double *coordy,double *gradcor,double  *grad, double *ww,int *local_wi, int *dev,char **fname,int *nparc,int *nparnuis);

char * getKernelSource(char *filename);
double sum_total(double *arr, int ngrid);

int DevOpenCL();
void create_binary_kernel(int *dev, char **fname);
void SubSamp_space_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *a, double *b,double *block_mean,int *weigthed, int *local_wi, int *dev,char **fname);
void SubSamp_time_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *a, double *b,double *winc, double *winstp,double *block_mean,int *weigthed,int *local_wi, int *dev,char **fname);
void SubSamp_spacetime_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *winc_t, double *winstp_t,double *block_mean,int *weigthed, int *local_wi, int *dev,char **fname);
void create_binary_kernel_GPU(int *dev, char **fname);

/*----------------------------------------------------------------
 File name: Host.c
 Description: procedures for OCL computation.
 End
 ---------------------------------------------------------------*/
