#ifdef __APPLE__
#include <OpenCL/opencl.h>
#include <unistd.h>
#else
#include <CL/cl.h>
#endif



#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>


#include "err_code.h"
#include "device_picker.h"
#define LOW -1.0e15
#define MAXERR 1e-6
#define REARTH 6378.388
#define MAX_BINARY_SIZE (0x1000000)

#define MAX_SOURCE_SIZE (0x100000)
#define BIN_PATH "Kernel.clbin"

#define SEP printf("-----------------------------------------------------------\n")



//---------START GLOBAL VARIABLES-----------


//int *first;//vector of index in the bivariate case

//---------END GLOBAL VARIABLES-------------


//---------START DECLARING FUNCTIONS-----------
void DoubleExpOCL(int *npts, int *ntime,double *coordt, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *data,int *weigthed, double *mom_cond, int *dist, double *coordx,double *coordy,double *gradcor,double  *grad, double *ww,int *local_wi, int *dev);

void GneitingOCL(int *npts, int *ntime,double *coordt, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *data,int *weigthed, double *mom_cond, int *dist, double *coordx,double *coordy,double *gradcor,double  *grad, double *ww,int *local_wi, int *dev);

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
void SubSamp_space_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *a, double *b,double *block_mean,int *weigthed, int *local_wi, int *dev);
void scalar_time(int *ncoord,int *nstime,double *sublagt,int *cormod,double *parcor,int *flagcor, double *gradcor,int *flagnuis, double *grad,int *npar,double *nuis, double *sdata,int *weigthed,double *maxtime, double *ww, double *mom_cond,int *dist, double *coordx, double *coordy,double *maxdist);
void SubSamp_time(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *a, double *b,double *winc, double *winstp,double *block_mean,int *weigthed, int *local_wi, int *dev);
void SubSamp_time_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *a, double *b,double *winc, double *winstp,double *block_mean,int *weigthed, int *local_wi, int *dev);
void scalar_spacetime(int *npts,int *nstime,double *sublagt,double *maxtime,int *cormod,double *parcor,int *flagcor, double *gradcor,int *flagnuis, double *grad,int *npar,double *nuis,double *s2data,int *weigthed, double *ww, double *mom_cond,double *maxdist,int *dist, double *scoordx, double *scoordy);
void SubSamp_spacetime(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *winc_t, double *winstp_t,double *block_mean,int *weigthed, int *local_wi, int *dev);
void SubSamp_spacetime_ocl(double *coordx, double *coordy,double *coordt, int *ncoord,int *ntime,int *cormod,double *data,int *dist, double *maxdist,double *maxtime,int *npar,double *parcor,int *nparc,double *nuis,int *nparnuis, int *flagcor, int *flagnuis,double *vari,double *winc, double *winstp,double *winc_t, double *winstp_t,double *block_mean,int *weigthed, int *local_wi, int *dev);
char * getKernelSource(char *filename);
float sum_total(float *arr, int ngrid);
void DoubleExpOCL(int *npts, int *ntime,double *coordt, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *data,int *weigthed, double *mom_cond, int *dist, double *coordx,double *coordy,double *gradcor,double  *grad, double *ww,int *local_wi, int *dev);
void GneitingOCL(int *npts, int *ntime,double *coordt, double *maxtime,double *maxdist,int *cormod, double *parcor,int *flagcor,int *flagnuis,int *npar,double *nuis,double *data,int *weigthed, double *mom_cond, int *dist, double *coordx,double *coordy,double *gradcor,double  *grad, double *ww,int *local_wi, int *dev);
int DevOpenCL();
void create_binary_kernel(int *dev, char **fname);




//---------END DECLARING FUNCTIONS-----------
