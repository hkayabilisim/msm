#ifndef msmLibrary_h
#define msmLibrary_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXDIM 3
#define LMAX 10
#define DIMX 0
#define DIMY 1
#define DIMZ 2
#define MAX_POLY_DEGREE 6
#define MAX_MU 14
#define PI  3.14159265358979323846

#define TAUP4(X) ((35. + (X*X) * (-35. + (21. - 5.*(X*X))*(X*X)) )/16.0)
#define ERFBETA(BETA,X) (X < DBL_EPSILON ? (2.0*BETA/sqrt(M_PI)) : (myerf(BETA*X)/X))

typedef struct LinkedListElement
{
    int    data;
    double value ;
    struct LinkedListElement *next;
    struct LinkedListElement *prev;
} LinkedListElement;

typedef struct LinkedList
{
    struct LinkedListElement *head;
    struct LinkedListElement *tail;
} LinkedList;

typedef struct msm_stats
{
    int    N;
    int    L;
    int    mu;
    double h;
    double hx, hy, hz;
    double a, aL;
    double alpha;
    double beta ;
    int    p;
    int    s;
    double csr;
    int chinx,chiny,chinz,chisize;
    double *cvecx,*cvecy,*cvecz;
    int    clenx,cleny,clenz;
    double chisum,*chi;
    
    int sLsize, sLnx, sLny,sLnz ;
    int    M0[MAXDIM]; /* TODO: to be removed */
    int    Mx, My, Mz ;
    int    MLx, MLy, MLz ;
    int    ML2x, ML2y, ML2z ;
    int    Mlx[LMAX],Mly[LMAX],Mlz[LMAX];
    int    slsizex[LMAX],slsizey[LMAX],slsizez[LMAX];
    int    Ax, Ay, Az;
    int    size[LMAX][MAXDIM];
    double rminx;
    double rminy;
    double rminz;
    double rmaxx;
    double rmaxy;
    double rmaxz;
    double volume;
    double density;
    long int short_range_interactions;
    int    storage;
    
    double netcharge ;
    double netgridcharge[LMAX] ;
    double cputime_initialization;
    double cputime_anterpolation;
    double cputime_restriction;
    double cputime_stencill;
    double cputime_stencilL;
    double cputime_grid_to_grid_mapping;
    double cputime_top_level;
    double cputime_prolongation;
    double cputime_interpolation;
    double cputime_potential_long;
    double cputime_short_range;
    double cputime_self;
    double cputime_total;
    
    double potential_total ;
    double potential_short_total;
    double potential_short_real;
    double potential_short_csr;
    double potential_short_self;
    double potential_long_total;
    double potential_long_real;
    double potential_long_self;
    double potential_long_fourier ;
} msm_stats;


typedef struct MasterBaseFunction
{
    const int p;
    double (*region[MAX_POLY_DEGREE]) (double x);
} MasterBaseFunction;


static double cubic_bspline0 (double t) { return  (1./6) * (2-t) * (2-t) * (2-t);}
static double cubic_bspline1 (double t) { return  (2./3) + t*t*(-1 + 0.5*t);     }
static double cubic_bspline2 (double t) { return  (2./3) + t*t*(-1 - 0.5*t);     }
static double cubic_bspline3 (double t) { return  (1./6) * (2+t) * (2+t) * (2+t);}

static double cubic_bspline0d(double t) { return -0.5 * (2-t) * (2-t);}
static double cubic_bspline1d(double t) { return  t*(-2 + 1.5*t)     ;}
static double cubic_bspline2d(double t) { return  t*(-2 - 1.5*t)     ;}
static double cubic_bspline3d(double t) { return  0.5 * (2+t) * (2+t);}

static double quintic_bspline0 (double t) { return 81./40+t*(-27./8+t*(9./4+t*(-3./4+t*(1./8+t*(-1./120)))));}
static double quintic_bspline1 (double t) { return 17./40+t*(5./8+t*(-7./4+t*(5./4+t*(-3./8+t*(1./24)))));   }
static double quintic_bspline2 (double t) { return 11./20+t*t*(-1./2+t*t*(1./4+t*(-1./12)));}
static double quintic_bspline3 (double t) { return 11./20+t*t*(-1./2+t*t*(1./4+t*(1./12)));}
static double quintic_bspline4 (double t) { return 17./40+t*(-5./8+t*(-7./4+t*(-5./4+t*(-3./8+t*(-1./24)))));}
static double quintic_bspline5 (double t) { return 81./40+t*(27./8+t*(9./4+t*(3./4+t*(1./8+t*(1./120)))));}

static double quintic_bspline0d (double t) { return (-27./8+t*(9./2+t*(-9./4+t*(1./2+t*(-1./24))))) ;}
static double quintic_bspline1d (double t) { return (5./8+t*(-7./2+t*(15./4+t*(-3./2+t*(5./24))))); }
static double quintic_bspline2d (double t) { return (t*(-1+t*t*(1+t*(-5./12))));}
static double quintic_bspline3d (double t) { return (t*(-1+t*t*(1+t*(5./12))));}
static double quintic_bspline4d (double t) { return (-5./8+t*(-7./2+t*(-15./4+t*(-3./2+t*(-5./24))))) ;}
static double quintic_bspline5d (double t) { return (27./8+t*(9./2+t*(9./4+t*(1./2+t*(1./24)))));}


static const struct MasterBaseFunction CubicBSpline  = { 4, {cubic_bspline0,  cubic_bspline1,  cubic_bspline2,  cubic_bspline3} };
static const struct MasterBaseFunction CubicBSplineD = { 4, {cubic_bspline0d, cubic_bspline1d, cubic_bspline2d, cubic_bspline3d} };

static const struct MasterBaseFunction QuinticBSpline  = { 6, {quintic_bspline0, quintic_bspline1,quintic_bspline2,quintic_bspline3,quintic_bspline4,quintic_bspline5}};
static const struct MasterBaseFunction QuinticBSplineD  = { 6, {quintic_bspline0d, quintic_bspline1d,quintic_bspline2d,quintic_bspline3d,quintic_bspline4d,quintic_bspline5d}};


static const double JN[2][MAX_POLY_DEGREE] = {{6./8,  4./8 , 1./8 , 0.   , 0., 0.},
    {20./32, 15./32, 6./32, 1./32, 0., 0.}};


static const double OMEGA[2][MAX_MU] = {
    {3.4641016151377539e+00, /* mu = 0 */
        -1.7320508075688767e+00,
        6.7949192431122685e-01,
        -2.3978297178318450e-01,
        7.9713982076653617e-02,
        -2.5502951436845282e-02,
        7.9437840692459776e-03,
        -2.4260315207973015e-03,
        7.2976833805943851e-04,
        -2.1633348486544432e-04,
        6.1410409127506203e-05, /* mu = 10 */
        -1.4537930745802521e-05, /* mu = 10 + 1 */
        1.9569521248038610e-06, /* mu = 10 + 2 (10 + p/2) */
        0},
    {1.2379121245354655e+01,
        -9.3765192523414775e+00,
        5.8087547854721668e+00,
        -3.2655236531195819e+00,
        1.7352829708690274e+00,
        -8.8893259491018195e-01,
        4.4379214314187693e-01,
        -2.1736815004215146e-01,
        1.0479655626770509e-01,
        -4.9312127870361926e-02,
        2.1608174371344265e-02,
        -7.8774221940891460e-03,
        1.9814890638090049e-03,
        -2.4354138541380615e-04}
};


msm_stats msmperiodic(int N,double *r,double *q,int p,double abar,int mu,double *uappr,double *fappr,
                      double Ax, double Ay, double Az);
double choosebeta(double r);
double choosebetaNew(double epsilon,double aL,double hL);
double choosekmax(double epsilon,double beta,double hL);

int nist_read_water(FILE *fp,double **q,double **r,int *N,double *Lx,double *Ly,double *Lz);
int nist_read_generic(FILE *fp,double **q,double **r,int *N,double *Lx,double *Ly,double *Lz);
int msm4g_math_factorial(int n);
int msm4g_math_cantor(int xx,int yy,int zz);
int msm4g_linkedlist_size(LinkedList *list);
LinkedList *msm4g_linkedlist_new();
void msm4g_linkedlist_add(LinkedList *list,int data);
int msm4g_linkedlist_search(LinkedList *list,int data);
void msm4g_linkedlist_addvalue(LinkedList *list,int data,double value);
double msm4g_linkedlist_getvalue(LinkedList *list,int index);
void msm4g_linkedlist_destroy(LinkedList *list);

double norm(int n, double *x);
double diffnorm(int n, double *x, double *y);
double gammarho(double p,int s);
double gammarhod(double p,int s);
double msm_gammal(int kx,int ky,int kz,double h,double a,int s);
double msm_gammaL(int kx,int ky,int kz,double h,double a,int s);
void direct(int N,double *r,double *q,double *ureal,double *freal);
void short_range(int i,int j, double *r,double *q,int s,double a,double *fshort,double *ushort);
double short_range_periodic(int i,int j, double *r,double *q,int s,double a,double *fshort,double Ax,double Ay,double Az);
void print_stats(msm_stats stats);
msm_stats msm3d(int N,double *r,double *q,int p,double h,double a,int mu,int s,double *uappr,double *fappr);
void testcube();
void testcube_ewald();
int teststep();
void scaling();
void g2gtest();
int testperiodic();
void msm_scalability();
void testperiodic_nist();
void report20171108();
void report20171114();
void report20171122();
void report20171122_NaCl();
void report20171129_ewald();
void report20171129_direct();
void report20171129();
void report20171221_changa();
void testperiodic_NaCl();
void test_nist(const char argv[],int isWater);
void testewald_nist();
void testewald();
void testMadelungSum();
int chooseM(double htilde,double Ax);
int minOf3Integers(int a, int b, int c);
int maxOf3Integers(int a, int b, int c);

int chooseL(int Mx,int My,int Mz,int N);
void report20171221_ewald();

int data_read(char *filename, double **q, double **r, int *N, double *Lx, double *Ly, double *Lz);

float myerfc(float x);
#endif /* msmLibrary_h */
