#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>

#define X 0
#define Y 1
#define Z 2
#define NMAX  1024
#define MMAX  10
#define MYPI  3.141592653589793238462643


int N ;
int Mx, My, Mz ;
int v;
int mu;
int kmax;
int pmax;
double Ax, Ay, Az, detA ;
double abar ;
double a,aL;
double h,hL,hx,hy,hz;
double beta;

double ushort_real;
double ushort_self;
double ushort_csr ;
double ulong_real ; double ulong_real_expected;
double ulong_four ; 
double ulong_four_expected;
double ulong_four_expected_cossum;
double ulong_four_expected_sinsum;
double ulong_self ;
double utotal_expected;


double r[NMAX][3] ;
double q[NMAX];
double qm[MMAX][MMAX][MMAX];
double em[MMAX][MMAX][MMAX];
double em_fourier[MMAX][MMAX][MMAX];


double Bspline(int k,double u) {
    if (k==1) {
        if (u >= 0 && u < 1)
            return 1.0;
        else
            return 0.0;
    } else {
        return u/(double)(k-1) * Bspline(k-1,u) + (k-u)/(double)(k-1) * Bspline(k-1,u-1.0);
    }
}


double gama(double rho) {
    double out = 0.0;
    if (rho >= 1.0)
        out = 1.0 / rho ;
    else {
        double  rho2 = rho * rho ;
        for (int k = 0 ; k < v ; k++) {
            double binom = tgamma(0.5)/(tgamma(k+1)*tgamma(0.5-k)) ;
            out += binom * pow(rho2-1.0 , k);
        }
    }
    return out;
}

double g1(double r) {
    return (1.0 / a) * gama(r/a) ;
}
double g0(double r) {
    return 1.0 / r - g1(r) ;
}
double gLR(double r) {
    if (r < DBL_EPSILON)
        return 2.0 * beta / sqrt(MYPI) ;
    else
        return erf(beta * r) / r ;
}
double g1star(double r) {
    return g1(r)  - gLR(r) ;
}

void load_benchmark_NaNaN8() {
    N    = 8;
    Ax   = Ay = Az = 2;
    detA = Ax * Ay * Az ;
    v    = 4;
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 1;
    Mx   = My = Mz = 2;
    a = abar * h ;
    aL = 2 * a ;
    beta = 3.9952779973449493e-01;
    r[0][X] = 0 ; r[0][Y] = 0 ; r[0][Z] = 0; q[0] = 1;
    r[1][X] = 0 ; r[1][Y] = 0 ; r[1][Z] = 1; q[1] = 1;
    r[2][X] = 0 ; r[2][Y] = 1 ; r[2][Z] = 0; q[2] = 1;
    r[3][X] = 1 ; r[3][Y] = 0 ; r[3][Z] = 0; q[3] = 1;
    r[4][X] = 1 ; r[4][Y] = 1 ; r[4][Z] = 0; q[4] = 1;
    r[5][X] = 1 ; r[5][Y] = 0 ; r[5][Z] = 1; q[5] = 1;
    r[6][X] = 0 ; r[6][Y] = 1 ; r[6][Z] = 1; q[6] = 1;
    r[7][X] = 1 ; r[7][Y] = 1 ; r[7][Z] = 1; q[7] = 1;
}


void calculate_ushort_real() {
    ushort_real = 0.0 ;
    for (int i = 0 ; i < N ; i++) {
        for (int j = 0 ; j < N ; j++) {
            if (i == j) continue;
            for (int px = -pmax ; px <= pmax ; px++) {
                for (int py = -pmax ; py <= pmax ; py++) {
                    for (int pz = -pmax ; pz <= pmax ; pz++) {
                        double rx = r[i][X] - r[j][X] - Ax * px ;
                        double ry = r[i][Y] - r[j][Y] - Ay * py ;
                        double rz = r[i][Z] - r[j][Z] - Az * pz ;
                        double rlen2 = rx * rx + ry * ry + rz * rz ;
                        double rlen = sqrt(rlen2);
                        ushort_real += 0.5 * q[i] * q[j] * g0(rlen);
                    }
                }
            }
        }
    }
}

void calculate_ushort_self() {
    ushort_self = 0.0;
    for (int i = 0 ; i < N ; i++) {
        for (int px = -pmax ; px <= pmax ; px++) {
            for (int py = -pmax ; py <= pmax ; py++) {
                for (int pz = -pmax ; pz <= pmax ; pz++) {
                    if (px == 0 && py == 0 && pz == 0) continue;
                    double rlen2 = Ax * px * Ax * px + 
                                   Ay * py * Ay * py +
                                   Az * pz * Az * pz ;
                    double rlen = sqrt(rlen2);
                    ushort_self += 0.5 * q[i] * q[i] * g0(rlen);
                }
            }
        }
    }
}

void calculate_ushort_csr() {
    double csr = MYPI / (beta * beta * detA );
    double qsum = 0.0;
    for (int i = 0 ; i < N ; i++ ) 
        qsum += q[i] ;
    ushort_csr = - 0.5 * qsum * qsum * csr ;
}

void calculate_ulong_real_expected() {
    ulong_real_expected = 0.0;
    for (int i = 0 ; i < N ; i++) {
        for (int j = 0 ; j < N ; j++) {
            for (int px = -pmax ; px <= pmax ; px++) {
                for (int py = -pmax ; py <= pmax ; py++) {
                    for (int pz = -pmax ; pz <= pmax ; pz++) {
                        double rx = r[i][X] - r[j][X] - Ax * px ;
                        double ry = r[i][Y] - r[j][Y] - Ay * py ;
                        double rz = r[i][Z] - r[j][Z] - Az * pz ;
                        double rlen2 = rx * rx + ry * ry + rz * rz ;
                        double rlen = sqrt(rlen2);
                        ulong_real_expected += 0.5 * q[i] * q[j] * g1star(rlen);
                    }
                }
            }
        }
    }
}

void calculate_ulong_four_expected() {
    ulong_four_expected = 0.0;
    ulong_four_expected_cossum = 0.0;
    ulong_four_expected_sinsum = 0.0;
    for (int i = 0 ; i < N ; i++) {
        for (int j = 0 ; j < N ; j++) {
            for (int kx = -kmax ; kx <= kmax ; kx++) {
                for (int ky = -kmax ; ky <= kmax ; ky++) {
                    for (int kz = -kmax ; kz <= kmax ; kz++) {
                        if (kx == 0 && ky == 0 && kz == 0)
                            continue;
                        double kvecx = kx / Ax ;
                        double kvecy = ky / Ay ;
                        double kvecz = kz / Az ;
                        double k2 = kvecx * kvecx + kvecy * kvecy + kvecz * kvecz ;
                        double chi = (1.0/(MYPI * k2 * detA)) * exp(-MYPI * MYPI * k2 / (beta * beta));
                        double dotprod = kvecx * (r[i][X] - r[j][X]) +
                                         kvecy * (r[i][Y] - r[j][Y]) +
                                         kvecz * (r[i][Z] - r[j][Z]) ;
                        ulong_four_expected_cossum += 0.5 * q[i] * q[j] * chi * cos(2.0 * MYPI * dotprod) ;
                        ulong_four_expected_sinsum += 0.5 * q[i] * q[j] * chi * sin(2.0 * MYPI * dotprod) ;
                    }
                }
            }
        }
    }
    ulong_four_expected = ulong_four_expected_cossum ;
}

void calculate_ulong_self() {
    double q2sum = 0.0;
    for (int i = 0 ; i < N ; i++)
        q2sum += q[i] * q[i] ;
    ulong_self = - 0.5 * q2sum * g1(0) ;
}

void do_anterpolation() {
    for (int mx = 0 ; mx < Mx ; mx++) {
        for (int my = 0 ; my < My ; my++) {
            for (int mz = 0 ; mz < Mz ; mz++) {
                qm[mx][my][mz] = 0.0;
                for (int px = -pmax ; px <= pmax ; px++) {
                    for (int py = -pmax ; py <= pmax ; py++) {
                        for (int pz = -pmax ; pz <= pmax ; pz++) {
                            for (int i = 0 ; i < N ; i++) {
                                double sx = r[i][X] / Ax ;
                                double sy = r[i][Y] / Ay ;
                                double sz = r[i][Z] / Az ;
                                double phix = Bspline(v, Mx * (sx - px) - mx + v/2);
                                double phiy = Bspline(v, My * (sy - py) - my + v/2);
                                double phiz = Bspline(v, Mz * (sz - pz) - mz + v/2);
                                qm[mx][my][mz] += phix * phiy * phiz * q[i];
                            }
                        }
                    }
                }
            }
        }
    }
}

void display_results(){
     printf("%-28s : %5.2f %5.2f %5.2f\n","Ax Ay Az" ,Ax,Ay,Az);
     printf("%-28s : %25.16e\n","ushort_real" ,ushort_real);
     printf("%-28s : %25.16e\n","ushort_self" ,ushort_self);
     printf("%-28s : %25.16e\n","ushort_csr"  ,ushort_csr);
     printf("%-28s : %25.16e\n","ulong_self"  ,ulong_self);
     printf("%-28s : %25.16e\n","ulong_real_expected"  ,ulong_real_expected);
     printf("%-28s : %25.16e\n","ulong_four_expected_cossum"  ,ulong_four_expected_cossum);
     printf("%-28s : %25.16e\n","ulong_four_expected_sinsum"  ,ulong_four_expected_sinsum);
     printf("%-28s : %25.16e\n","utotal_expected",utotal_expected);

}
void display_grid_charges() {
    for (int mx = 0 ; mx < Mx ; mx++) {
        for (int my = 0 ; my < My ; my++) {
            for (int mz = 0 ; mz < Mz ; mz++) {
                printf("qm[%d][%d][%d] = %10.5f\n",mx,my,mz,qm[mx][my][mz]);
            }
        }
    }
}

int main(int argc, char *argv[]) {
  
    load_benchmark_NaNaN8();
    
    calculate_ushort_real();   
    calculate_ushort_self();   
    calculate_ushort_csr();   
    
    calculate_ulong_self();   
    calculate_ulong_real_expected();
    calculate_ulong_four_expected();

    do_anterpolation();
    display_grid_charges();
    
    utotal_expected = ushort_real + ushort_self + ushort_csr + ulong_real_expected +
                      ulong_four_expected + ulong_self ;
    display_results();
    
    return 0;
}

