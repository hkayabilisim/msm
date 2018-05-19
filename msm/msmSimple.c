#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>

#define X 0
#define Y 1
#define Z 2
#define NMAX  1000
#define MMAX  100
#define SMAX  100
#define MYPI  3.141592653589793238462643

// Pre-preprocessed for v=4 and mu=10
double wprime[13] = {
     3.4641016151377539e+00, // 0
    -1.7320508075688767e+00,
     6.7949192431122685e-01,
    -2.3978297178318450e-01,
     7.9713982076653617e-02,
    -2.5502951436845282e-02,
     7.9437840692459776e-03,
    -2.4260315207973015e-03,
     7.2976833805943851e-04,
    -2.1633348486544432e-04,
     6.1410409127506203e-05,  // 10
    -1.4537930745802521e-05,
     1.9569521248038610e-06}; // 10 + v/2 (v=4)

int N ;                    // # particles
int Mx, My, Mz ;           // # grid points
int v;                     // B-spline (v=4)
int mu;                    // Quasi-interpolation (mu=10) 
int kmax;                  // # wavenumbers in Fourier sum
int pmax;                  // # replications in direct-sum

char dataname[60];         // benchmark name 
double Ax, Ay, Az;         // Width of periodic cell
double detA;               // volume (Ax * Ay * Az)
double abar ;              // relative cuttoff
double a,aL;               // cutoffs
double h,hL,hx,hy,hz;      // grid spacings
double beta;               // Ewald splitting parameter

double ushort_real;     // Eq.11 first component
double ushort_self;     // Eq.11 second component
double ushort_csr ;     // Eq.11 third component     
double ulong_real ;     // Eq.12 first component when l=1
double ulong_real_expected; // Eq.14 first component when l=1 
double ulong_four ;               // Eq.12 first component when l=2
double ulong_four_expected;       // Eq.14 first component when l=2 
double ulong_four_expected_cossum;// cos sum of Eq.14 first component when l=2 
double ulong_four_expected_sinsum;// sin sum of Eq.14 first component when l=2 
double ulong_self ;               // Eq.12 second component
double utotal;                    // Eq.11 + Eq.14
double utotal_expected;           // Eq.11 + Eq.12

double r[NMAX][3];                   // Particle positions (0 <= r <= Ax)
double q[NMAX];                      // Particle charges
double qm[MMAX][MMAX][MMAX];         // Grid charges
double em[MMAX][MMAX][MMAX];         // Grid potentials for l=1
double em_fourier[MMAX][MMAX][MMAX]; // Grid potentials for l=2
double Kl[SMAX][SMAX][SMAX];         // Stencil for l=1
double KL[SMAX][SMAX][SMAX];         // Stencil for l=2

double Bspline(int k,double u) ; // Recursive definition in 2.2.1
double gama(double rho);         // Eq.31 and the one before
double gLR(double r);            // Eq.1 and 2
double g0(double r);             // Eq.15 first part
double g1(double r);             // Eq.15 last part

void calculate_ushort_real();
void calculate_ushort_self();
void calculate_ushort_csr();
void calculate_ulong_self();
void calculate_ulong_real_expected();
void calculate_ulong_four_expected();
void do_anterpolation();
void calculate_ulong_real();
void calculate_ulong_four();

void load_benchmark_NaNaN8();
void load_benchmark_NaNaN64();
void load_benchmark_NaNaN512();
void load_benchmark_NaClN8();
void load_benchmark_NaClN64();
void load_benchmark_NaClN512();
void load_benchmark_CsClN2();
void load_benchmark_CsClN16();
void load_benchmark_CsClN128();
void load_benchmark_changaN8();
void load_benchmark_changaN64();
void load_benchmark_changaN512();
void run_benchmark(int id);

int main(int argc, char *argv[]) {    
    for (int i = 1 ; i <= 12 ; i++ )
        run_benchmark(i);
    return 0;
}

// Appendix A: Eq. 31 and the one before.
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

// Eq.1 and 2 --> gLR(r) = erf(beta * r)/r
double gLR(double r) {
    if (r < DBL_EPSILON)
        return 2.0 * beta / sqrt(MYPI) ;
    else
        return erf(beta * r) / r ;
}

// Eq. 5 first part
double g0(double r) {
    return 1/r - (1/a) * gama(r/a)  ;
}

// Eq. 15 last part when L=1
double g1(double r) {
    return (1/a) * gama(r/a) - gLR(r) ;
}

// Eq. 11 first component
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

// Eq. 11 second component
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

// Eq. 11 third component
void calculate_ushort_csr() {
    double csr = MYPI / (beta * beta * detA );
    double qsum = 0.0;
    for (int i = 0 ; i < N ; i++ )
        qsum += q[i] ;
    ushort_csr = - 0.5 * qsum * qsum * csr ;
}

// Eq. 12 first component for l = 1
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
                        ulong_real_expected += 0.5 * q[i] * q[j] * g1(rlen);
                    }
                }
            }
        }
    }
}

// Eq. 12 first component for l = 2 
// coupled with Eq. 4 and 5.
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
                        double kvecx = kx / Ax ; // Eq. 5 (kvec = A^{-1}k)
                        double kvecy = ky / Ay ;
                        double kvecz = kz / Az ;                        
                        double k2 = kvecx * kvecx + kvecy * kvecy + kvecz * kvecz ;
                        // Eq. 5
                        double chi = (1.0/(MYPI * k2 * detA)) * exp(-MYPI * MYPI * k2 / (beta * beta));
                        double dotprod = kvecx * (r[i][X] - r[j][X]) +
                                         kvecy * (r[i][Y] - r[j][Y]) +
                                         kvecz * (r[i][Z] - r[j][Z]) ;
                        // cos and sin components in Eq. 4
                        ulong_four_expected_cossum += 0.5 * q[i] * q[j] * chi * cos(2.0 * MYPI * dotprod) ;
                        ulong_four_expected_sinsum += 0.5 * q[i] * q[j] * chi * sin(2.0 * MYPI * dotprod) ;
                    }
                }
            }
        }
    }
    // sin sum is practically zero so I take only cos sum.
    ulong_four_expected = ulong_four_expected_cossum ;
}

// Eq. 12 second component
// Please note that upper limit for l-sum is L+1
// Then the sum boils down to (1/a) * gamma(0)
void calculate_ulong_self() {
    double q2sum = 0.0;
    for (int i = 0 ; i < N ; i++)
        q2sum += q[i] * q[i] ;
    ulong_self = - 0.5 * q2sum * (1/a) * gama(0) ;
}

// The equation just after Eq. 14
// The range for m obeys Eq. 23
void do_anterpolation() {
    for (int mx = - Mx/2 + 1 ; mx <= Mx/2 ; mx++) {
        for (int my = - My/2 + 1 ; my <= My/2 ; my++) {
            for (int mz = - Mz/2 + 1 ; mz <= Mz/2 ; mz++) {
                double sum = 0.0 ;
                // p-sum in Eq.22
                for (int px = -pmax ; px <= pmax ; px++) {
                    for (int py = -pmax ; py <= pmax ; py++) {
                        for (int pz = -pmax ; pz <= pmax ; pz++) {
                            for (int i = 0 ; i < N ; i++) {
                            	// s = inv(A) * r
                                double sx = r[i][X] / Ax ;
                                double sy = r[i][Y] / Ay ;
                                double sz = r[i][Z] / Az ;
                                double phix = Bspline(v, Mx * (sx - px) - mx + v/2);
                                double phiy = Bspline(v, My * (sy - py) - my + v/2);
                                double phiz = Bspline(v, Mz * (sz - pz) - mz + v/2);
                                sum += phix * phiy * phiz * q[i];
                            }
                        }
                    }
                }
                qm[mx + Mx/2 - 1][my + My/2 - 1][mz + Mz/2 - 1] = sum ;
            }
        }
    }
}

// Section 2.2.3
// The range for m obeys Eq. 23
void calculate_stencil_Kl() {
    for (int mx = - Mx/2 + 1 ; mx <= Mx/2 ; mx++) {
        for (int my = - My/2 + 1 ; my <= My/2 ; my++) {
            for (int mz = - Mz/2 + 1 ; mz <= Mz/2 ; mz++) {                
                double sum = 0.0 ;                
                for (int kx = - mu - v/2 ; kx <= mu + v/2 ; kx++) {
                    for (int ky = - mu - v/2 ; ky <= mu + v/2 ; ky++) {
                        for (int kz = - mu - v/2 ; kz <= mu + v/2 ; kz++) {
                        	// omega prime is pre-precalculated
                            double w = wprime[abs(kx)] * wprime[abs(ky)] * wprime[abs(kz)] ;
                            // The sum in Eq. 24
                            for (int px = -pmax ; px <= pmax ; px++) {
                                for (int py = -pmax ; py <= pmax ; py++ ) {
                                    for (int pz = -pmax ; pz <= pmax ; pz++) {
                                    	// Hl in Eq. 24 is 1/M because of Eq. 27
                                        double rx = Ax * ((mx + kx)/(double)Mx - px) ;
                                        double ry = Ay * ((my + ky)/(double)My - py) ;
                                        double rz = Az * ((mz + kz)/(double)Mz - pz) ;
                                        double rlen2 = rx * rx + ry * ry + rz * rz ;
                                        double rlen = sqrt(rlen2);
                                        sum +=  w * g1(rlen);
                                    }
                                }
                            }
                        }
                    }
                }
                Kl[mx + Mx/2 - 1][my + My/2 - 1][mz + Mz/2 - 1]  = sum ;
            }
        }
    }
}

// The grid potential e_m calculated for l=1
// There is no specific equation for this in the article but it is:
// e_m = \sum_n Kl_{m-n} q_n 
// The range for m and n obeys Eq. 23
void do_grid_to_grid_mapping() {
    for (int mx = - Mx/2 + 1 ; mx <= Mx/2 ; mx++) {
        for (int my = - My/2 + 1 ; my <= My/2 ; my++) {
            for (int mz = - Mz/2 + 1 ; mz <= Mz/2 ; mz++) {
                
                double sum = 0.0 ;                
                for (int nx = - Mx/2 + 1 ; nx <= Mx/2 ; nx++) {
                    for (int ny = - My/2 + 1 ; ny <= My/2 ; ny++) {
                        for (int nz = - Mz/2 + 1 ; nz <= Mz/2 ; nz++) {
                            
                        	// Kl(m-n) is needed 
                            int mnx = mx - nx ;
                            int mny = my - ny ;
                            int mnz = mz - nz ;
                            
                            // Kl is periodic 
                            if (mnx < - Mx/2 + 1) { mnx += Mx ; } while (mnx < - Mx/2 + 1) ;
                            if (mnx >   Mx/2    ) { mnx -= Mx ; } while (mnx >   Mx/2    ) ;
                            
                            if (mny < - My/2 + 1) { mny += My ; } while (mny < - My/2 + 1) ;
                            if (mny >   My/2    ) { mny -= My ; } while (mny >   My/2    ) ;
                            
                            if (mnz < - Mz/2 + 1) { mnz += Mz ; } while (mnz < - Mz/2 + 1) ;
                            if (mnz >   Mz/2    ) { mnz -= Mz ; } while (mnz >   Mz/2    ) ;
                            
                            sum += Kl[mnx + Mx/2 - 1][mny + My/2 - 1][mnz + Mz/2 - 1] *
                                   qm[nx  + Mx/2 - 1][ny  + My/2 - 1][nz  + Mz/2 - 1] ;                           
                        }
                    }
                }                
                em[mx + Mx/2 - 1][my + My/2 - 1][mz + Mz/2 - 1] = sum ;
            }
        }
    }
}

// The grid potential e_m calculated for l=2 (Fourier)
// There is no specific equation for this in the article but it is:
// e_m = \sum_n KL_{m-n} q_n 
// The range for m and n obeys Eq. 23
void do_grid_to_grid_mapping_fourier() {
    for (int mx = - Mx/2 + 1 ; mx <= Mx/2 ; mx++) {
        for (int my = - My/2 + 1 ; my <= My/2 ; my++) {
            for (int mz = - Mz/2 + 1 ; mz <= Mz/2 ; mz++) {
                
                double sum = 0.0 ;                
                for (int nx = - Mx/2 + 1 ; nx <= Mx/2 ; nx++) {
                    for (int ny = - My/2 + 1 ; ny <= My/2 ; ny++) {
                        for (int nz = - Mz/2 + 1 ; nz <= Mz/2 ; nz++) {
                            
                        	// KL(m-n) needed
                            int mnx = mx - nx ;
                            int mny = my - ny ;
                            int mnz = mz - nz ;
                            
                            // KL is periodic
                            if (mnx < - Mx/2 + 1) { mnx += Mx ; } while (mnx < - Mx/2 + 1) ;
                            if (mnx >   Mx/2    ) { mnx -= Mx ; } while (mnx >   Mx/2    ) ;
                            
                            if (mny < - My/2 + 1) { mny += My ; } while (mny < - My/2 + 1) ;
                            if (mny >   My/2    ) { mny -= My ; } while (mny >   My/2    ) ;
                            
                            if (mnz < - Mz/2 + 1) { mnz += Mz ; } while (mnz < - Mz/2 + 1) ;
                            if (mnz >   Mz/2    ) { mnz -= Mz ; } while (mnz >   Mz/2    ) ;
                            
                            sum += KL[mnx + Mx/2 - 1][mny + My/2 - 1][mnz + Mz/2 - 1] *
                                   qm[nx  + Mx/2 - 1][ny  + My/2 - 1][nz  + Mz/2 - 1] ;
                        }
                    }
                }
                em_fourier[mx + Mx/2 - 1][my + My/2 - 1][mz + Mz/2 - 1] = sum ;
            }
        }
    }
}

// Long-range real-sum is just 0.5 * qm^T * em
void calculate_ulong_real() {
	// I need stencil
    calculate_stencil_Kl();
    // then grid potentials 
    do_grid_to_grid_mapping();
    
    ulong_real = 0.0;
    for (int mx = 0 ; mx < Mx ; mx++) {
        for (int my = 0 ; my < My ; my++) {
            for (int mz = 0 ; mz < Mz ; mz++) {
                ulong_real += 0.5 * qm[mx][my][mz] * em[mx][my][mz] ;
            }
        }
    }
}

// Centered B-spline with order v - 1
// Section 2.2.1
double Phi(double t) {
    return Bspline(v,t + v/2);
}

/* Eq. 35 */
double calculate_c(double k,double M) {
    double c = Phi(0);
    // Use the fact that sin components cancel
    for (int m = 1 ; m <= v/2 - 1 ; m++) {
        c += 2 * cos(2*MYPI * k * m / M) * Phi(m) ;
    }
    
    return 1/c ;
}

// Section 2.3, Eq. 29
// The range for m obeys Eq. 23
void calculate_stencil_KL() {
    for (int mx = - Mx/2 + 1 ; mx <= Mx/2 ; mx++) {
        for (int my = - My/2 + 1 ; my <= My/2 ; my++) {
            for (int mz = - Mz/2 + 1 ; mz <= Mz/2 ; mz++) {
                
                double sum = 0.0 ;
                for (int kx = -kmax ; kx <= kmax ; kx++) {
                    for (int ky = -kmax ; ky <= kmax ; ky++) {
                        for (int kz = -kmax ; kz <= kmax ; kz++) {
                            if (kx == 0 && ky == 0 && kz == 0)
                                continue;
                            double kvecx = kx / Ax ; // kvec = inv(A) * k
                            double kvecy = ky / Ay ;
                            double kvecz = kz / Az ;
                            double k2 = kvecx * kvecx + kvecy * kvecy + kvecz * kvecz ;
                            // Eq. 5
                            double chi = (1.0/(MYPI * k2 * detA)) * exp(-MYPI * MYPI * k2 / (beta * beta));
                            double dotprod = kx * mx / (double) Mx +
                                             ky * my / (double) My +
                                             kz * mz / (double) Mz ;
                            
                            double cx = calculate_c(kx, Mx);
                            double cy = calculate_c(ky, My);
                            double cz = calculate_c(kz, Mz);
                            
                            double c2 = cx * cx * cy * cy * cz * cz ;
                            sum += chi * c2 * cos(2*MYPI*dotprod) ;
                        }
                    }
                }
                KL[mx + Mx/2 - 1][my + My/2 - 1][mz + Mz/2 - 1]  = sum ;
            }
        }
    }
}

// Potential energy for Long-range recp-sum is
// just 0.5 * em_fourier^T * qm where
void calculate_ulong_four() {
	
	// Need stencil l=2 
    calculate_stencil_KL();
    // Need grid potentials when l=2
    do_grid_to_grid_mapping_fourier();
    
    ulong_four = 0.0;
    for (int mx = 0 ; mx < Mx ; mx++) {
        for (int my = 0 ; my < My ; my++) {
            for (int mz = 0 ; mz < Mz ; mz++) {
                ulong_four += 0.5 * qm[mx][my][mz] * em_fourier[mx][my][mz] ;
            }
        }
    }
}

// Self-explanatory
void display_results(){
    printf("%-28s : %-20s\n","Testing...",dataname);
    printf("%-28s : %-d\n","N" ,N);
    printf("%-28s : %-d\n","mu",mu);
    printf("%-28s : %-d\n","v" ,v);
    printf("%-28s : %-d\n","kmax" ,kmax);
    printf("%-28s : %-d\n","pmax" ,pmax);
    printf("%-28s : %-5.2f %-5.2f %-5.2f\n","hx hy hz" ,hx,hy,hz);
    printf("%-28s : %-5.2f %-5.2f %-5.2f\n","Ax Ay Az" ,Ax,Ay,Az);
    printf("%-28s : %-5d %-5d %-5d\n","Mx My Mz" ,Mx,My,Mz);
    printf("%-28s : %-25.16e\n","beta" ,beta);
    printf("%-28s : %-25.16e\n","detA" ,detA);
    printf("%-28s : %-25.16e\n","abar" ,abar);
    printf("%-28s : %-25.16e\n","a" ,a);
    printf("%-28s : %-25.16e\n","aL" ,aL);
    printf("%-28s : %25.16e\n","ushort_real" ,ushort_real);
    printf("%-28s : %25.16e\n","ushort_self" ,ushort_self);
    printf("%-28s : %25.16e\n","ushort_csr"  ,ushort_csr);
    printf("%-28s : %25.16e\n","ulong_self"  ,ulong_self);
    printf("%-28s : %25.16e\n","ulong_real"  ,ulong_real);
    printf("%-28s : %25.16e\n","ulong_real_expected"  ,ulong_real_expected);
    printf("%-28s : %25.16e\n","ulong_four"  ,ulong_four);
    printf("%-28s : %25.16e\n","ulong_four_expected_cossum"  ,ulong_four_expected_cossum);
    printf("%-28s : %25.16e\n","ulong_four_expected_sinsum"  ,ulong_four_expected_sinsum);
    printf("%-28s : %25.16e\n","utotal",utotal);
    printf("%-28s : %25.16e\n","utotal_expected",utotal_expected);   
}

// I numbered all 12 benchmarks with an integer
// I load the related data and apply MSM
void run_benchmark(int id) {
    switch (id) {
        case 1:
            load_benchmark_NaNaN8();
            break;
        case 2:
            load_benchmark_NaNaN64();
            break;
        case 3:
            load_benchmark_NaNaN512();
            break;
        case 4:
            load_benchmark_NaClN8();
            break;
        case 5:
            load_benchmark_NaClN64();
            break;
        case 6:
            load_benchmark_NaClN512();
            break;
        case 7:
            load_benchmark_CsClN2();
            break;
        case 8:
            load_benchmark_CsClN16();
            break;
        case 9:
            load_benchmark_CsClN128();
            break;
        case 10:
            load_benchmark_changaN8();
            break;
        case 11:
            load_benchmark_changaN64();
            break;
        case 12:
            load_benchmark_changaN512();
            break;
        default:
            break;
    }

    do_anterpolation();

    // There is no approximated versions of these four quantities
    calculate_ushort_real();
    calculate_ushort_self();
    calculate_ushort_csr();    
    calculate_ulong_self();
    
    // ulong_real and ulong_four are calculated by using interpolation
    calculate_ulong_real();
    calculate_ulong_four();
    
    // For comparioson there two are direct calculations
    calculate_ulong_real_expected();
    calculate_ulong_four_expected();
    
    // Total potential energy when using interpolation for ulong_real and ulong_fourier
    utotal  = ushort_real + ushort_self + ushort_csr +  ulong_self + 
    		ulong_real + ulong_four ;
    
    // Total potential energy without interpolation
    utotal_expected = ushort_real + ushort_self + ushort_csr + ulong_self +
    				  ulong_real_expected + ulong_four_expected  ;

    display_results();
}

// A very straightforward recursive implementation of B-splines in 2.2.1
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

// N=8 case of Figure 1 in report20180327
void load_benchmark_NaNaN8() {
    sprintf(dataname, "NaNaN8");
    N    = 8;
    Ax   = Ay = Az = 2;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10; 
    pmax = 10;
    h    = hL = hx = hy = hz = 1;
    Mx   = My = Mz = 2;
    a = abar * h ;
    aL = 2 * a ;
    beta = 3.9952779973449493e-01; // precalculated w.r.t Appendix B.1
    
    int nx = 2; int ny = 2 ; int nz = 2 ;
    int bodyindex = 0 ;
    for (int i = 0 ; i < nx ; i++) {
        for (int j = 0 ; j < ny ; j++) {
            for (int k = 0 ; k < nz ; k++) {
                r[bodyindex][X] = i ;
                r[bodyindex][Y] = j ;
                r[bodyindex][Z] = k ;
                q[bodyindex]    = 1 ;
                bodyindex++;
            }
        }
    }
}

// N=64 case of Figure 1 in report20180327
void load_benchmark_NaNaN64() {
    sprintf(dataname, "NaNaN64");
    N    = 64;
    Ax   = Ay = Az = 4;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 1;
    Mx   = My = Mz = 4;
    a = abar * h ;
    aL = 2 * a ;
    beta = 3.9952779973449493e-01;
    
    int nx = 4; int ny = 4 ; int nz = 4 ;
    int bodyindex = 0 ;
    for (int i = 0 ; i < nx ; i++) {
        for (int j = 0 ; j < ny ; j++) {
            for (int k = 0 ; k < nz ; k++) {
                r[bodyindex][X] = i ;
                r[bodyindex][Y] = j ;
                r[bodyindex][Z] = k ;
                q[bodyindex]    = 1 ;
                bodyindex++;
            }
        }
    }
}

// N=512 case of Figure 1 in report20180327
void load_benchmark_NaNaN512() {
    sprintf(dataname, "NaNaN512");
    N    = 512;
    Ax   = Ay = Az = 8;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 1;
    Mx   = My = Mz = 8;
    a = abar * h ;
    aL = 2 * a ;
    beta = 3.9952779973449493e-01;
    
    int nx = 8; int ny = 8 ; int nz = 8 ;
    int bodyindex = 0 ;
    for (int i = 0 ; i < nx ; i++) {
        for (int j = 0 ; j < ny ; j++) {
            for (int k = 0 ; k < nz ; k++) {
                r[bodyindex][X] = i ;
                r[bodyindex][Y] = j ;
                r[bodyindex][Z] = k ;
                q[bodyindex]    = 1 ;
                bodyindex++;
            }
        }
    }
}

// N=8 case of Figure 2 in report20180327
void load_benchmark_NaClN8() {
    sprintf(dataname, "NaClN8");
    N    = 8;
    Ax   = Ay = Az = 2;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 1;
    Mx   = My = Mz = 2;
    a = abar * h ;
    aL = 2 * a ;
    beta = 3.9952779973449493e-01;
    
    int nx = 2; int ny = 2 ; int nz = 2 ;
    int bodyindex = 0 ;
    for (int i = 0 ; i < nx ; i++) {
        for (int j = 0 ; j < ny ; j++) {
            for (int k = 0 ; k < nz ; k++) {
                r[bodyindex][X] = i ;
                r[bodyindex][Y] = j ;
                r[bodyindex][Z] = k ;
                q[bodyindex]    = pow(-1,i+j+k) ;
                bodyindex++;
            }
        }
    }
}

// N=64 case of Figure 2 in report20180327
void load_benchmark_NaClN64() {
    sprintf(dataname, "NaClN64");
    N    = 64;
    Ax   = Ay = Az = 4;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 1;
    Mx   = My = Mz = 4;
    a = abar * h ;
    aL = 2 * a ;
    beta = 3.9952779973449493e-01;
    
    int nx = 4; int ny = 4 ; int nz = 4 ;
    int bodyindex = 0 ;
    for (int i = 0 ; i < nx ; i++) {
        for (int j = 0 ; j < ny ; j++) {
            for (int k = 0 ; k < nz ; k++) {
                r[bodyindex][X] = i ;
                r[bodyindex][Y] = j ;
                r[bodyindex][Z] = k ;
                q[bodyindex]    = pow(-1,i+j+k) ;
                bodyindex++;
            }
        }
    }
}

// N=512 case of Figure 2 in report20180327
void load_benchmark_NaClN512() {
    sprintf(dataname, "NaClN512");
    N    = 512;
    Ax   = Ay = Az = 8;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 1;
    Mx   = My = Mz = 8;
    a = abar * h ;
    aL = 2 * a ;
    beta = 3.9952779973449493e-01;
    
    int nx = 8; int ny = 8 ; int nz = 8 ;
    int bodyindex = 0 ;
    for (int i = 0 ; i < nx ; i++) {
        for (int j = 0 ; j < ny ; j++) {
            for (int k = 0 ; k < nz ; k++) {
                r[bodyindex][X] = i ;
                r[bodyindex][Y] = j ;
                r[bodyindex][Z] = k ;
                q[bodyindex]    = pow(-1,i+j+k) ;
                bodyindex++;
            }
        }
    }
}

// N=2 case of Figure 3 in report20180327
void load_benchmark_CsClN2() {
    sprintf(dataname, "CsClN2");
    N    = 2;
    Ax   = Ay = Az = 1;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 0.5;
    Mx   = My = Mz = 2;
    a = abar * h ;
    aL = 2 * a ;
    beta = 5.1015658840559586e-01;
    
    int nx = 1; int ny = 1 ; int nz = 1 ;
    int bodyindex = 0 ;
    for (int i = 0 ; i < nx ; i++) {
        for (int j = 0 ; j < ny ; j++) {
            for (int k = 0 ; k < nz ; k++) {
                // Cs Atom
                r[bodyindex][X] = i ;
                r[bodyindex][Y] = j ;
                r[bodyindex][Z] = k ;
                q[bodyindex]    = 1 ;
                bodyindex++;
                // Cl Atom
                r[bodyindex][X] = i + 0.5 ;
                r[bodyindex][Y] = j + 0.5 ;
                r[bodyindex][Z] = k + 0.5 ;
                q[bodyindex]    = -1 ;
                bodyindex++;
            }
        }
    }
}

// N=16 case of Figure 3 in report20180327
void load_benchmark_CsClN16() {
    sprintf(dataname, "CsClN16");

    N    = 16;
    Ax   = Ay = Az = 2;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 0.5;
    Mx   = My = Mz = 4;
    a = abar * h ;
    aL = 2 * a ;
    beta = 5.1015658840559586e-01;
    
    int nx = 2; int ny = 2 ; int nz = 2 ;
    int bodyindex = 0 ;
    for (int i = 0 ; i < nx ; i++) {
        for (int j = 0 ; j < ny ; j++) {
            for (int k = 0 ; k < nz ; k++) {
                // Cs Atom
                r[bodyindex][X] = i ;
                r[bodyindex][Y] = j ;
                r[bodyindex][Z] = k ;
                q[bodyindex]    = 1 ;
                bodyindex++;
                // Cl Atom
                r[bodyindex][X] = i + 0.5 ;
                r[bodyindex][Y] = j + 0.5 ;
                r[bodyindex][Z] = k + 0.5 ;
                q[bodyindex]    = -1 ;
                bodyindex++;
            }
        }
    }
}

// N=128 case of Figure 3 in report20180327
void load_benchmark_CsClN128() {
    sprintf(dataname, "CsClN128");
    N    = 128;
    Ax   = Ay = Az = 4;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 0.5;
    Mx   = My = Mz = 8;
    a = abar * h ;
    aL = 2 * a ;
    beta = 5.1015658840559586e-01;
    
    int nx = 4; int ny = 4 ; int nz = 4 ;
    int bodyindex = 0 ;
    for (int i = 0 ; i < nx ; i++) {
        for (int j = 0 ; j < ny ; j++) {
            for (int k = 0 ; k < nz ; k++) {
                // Cs Atom
                r[bodyindex][X] = i ;
                r[bodyindex][Y] = j ;
                r[bodyindex][Z] = k ;
                q[bodyindex]    = 1 ;
                bodyindex++;
                // Cl Atom
                r[bodyindex][X] = i + 0.5 ;
                r[bodyindex][Y] = j + 0.5 ;
                r[bodyindex][Z] = k + 0.5 ;
                q[bodyindex]    = -1 ;
                bodyindex++;
            }
        }
    }
}

// N=8 randomly distributed particles in unit-cube
void load_benchmark_changaN8() {
    sprintf(dataname, "changaN8");
    N    = 8;
    Ax   = Ay = Az = 1;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 0.5 ;
    Mx   = My = Mz = 2;
    a = abar * h ;
    aL = 2 * a ;
    beta = 8.2233267459052750e-01 ;
    r[ 0][X]= 7.892292e-01; r[ 0][Y]= 5.422255e-01; r[ 0][Z]= 1.832440e-01; q[  0]= 9.993093e-06;
    r[ 1][X]= 3.162612e-01; r[ 1][Y]= 9.933904e-01; r[ 1][Z]= 5.268250e-02; q[  1]= 9.993093e-06;
    r[ 2][X]= 6.496070e-02; r[ 2][Y]= 2.662055e-01; r[ 2][Z]= 4.990308e-01; q[  2]= 9.993093e-06;
    r[ 3][X]= 2.377985e-01; r[ 3][Y]= 7.406127e-01; r[ 3][Z]= 9.129351e-01; q[  3]= 9.993093e-06;
    r[ 4][X]= 9.024968e-01; r[ 4][Y]= 1.365204e-01; r[ 4][Z]= 8.475435e-01; q[  4]= 9.993093e-06;
    r[ 5][X]= 6.892704e-01; r[ 5][Y]= 2.344152e-01; r[ 5][Z]= 5.649641e-01; q[  5]= 9.993093e-06;
    r[ 6][X]= 2.313120e-01; r[ 6][Y]= 7.977766e-01; r[ 6][Z]= 4.078127e-01; q[  6]= 9.993093e-06;
    r[ 7][X]= 4.682633e-01; r[ 7][Y]= 8.296577e-01; r[ 7][Z]= 7.368002e-01; q[  7]= 9.993093e-06;
}

// N=64 randomly distributed particles in unit-cube
void load_benchmark_changaN64() {
    sprintf(dataname, "changaN64");

    N    = 64;
    Ax   = Ay = Az = 1;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 0.25 ;
    Mx   = My = Mz = 4;
    a = abar * h ;
    aL = 2 * a ;
    beta = 1.6899842461315719e+00 ;
    r[ 0][X]= 6.851051e-01; r[ 0][Y]= 3.284786e-01; r[ 0][Z]= 8.348171e-01; q[  0]= 9.993093e-06;
    r[ 1][X]= 1.227815e-01; r[ 1][Y]= 3.274495e-01; r[ 1][Z]= 9.072554e-01; q[  1]= 9.993093e-06;
    r[ 2][X]= 1.565225e-01; r[ 2][Y]= 9.730324e-01; r[ 2][Z]= 3.278511e-01; q[  2]= 9.993093e-06;
    r[ 3][X]= 8.428839e-01; r[ 3][Y]= 9.244588e-01; r[ 3][Z]= 9.955087e-01; q[  3]= 9.993093e-06;
    r[ 4][X]= 3.789962e-01; r[ 4][Y]= 6.495880e-02; r[ 4][Z]= 7.120740e-02; q[  4]= 9.993093e-06;
    r[ 5][X]= 3.670836e-01; r[ 5][Y]= 1.373245e-01; r[ 5][Z]= 9.765455e-01; q[  5]= 9.993093e-06;
    r[ 6][X]= 5.134624e-01; r[ 6][Y]= 5.680366e-01; r[ 6][Z]= 4.476654e-01; q[  6]= 9.993093e-06;
    r[ 7][X]= 4.656702e-01; r[ 7][Y]= 9.208753e-01; r[ 7][Z]= 5.902060e-01; q[  7]= 9.993093e-06;
    r[ 8][X]= 4.627095e-01; r[ 8][Y]= 2.438240e-01; r[ 8][Z]= 8.246222e-01; q[  8]= 9.993093e-06;
    r[ 9][X]= 4.467418e-01; r[ 9][Y]= 8.718548e-01; r[ 9][Z]= 3.926372e-01; q[  9]= 9.993093e-06;
    r[10][X]= 7.892292e-01; r[10][Y]= 5.422255e-01; r[10][Z]= 1.832440e-01; q[ 10]= 9.993093e-06;
    r[11][X]= 4.522532e-01; r[11][Y]= 6.574540e-02; r[11][Z]= 6.650430e-01; q[ 11]= 9.993093e-06;
    r[12][X]= 4.707406e-01; r[12][Y]= 2.364980e-01; r[12][Z]= 8.156511e-01; q[ 12]= 9.993093e-06;
    r[13][X]= 8.840847e-01; r[13][Y]= 4.808207e-01; r[13][Z]= 8.363834e-01; q[ 13]= 9.993093e-06;
    r[14][X]= 3.758621e-01; r[14][Y]= 5.297577e-01; r[14][Z]= 2.964660e-02; q[ 14]= 9.993093e-06;
    r[15][X]= 4.907990e-01; r[15][Y]= 7.365370e-02; r[15][Z]= 4.241293e-01; q[ 15]= 9.993093e-06;
    r[16][X]= 1.887263e-01; r[16][Y]= 4.982793e-01; r[16][Z]= 6.025294e-01; q[ 16]= 9.993093e-06;
    r[17][X]= 5.530562e-01; r[17][Y]= 8.674165e-01; r[17][Z]= 4.477672e-01; q[ 17]= 9.993093e-06;
    r[18][X]= 4.670775e-01; r[18][Y]= 2.087175e-01; r[18][Z]= 8.113433e-01; q[ 18]= 9.993093e-06;
    r[19][X]= 4.987226e-01; r[19][Y]= 6.832300e-03; r[19][Z]= 3.617589e-01; q[ 19]= 9.993093e-06;
    r[20][X]= 3.533229e-01; r[20][Y]= 8.994851e-01; r[20][Z]= 2.276869e-01; q[ 20]= 9.993093e-06;
    r[21][X]= 7.121460e-01; r[21][Y]= 1.402387e-01; r[21][Z]= 4.517867e-01; q[ 21]= 9.993093e-06;
    r[22][X]= 7.694070e-02; r[22][Y]= 5.786474e-01; r[22][Z]= 9.419320e-01; q[ 22]= 9.993093e-06;
    r[23][X]= 2.960440e-01; r[23][Y]= 5.241426e-01; r[23][Z]= 8.570910e-01; q[ 23]= 9.993093e-06;
    r[24][X]= 9.672149e-01; r[24][Y]= 3.494468e-01; r[24][Z]= 6.679157e-01; q[ 24]= 9.993093e-06;
    r[25][X]= 1.285045e-01; r[25][Y]= 9.609555e-01; r[25][Z]= 6.623552e-01; q[ 25]= 9.993093e-06;
    r[26][X]= 6.459974e-01; r[26][Y]= 4.625895e-01; r[26][Z]= 5.478487e-01; q[ 26]= 9.993093e-06;
    r[27][X]= 7.286598e-01; r[27][Y]= 7.467237e-01; r[27][Z]= 6.556760e-02; q[ 27]= 9.993093e-06;
    r[28][X]= 7.728232e-01; r[28][Y]= 7.686914e-01; r[28][Z]= 9.938379e-01; q[ 28]= 9.993093e-06;
    r[29][X]= 3.586270e-02; r[29][Y]= 3.116630e-02; r[29][Z]= 7.758842e-01; q[ 29]= 9.993093e-06;
    r[30][X]= 3.902618e-01; r[30][Y]= 4.651819e-01; r[30][Z]= 9.060126e-01; q[ 30]= 9.993093e-06;
    r[31][X]= 1.205499e-01; r[31][Y]= 5.451384e-01; r[31][Z]= 9.054577e-01; q[ 31]= 9.993093e-06;
    r[32][X]= 4.931200e-01; r[32][Y]= 4.095510e-01; r[32][Z]= 9.231425e-01; q[ 32]= 9.993093e-06;
    r[33][X]= 2.463758e-01; r[33][Y]= 7.901024e-01; r[33][Z]= 3.987949e-01; q[ 33]= 9.993093e-06;
    r[34][X]= 5.832959e-01; r[34][Y]= 6.133283e-01; r[34][Z]= 3.026630e-01; q[ 34]= 9.993093e-06;
    r[35][X]= 3.162612e-01; r[35][Y]= 9.933904e-01; r[35][Z]= 5.268250e-02; q[ 35]= 9.993093e-06;
    r[36][X]= 2.351270e-01; r[36][Y]= 6.982184e-01; r[36][Z]= 8.515429e-01; q[ 36]= 9.993093e-06;
    r[37][X]= 6.496070e-02; r[37][Y]= 2.662055e-01; r[37][Z]= 4.990308e-01; q[ 37]= 9.993093e-06;
    r[38][X]= 1.187146e-01; r[38][Y]= 2.162690e-02; r[38][Z]= 4.467074e-01; q[ 38]= 9.993093e-06;
    r[39][X]= 8.478030e-01; r[39][Y]= 9.851480e-02; r[39][Z]= 9.050156e-01; q[ 39]= 9.993093e-06;
    r[40][X]= 1.794536e-01; r[40][Y]= 1.444035e-01; r[40][Z]= 9.003375e-01; q[ 40]= 9.993093e-06;
    r[41][X]= 7.776713e-01; r[41][Y]= 7.199529e-01; r[41][Z]= 4.833861e-01; q[ 41]= 9.993093e-06;
    r[42][X]= 8.906200e-03; r[42][Y]= 2.773881e-01; r[42][Z]= 8.253802e-01; q[ 42]= 9.993093e-06;
    r[43][X]= 3.691860e-02; r[43][Y]= 9.881868e-01; r[43][Z]= 8.241665e-01; q[ 43]= 9.993093e-06;
    r[44][X]= 6.536385e-01; r[44][Y]= 9.965687e-01; r[44][Z]= 8.599800e-03; q[ 44]= 9.993093e-06;
    r[45][X]= 3.916615e-01; r[45][Y]= 2.249520e-01; r[45][Z]= 6.914474e-01; q[ 45]= 9.993093e-06;
    r[46][X]= 7.047238e-01; r[46][Y]= 3.997731e-01; r[46][Z]= 6.034167e-01; q[ 46]= 9.993093e-06;
    r[47][X]= 7.897092e-01; r[47][Y]= 9.386740e-01; r[47][Z]= 2.204521e-01; q[ 47]= 9.993093e-06;
    r[48][X]= 4.507226e-01; r[48][Y]= 2.807021e-01; r[48][Z]= 8.755396e-01; q[ 48]= 9.993093e-06;
    r[49][X]= 5.852413e-01; r[49][Y]= 1.264235e-01; r[49][Z]= 2.852339e-01; q[ 49]= 9.993093e-06;
    r[50][X]= 2.377985e-01; r[50][Y]= 7.406127e-01; r[50][Z]= 9.129351e-01; q[ 50]= 9.993093e-06;
    r[51][X]= 9.024968e-01; r[51][Y]= 1.365204e-01; r[51][Z]= 8.475435e-01; q[ 51]= 9.993093e-06;
    r[52][X]= 3.954231e-01; r[52][Y]= 5.172741e-01; r[52][Z]= 2.190200e-03; q[ 52]= 9.993093e-06;
    r[53][X]= 3.673120e-01; r[53][Y]= 9.115603e-01; r[53][Z]= 4.667580e-02; q[ 53]= 9.993093e-06;
    r[54][X]= 6.892704e-01; r[54][Y]= 2.344152e-01; r[54][Z]= 5.649641e-01; q[ 54]= 9.993093e-06;
    r[55][X]= 3.916906e-01; r[55][Y]= 3.384480e-02; r[55][Z]= 2.363765e-01; q[ 55]= 9.993093e-06;
    r[56][X]= 2.313120e-01; r[56][Y]= 7.977766e-01; r[56][Z]= 4.078127e-01; q[ 56]= 9.993093e-06;
    r[57][X]= 9.807751e-01; r[57][Y]= 5.710414e-01; r[57][Z]= 9.844820e-01; q[ 57]= 9.993093e-06;
    r[58][X]= 2.856527e-01; r[58][Y]= 6.308636e-01; r[58][Z]= 1.678970e-02; q[ 58]= 9.993093e-06;
    r[59][X]= 4.030789e-01; r[59][Y]= 5.240499e-01; r[59][Z]= 9.669881e-01; q[ 59]= 9.993093e-06;
    r[60][X]= 4.682633e-01; r[60][Y]= 8.296577e-01; r[60][Z]= 7.368002e-01; q[ 60]= 9.993093e-06;
    r[61][X]= 4.733331e-01; r[61][Y]= 3.019186e-01; r[61][Z]= 8.890947e-01; q[ 61]= 9.993093e-06;
    r[62][X]= 7.012500e-02; r[62][Y]= 4.189419e-01; r[62][Z]= 8.550628e-01; q[ 62]= 9.993093e-06;
    r[63][X]= 6.721672e-01; r[63][Y]= 1.791827e-01; r[63][Z]= 2.731031e-01; q[ 63]= 9.993093e-06;
}

// N=512 randomly distributed particles in unit-cube
void load_benchmark_changaN512() {
    sprintf(dataname, "changaN512");
    N    = 512;
    Ax   = Ay = Az = 1;
    detA = Ax * Ay * Az ;
    v    = 4; // Please don't change v=4 because wprime
              // in the code is only valid for v=4
    abar = 6;
    mu   = 10;
    kmax = 10;
    pmax = 10;
    h    = hL = hx = hy = hz = 0.125 ;
    Mx   = My = Mz = 8;
    a = abar * h ;
    aL = 2 * a ;
    beta = 3.4683206143698224e+00;
    r[ 0][X]= 6.911842e-01; r[ 0][Y]= 7.251611e-01; r[ 0][Z]= 2.556736e-01; q[  0]= 9.993093e-06;
    r[ 1][X]= 4.153972e-01; r[ 1][Y]= 1.923024e-01; r[ 1][Z]= 6.005900e-02; q[  1]= 9.993093e-06;
    r[ 2][X]= 7.012648e-01; r[ 2][Y]= 7.103131e-01; r[ 2][Z]= 3.901396e-01; q[  2]= 9.993093e-06;
    r[ 3][X]= 6.227630e-01; r[ 3][Y]= 8.576262e-01; r[ 3][Z]= 3.998055e-01; q[  3]= 9.993093e-06;
    r[ 4][X]= 8.350638e-01; r[ 4][Y]= 2.989518e-01; r[ 4][Z]= 7.021725e-01; q[  4]= 9.993093e-06;
    r[ 5][X]= 4.619036e-01; r[ 5][Y]= 3.093388e-01; r[ 5][Z]= 8.828159e-01; q[  5]= 9.993093e-06;
    r[ 6][X]= 6.851051e-01; r[ 6][Y]= 3.284786e-01; r[ 6][Z]= 8.348171e-01; q[  6]= 9.993093e-06;
    r[ 7][X]= 7.917247e-01; r[ 7][Y]= 3.678188e-01; r[ 7][Z]= 6.427480e-01; q[  7]= 9.993093e-06;
    r[ 8][X]= 7.004844e-01; r[ 8][Y]= 1.058245e-01; r[ 8][Z]= 8.226223e-01; q[  8]= 9.993093e-06;
    r[ 9][X]= 2.207296e-01; r[ 9][Y]= 6.409511e-01; r[ 9][Z]= 9.236183e-01; q[  9]= 9.993093e-06;
    r[10][X]= 1.227815e-01; r[10][Y]= 3.274495e-01; r[10][Z]= 9.072554e-01; q[ 10]= 9.993093e-06;
    r[11][X]= 9.699941e-01; r[11][Y]= 3.776225e-01; r[11][Z]= 1.059234e-01; q[ 11]= 9.993093e-06;
    r[12][X]= 7.812613e-01; r[12][Y]= 6.784635e-01; r[12][Z]= 1.438599e-01; q[ 12]= 9.993093e-06;
    r[13][X]= 8.264570e-01; r[13][Y]= 4.268751e-01; r[13][Z]= 2.837529e-01; q[ 13]= 9.993093e-06;
    r[14][X]= 1.933017e-01; r[14][Y]= 8.586112e-01; r[14][Z]= 7.379635e-01; q[ 14]= 9.993093e-06;
    r[15][X]= 2.407870e-02; r[15][Y]= 2.565242e-01; r[15][Z]= 8.870084e-01; q[ 15]= 9.993093e-06;
    r[16][X]= 9.819346e-01; r[16][Y]= 3.903937e-01; r[16][Z]= 4.563230e-01; q[ 16]= 9.993093e-06;
    r[17][X]= 1.365981e-01; r[17][Y]= 7.537125e-01; r[17][Z]= 6.112866e-01; q[ 17]= 9.993093e-06;
    r[18][X]= 7.600073e-01; r[18][Y]= 4.773000e-02; r[18][Z]= 9.509733e-01; q[ 18]= 9.993093e-06;
    r[19][X]= 2.838080e-01; r[19][Y]= 6.396568e-01; r[19][Z]= 2.564560e-02; q[ 19]= 9.993093e-06;
    r[20][X]= 2.333143e-01; r[20][Y]= 7.827117e-01; r[20][Z]= 8.730019e-01; q[ 20]= 9.993093e-06;
    r[21][X]= 4.594206e-01; r[21][Y]= 9.183779e-01; r[21][Z]= 5.843892e-01; q[ 21]= 9.993093e-06;
    r[22][X]= 1.565225e-01; r[22][Y]= 9.730324e-01; r[22][Z]= 3.278511e-01; q[ 22]= 9.993093e-06;
    r[23][X]= 7.195269e-01; r[23][Y]= 2.134463e-01; r[23][Z]= 6.175054e-01; q[ 23]= 9.993093e-06;
    r[24][X]= 8.168274e-01; r[24][Y]= 6.851196e-01; r[24][Z]= 7.640840e-02; q[ 24]= 9.993093e-06;
    r[25][X]= 9.247192e-01; r[25][Y]= 5.539906e-01; r[25][Z]= 4.200297e-01; q[ 25]= 9.993093e-06;
    r[26][X]= 8.428839e-01; r[26][Y]= 9.244588e-01; r[26][Z]= 9.955087e-01; q[ 26]= 9.993093e-06;
    r[27][X]= 8.047912e-01; r[27][Y]= 8.214058e-01; r[27][Z]= 7.098020e-01; q[ 27]= 9.993093e-06;
    r[28][X]= 6.420638e-01; r[28][Y]= 7.244584e-01; r[28][Z]= 2.613080e-01; q[ 28]= 9.993093e-06;
    r[29][X]= 5.250883e-01; r[29][Y]= 4.104060e-02; r[29][Z]= 9.524418e-01; q[ 29]= 9.993093e-06;
    r[30][X]= 6.475008e-01; r[30][Y]= 7.320547e-01; r[30][Z]= 2.544812e-01; q[ 30]= 9.993093e-06;
    r[31][X]= 7.879515e-01; r[31][Y]= 3.560017e-01; r[31][Z]= 6.422039e-01; q[ 31]= 9.993093e-06;
    r[32][X]= 5.652663e-01; r[32][Y]= 6.056150e-01; r[32][Z]= 7.031073e-01; q[ 32]= 9.993093e-06;
    r[33][X]= 5.761917e-01; r[33][Y]= 7.729439e-01; r[33][Z]= 5.224199e-01; q[ 33]= 9.993093e-06;
    r[34][X]= 4.911675e-01; r[34][Y]= 1.808594e-01; r[34][Z]= 8.748586e-01; q[ 34]= 9.993093e-06;
    r[35][X]= 3.789962e-01; r[35][Y]= 6.495880e-02; r[35][Z]= 7.120740e-02; q[ 35]= 9.993093e-06;
    r[36][X]= 1.050423e-01; r[36][Y]= 1.896743e-01; r[36][Z]= 4.296410e-02; q[ 36]= 9.993093e-06;
    r[37][X]= 3.157304e-01; r[37][Y]= 5.648595e-01; r[37][Z]= 3.219040e-02; q[ 37]= 9.993093e-06;
    r[38][X]= 3.697639e-01; r[38][Y]= 3.220214e-01; r[38][Z]= 7.809486e-01; q[ 38]= 9.993093e-06;
    r[39][X]= 6.725738e-01; r[39][Y]= 3.606801e-01; r[39][Z]= 5.641725e-01; q[ 39]= 9.993093e-06;
    r[40][X]= 3.670836e-01; r[40][Y]= 1.373245e-01; r[40][Z]= 9.765455e-01; q[ 40]= 9.993093e-06;
    r[41][X]= 2.461173e-01; r[41][Y]= 4.916917e-01; r[41][Z]= 7.098064e-01; q[ 41]= 9.993093e-06;
    r[42][X]= 5.561610e-02; r[42][Y]= 2.658020e-01; r[42][Z]= 6.483387e-01; q[ 42]= 9.993093e-06;
    r[43][X]= 8.762310e-01; r[43][Y]= 6.076758e-01; r[43][Z]= 8.423160e-01; q[ 43]= 9.993093e-06;
    r[44][X]= 2.874427e-01; r[44][Y]= 9.726147e-01; r[44][Z]= 3.774505e-01; q[ 44]= 9.993093e-06;
    r[45][X]= 5.134624e-01; r[45][Y]= 5.680366e-01; r[45][Z]= 4.476654e-01; q[ 45]= 9.993093e-06;
    r[46][X]= 3.948462e-01; r[46][Y]= 2.404915e-01; r[46][Z]= 6.830163e-01; q[ 46]= 9.993093e-06;
    r[47][X]= 9.944960e-02; r[47][Y]= 4.820007e-01; r[47][Z]= 9.236963e-01; q[ 47]= 9.993093e-06;
    r[48][X]= 8.579296e-01; r[48][Y]= 5.242882e-01; r[48][Z]= 8.531370e-02; q[ 48]= 9.993093e-06;
    r[49][X]= 7.796355e-01; r[49][Y]= 6.514324e-01; r[49][Z]= 2.923591e-01; q[ 49]= 9.993093e-06;
    r[50][X]= 6.665321e-01; r[50][Y]= 3.061855e-01; r[50][Z]= 7.032596e-01; q[ 50]= 9.993093e-06;
    r[51][X]= 4.048947e-01; r[51][Y]= 5.065629e-01; r[51][Z]= 1.434100e-03; q[ 51]= 9.993093e-06;
    r[52][X]= 3.790896e-01; r[52][Y]= 3.455354e-01; r[52][Z]= 6.540263e-01; q[ 52]= 9.993093e-06;
    r[53][X]= 1.882640e-02; r[53][Y]= 5.139580e-01; r[53][Z]= 5.481883e-01; q[ 53]= 9.993093e-06;
    r[54][X]= 9.265717e-01; r[54][Y]= 9.553683e-01; r[54][Z]= 8.062507e-01; q[ 54]= 9.993093e-06;
    r[55][X]= 9.216734e-01; r[55][Y]= 6.109082e-01; r[55][Z]= 1.586790e-02; q[ 55]= 9.993093e-06;
    r[56][X]= 2.682213e-01; r[56][Y]= 2.084542e-01; r[56][Z]= 9.856417e-01; q[ 56]= 9.993093e-06;
    r[57][X]= 8.895800e-03; r[57][Y]= 8.507967e-01; r[57][Z]= 5.898781e-01; q[ 57]= 9.993093e-06;
    r[58][X]= 4.376255e-01; r[58][Y]= 7.002272e-01; r[58][Z]= 5.845311e-01; q[ 58]= 9.993093e-06;
    r[59][X]= 2.804805e-01; r[59][Y]= 8.104540e-01; r[59][Z]= 2.804095e-01; q[ 59]= 9.993093e-06;
    r[60][X]= 4.656702e-01; r[60][Y]= 9.208753e-01; r[60][Z]= 5.902060e-01; q[ 60]= 9.993093e-06;
    r[61][X]= 4.627095e-01; r[61][Y]= 2.438240e-01; r[61][Z]= 8.246222e-01; q[ 61]= 9.993093e-06;
    r[62][X]= 4.467418e-01; r[62][Y]= 8.718548e-01; r[62][Z]= 3.926372e-01; q[ 62]= 9.993093e-06;
    r[63][X]= 4.333463e-01; r[63][Y]= 9.881022e-01; r[63][Z]= 6.078659e-01; q[ 63]= 9.993093e-06;
    r[64][X]= 2.365042e-01; r[64][Y]= 4.823340e-01; r[64][Z]= 7.061769e-01; q[ 64]= 9.993093e-06;
    r[65][X]= 7.723974e-01; r[65][Y]= 2.823441e-01; r[65][Z]= 7.073730e-01; q[ 65]= 9.993093e-06;
    r[66][X]= 2.605053e-01; r[66][Y]= 2.940750e-01; r[66][Z]= 2.300842e-01; q[ 66]= 9.993093e-06;
    r[67][X]= 7.892292e-01; r[67][Y]= 5.422255e-01; r[67][Z]= 1.832440e-01; q[ 67]= 9.993093e-06;
    r[68][X]= 2.446588e-01; r[68][Y]= 7.889431e-01; r[68][Z]= 4.026880e-01; q[ 68]= 9.993093e-06;
    r[69][X]= 2.972952e-01; r[69][Y]= 8.223604e-01; r[69][Z]= 4.702353e-01; q[ 69]= 9.993093e-06;
    r[70][X]= 4.522532e-01; r[70][Y]= 6.574540e-02; r[70][Z]= 6.650430e-01; q[ 70]= 9.993093e-06;
    r[71][X]= 2.760394e-01; r[71][Y]= 2.281119e-01; r[71][Z]= 5.609820e-02; q[ 71]= 9.993093e-06;
    r[72][X]= 2.678529e-01; r[72][Y]= 6.595889e-01; r[72][Z]= 2.658800e-02; q[ 72]= 9.993093e-06;
    r[73][X]= 2.446393e-01; r[73][Y]= 6.910035e-01; r[73][Z]= 8.751380e-01; q[ 73]= 9.993093e-06;
    r[74][X]= 1.055041e-01; r[74][Y]= 3.931574e-01; r[74][Z]= 2.519061e-01; q[ 74]= 9.993093e-06;
    r[75][X]= 4.707406e-01; r[75][Y]= 2.364980e-01; r[75][Z]= 8.156511e-01; q[ 75]= 9.993093e-06;
    r[76][X]= 8.840847e-01; r[76][Y]= 4.808207e-01; r[76][Z]= 8.363834e-01; q[ 76]= 9.993093e-06;
    r[77][X]= 9.225118e-01; r[77][Y]= 2.535498e-01; r[77][Z]= 5.312372e-01; q[ 77]= 9.993093e-06;
    r[78][X]= 1.793456e-01; r[78][Y]= 1.512681e-01; r[78][Z]= 5.523479e-01; q[ 78]= 9.993093e-06;
    r[79][X]= 6.160921e-01; r[79][Y]= 5.336020e-01; r[79][Z]= 5.599179e-01; q[ 79]= 9.993093e-06;
    r[80][X]= 6.019071e-01; r[80][Y]= 9.667958e-01; r[80][Z]= 6.692123e-01; q[ 80]= 9.993093e-06;
    r[81][X]= 3.810881e-01; r[81][Y]= 4.726975e-01; r[81][Z]= 6.629564e-01; q[ 81]= 9.993093e-06;
    r[82][X]= 3.421055e-01; r[82][Y]= 1.646141e-01; r[82][Z]= 6.295140e-02; q[ 82]= 9.993093e-06;
    r[83][X]= 3.023214e-01; r[83][Y]= 7.843910e-01; r[83][Z]= 9.475422e-01; q[ 83]= 9.993093e-06;
    r[84][X]= 4.480307e-01; r[84][Y]= 8.482220e-02; r[84][Z]= 6.822816e-01; q[ 84]= 9.993093e-06;
    r[85][X]= 4.176568e-01; r[85][Y]= 2.316016e-01; r[85][Z]= 8.901177e-01; q[ 85]= 9.993093e-06;
    r[86][X]= 3.730756e-01; r[86][Y]= 5.283381e-01; r[86][Z]= 8.008147e-01; q[ 86]= 9.993093e-06;
    r[87][X]= 2.291990e-02; r[87][Y]= 2.708348e-01; r[87][Z]= 9.259786e-01; q[ 87]= 9.993093e-06;
    r[88][X]= 5.079385e-01; r[88][Y]= 4.242180e-02; r[88][Z]= 9.507989e-01; q[ 88]= 9.993093e-06;
    r[89][X]= 9.748271e-01; r[89][Y]= 3.529035e-01; r[89][Z]= 5.067425e-01; q[ 89]= 9.993093e-06;
    r[90][X]= 6.707042e-01; r[90][Y]= 9.643491e-01; r[90][Z]= 1.015491e-01; q[ 90]= 9.993093e-06;
    r[91][X]= 3.758621e-01; r[91][Y]= 5.297577e-01; r[91][Z]= 2.964660e-02; q[ 91]= 9.993093e-06;
    r[92][X]= 4.975944e-01; r[92][Y]= 7.600470e-02; r[92][Z]= 4.291745e-01; q[ 92]= 9.993093e-06;
    r[93][X]= 6.153944e-01; r[93][Y]= 3.733612e-01; r[93][Z]= 6.643146e-01; q[ 93]= 9.993093e-06;
    r[94][X]= 6.774827e-01; r[94][Y]= 2.160210e-01; r[94][Z]= 5.669574e-01; q[ 94]= 9.993093e-06;
    r[95][X]= 4.360429e-01; r[95][Y]= 1.845264e-01; r[95][Z]= 7.204625e-01; q[ 95]= 9.993093e-06;
    r[96][X]= 6.455320e-01; r[96][Y]= 8.217722e-01; r[96][Z]= 3.704914e-01; q[ 96]= 9.993093e-06;
    r[97][X]= 4.833917e-01; r[97][Y]= 2.673869e-01; r[97][Z]= 8.130147e-01; q[ 97]= 9.993093e-06;
    r[98][X]= 2.856406e-01; r[98][Y]= 6.262542e-01; r[98][Z]= 1.537660e-02; q[ 98]= 9.993093e-06;
    r[99][X]= 3.676946e-01; r[99][Y]= 5.239510e-01; r[99][Z]= 3.575990e-02; q[ 99]= 9.993093e-06;
    r[100][X]= 7.467359e-01; r[100][Y]= 2.881302e-01; r[100][Z]= 7.045232e-01; q[100]= 9.993093e-06;
    r[101][X]= 2.899746e-01; r[101][Y]= 2.020350e-02; r[101][Z]= 9.806426e-01; q[101]= 9.993093e-06;
    r[102][X]= 1.615000e-01; r[102][Y]= 5.612388e-01; r[102][Z]= 5.564516e-01; q[102]= 9.993093e-06;
    r[103][X]= 6.521758e-01; r[103][Y]= 7.274274e-01; r[103][Z]= 2.439741e-01; q[103]= 9.993093e-06;
    r[104][X]= 3.009680e-02; r[104][Y]= 9.166857e-01; r[104][Z]= 8.742697e-01; q[104]= 9.993093e-06;
    r[105][X]= 2.929042e-01; r[105][Y]= 6.289206e-01; r[105][Z]= 2.708670e-02; q[105]= 9.993093e-06;
    r[106][X]= 3.356011e-01; r[106][Y]= 5.196220e-01; r[106][Z]= 8.170566e-01; q[106]= 9.993093e-06;
    r[107][X]= 6.385896e-01; r[107][Y]= 7.200064e-01; r[107][Z]= 2.629148e-01; q[107]= 9.993093e-06;
    r[108][X]= 2.419459e-01; r[108][Y]= 6.530228e-01; r[108][Z]= 9.189421e-01; q[108]= 9.993093e-06;
    r[109][X]= 4.907990e-01; r[109][Y]= 7.365370e-02; r[109][Z]= 4.241293e-01; q[109]= 9.993093e-06;
    r[110][X]= 1.622683e-01; r[110][Y]= 9.956923e-01; r[110][Z]= 9.145464e-01; q[110]= 9.993093e-06;
    r[111][X]= 8.326520e-01; r[111][Y]= 2.939376e-01; r[111][Z]= 7.140182e-01; q[111]= 9.993093e-06;
    r[112][X]= 1.887263e-01; r[112][Y]= 4.982793e-01; r[112][Z]= 6.025294e-01; q[112]= 9.993093e-06;
    r[113][X]= 5.583218e-01; r[113][Y]= 5.934413e-01; r[113][Z]= 3.698000e-03; q[113]= 9.993093e-06;
    r[114][X]= 1.299183e-01; r[114][Y]= 5.094689e-01; r[114][Z]= 1.716287e-01; q[114]= 9.993093e-06;
    r[115][X]= 4.672300e-01; r[115][Y]= 2.269067e-01; r[115][Z]= 7.892066e-01; q[115]= 9.993093e-06;
    r[116][X]= 6.050443e-01; r[116][Y]= 3.661129e-01; r[116][Z]= 6.605844e-01; q[116]= 9.993093e-06;
    r[117][X]= 5.640931e-01; r[117][Y]= 1.543770e-02; r[117][Z]= 3.851490e-02; q[117]= 9.993093e-06;
    r[118][X]= 8.721566e-01; r[118][Y]= 2.855216e-01; r[118][Z]= 7.500763e-01; q[118]= 9.993093e-06;
    r[119][X]= 9.815275e-01; r[119][Y]= 3.510772e-01; r[119][Z]= 5.031639e-01; q[119]= 9.993093e-06;
    r[120][X]= 1.625797e-01; r[120][Y]= 2.881325e-01; r[120][Z]= 2.083217e-01; q[120]= 9.993093e-06;
    r[121][X]= 2.565763e-01; r[121][Y]= 6.349358e-01; r[121][Z]= 9.798580e-01; q[121]= 9.993093e-06;
    r[122][X]= 6.395137e-01; r[122][Y]= 7.256236e-01; r[122][Z]= 2.565639e-01; q[122]= 9.993093e-06;
    r[123][X]= 1.470283e-01; r[123][Y]= 1.992442e-01; r[123][Z]= 9.215969e-01; q[123]= 9.993093e-06;
    r[124][X]= 6.218717e-01; r[124][Y]= 7.286442e-01; r[124][Z]= 2.761634e-01; q[124]= 9.993093e-06;
    r[125][X]= 5.530562e-01; r[125][Y]= 8.674165e-01; r[125][Z]= 4.477672e-01; q[125]= 9.993093e-06;
    r[126][X]= 4.972708e-01; r[126][Y]= 1.570250e-02; r[126][Z]= 3.560095e-01; q[126]= 9.993093e-06;
    r[127][X]= 4.670775e-01; r[127][Y]= 2.087175e-01; r[127][Z]= 8.113433e-01; q[127]= 9.993093e-06;
    r[128][X]= 9.715530e-02; r[128][Y]= 5.118178e-01; r[128][Z]= 4.529825e-01; q[128]= 9.993093e-06;
    r[129][X]= 4.987226e-01; r[129][Y]= 6.832300e-03; r[129][Z]= 3.617589e-01; q[129]= 9.993093e-06;
    r[130][X]= 2.292325e-01; r[130][Y]= 7.948825e-01; r[130][Z]= 6.162234e-01; q[130]= 9.993093e-06;
    r[131][X]= 1.719104e-01; r[131][Y]= 1.482067e-01; r[131][Z]= 8.955181e-01; q[131]= 9.993093e-06;
    r[132][X]= 6.280350e-01; r[132][Y]= 2.911300e-02; r[132][Z]= 2.307379e-01; q[132]= 9.993093e-06;
    r[133][X]= 9.083962e-01; r[133][Y]= 1.539924e-01; r[133][Z]= 8.650628e-01; q[133]= 9.993093e-06;
    r[134][X]= 5.972086e-01; r[134][Y]= 7.876410e-01; r[134][Z]= 2.744203e-01; q[134]= 9.993093e-06;
    r[135][X]= 2.368819e-01; r[135][Y]= 7.905487e-01; r[135][Z]= 3.931051e-01; q[135]= 9.993093e-06;
    r[136][X]= 2.650500e-03; r[136][Y]= 2.918860e-01; r[136][Z]= 8.154796e-01; q[136]= 9.993093e-06;
    r[137][X]= 7.781844e-01; r[137][Y]= 6.474049e-01; r[137][Z]= 2.908288e-01; q[137]= 9.993093e-06;
    r[138][X]= 4.108779e-01; r[138][Y]= 1.853235e-01; r[138][Z]= 8.308992e-01; q[138]= 9.993093e-06;
    r[139][X]= 5.381359e-01; r[139][Y]= 9.546368e-01; r[139][Z]= 6.567694e-01; q[139]= 9.993093e-06;
    r[140][X]= 8.599949e-01; r[140][Y]= 4.434327e-01; r[140][Z]= 1.804636e-01; q[140]= 9.993093e-06;
    r[141][X]= 4.890984e-01; r[141][Y]= 3.904202e-01; r[141][Z]= 9.291727e-01; q[141]= 9.993093e-06;
    r[142][X]= 9.382193e-01; r[142][Y]= 5.502525e-01; r[142][Z]= 4.273263e-01; q[142]= 9.993093e-06;
    r[143][X]= 6.683615e-01; r[143][Y]= 3.099863e-01; r[143][Z]= 7.031798e-01; q[143]= 9.993093e-06;
    r[144][X]= 6.230500e-02; r[144][Y]= 3.137600e-01; r[144][Z]= 3.367838e-01; q[144]= 9.993093e-06;
    r[145][X]= 6.660296e-01; r[145][Y]= 8.276648e-01; r[145][Z]= 1.959473e-01; q[145]= 9.993093e-06;
    r[146][X]= 5.155564e-01; r[146][Y]= 2.216245e-01; r[146][Z]= 1.227145e-01; q[146]= 9.993093e-06;
    r[147][X]= 3.674629e-01; r[147][Y]= 8.799456e-01; r[147][Z]= 5.073631e-01; q[147]= 9.993093e-06;
    r[148][X]= 3.158536e-01; r[148][Y]= 9.934777e-01; r[148][Z]= 5.327170e-02; q[148]= 9.993093e-06;
    r[149][X]= 9.817073e-01; r[149][Y]= 1.225943e-01; r[149][Z]= 8.448881e-01; q[149]= 9.993093e-06;
    r[150][X]= 3.533229e-01; r[150][Y]= 8.994851e-01; r[150][Z]= 2.276869e-01; q[150]= 9.993093e-06;
    r[151][X]= 3.938913e-01; r[151][Y]= 4.345825e-01; r[151][Z]= 9.868051e-01; q[151]= 9.993093e-06;
    r[152][X]= 8.100450e-01; r[152][Y]= 1.932200e-02; r[152][Z]= 8.783245e-01; q[152]= 9.993093e-06;
    r[153][X]= 4.130391e-01; r[153][Y]= 1.855438e-01; r[153][Z]= 8.282541e-01; q[153]= 9.993093e-06;
    r[154][X]= 8.060218e-01; r[154][Y]= 6.410370e-02; r[154][Z]= 9.287938e-01; q[154]= 9.993093e-06;
    r[155][X]= 7.121460e-01; r[155][Y]= 1.402387e-01; r[155][Z]= 4.517867e-01; q[155]= 9.993093e-06;
    r[156][X]= 2.869502e-01; r[156][Y]= 3.704047e-01; r[156][Z]= 5.246240e-01; q[156]= 9.993093e-06;
    r[157][X]= 2.729679e-01; r[157][Y]= 8.177800e-03; r[157][Z]= 9.832208e-01; q[157]= 9.993093e-06;
    r[158][X]= 2.959580e-01; r[158][Y]= 6.174165e-01; r[158][Z]= 4.073700e-03; q[158]= 9.993093e-06;
    r[159][X]= 7.694070e-02; r[159][Y]= 5.786474e-01; r[159][Z]= 9.419320e-01; q[159]= 9.993093e-06;
    r[160][X]= 2.637322e-01; r[160][Y]= 9.412896e-01; r[160][Z]= 1.988506e-01; q[160]= 9.993093e-06;
    r[161][X]= 2.301443e-01; r[161][Y]= 6.915857e-01; r[161][Z]= 7.713379e-01; q[161]= 9.993093e-06;
    r[162][X]= 3.966906e-01; r[162][Y]= 9.034274e-01; r[162][Z]= 5.791775e-01; q[162]= 9.993093e-06;
    r[163][X]= 3.863714e-01; r[163][Y]= 4.148750e-02; r[163][Z]= 9.069938e-01; q[163]= 9.993093e-06;
    r[164][X]= 9.356280e-01; r[164][Y]= 5.570658e-01; r[164][Z]= 4.219110e-01; q[164]= 9.993093e-06;
    r[165][X]= 6.154486e-01; r[165][Y]= 5.042370e-01; r[165][Z]= 9.735208e-01; q[165]= 9.993093e-06;
    r[166][X]= 4.938524e-01; r[166][Y]= 7.009363e-01; r[166][Z]= 5.629354e-01; q[166]= 9.993093e-06;
    r[167][X]= 2.827819e-01; r[167][Y]= 7.447601e-01; r[167][Z]= 1.712008e-01; q[167]= 9.993093e-06;
    r[168][X]= 5.446769e-01; r[168][Y]= 7.348392e-01; r[168][Z]= 3.662932e-01; q[168]= 9.993093e-06;
    r[169][X]= 1.698358e-01; r[169][Y]= 9.878158e-01; r[169][Z]= 3.886524e-01; q[169]= 9.993093e-06;
    r[170][X]= 6.327914e-01; r[170][Y]= 7.214460e-01; r[170][Z]= 2.496301e-01; q[170]= 9.993093e-06;
    r[171][X]= 3.723836e-01; r[171][Y]= 3.099490e-02; r[171][Z]= 9.203130e-02; q[171]= 9.993093e-06;
    r[172][X]= 5.931068e-01; r[172][Y]= 4.915625e-01; r[172][Z]= 5.700560e-01; q[172]= 9.993093e-06;
    r[173][X]= 5.432445e-01; r[173][Y]= 1.244113e-01; r[173][Z]= 5.486222e-01; q[173]= 9.993093e-06;
    r[174][X]= 2.226736e-01; r[174][Y]= 1.044256e-01; r[174][Z]= 4.662194e-01; q[174]= 9.993093e-06;
    r[175][X]= 2.942743e-01; r[175][Y]= 7.471143e-01; r[175][Z]= 1.857922e-01; q[175]= 9.993093e-06;
    r[176][X]= 2.351411e-01; r[176][Y]= 9.911055e-01; r[176][Z]= 9.575936e-01; q[176]= 9.993093e-06;
    r[177][X]= 4.904187e-01; r[177][Y]= 4.351109e-01; r[177][Z]= 9.300984e-01; q[177]= 9.993093e-06;
    r[178][X]= 6.412730e-02; r[178][Y]= 4.085905e-01; r[178][Z]= 8.517531e-01; q[178]= 9.993093e-06;
    r[179][X]= 5.733887e-01; r[179][Y]= 9.808827e-01; r[179][Z]= 4.849513e-01; q[179]= 9.993093e-06;
    r[180][X]= 1.713340e-01; r[180][Y]= 1.482114e-01; r[180][Z]= 8.954592e-01; q[180]= 9.993093e-06;
    r[181][X]= 5.012940e-02; r[181][Y]= 8.555209e-01; r[181][Z]= 9.577163e-01; q[181]= 9.993093e-06;
    r[182][X]= 3.059542e-01; r[182][Y]= 8.556743e-01; r[182][Z]= 2.138368e-01; q[182]= 9.993093e-06;
    r[183][X]= 6.407344e-01; r[183][Y]= 7.260650e-01; r[183][Z]= 2.560783e-01; q[183]= 9.993093e-06;
    r[184][X]= 2.960440e-01; r[184][Y]= 5.241426e-01; r[184][Z]= 8.570910e-01; q[184]= 9.993093e-06;
    r[185][X]= 1.681874e-01; r[185][Y]= 9.725490e-01; r[185][Z]= 3.739514e-01; q[185]= 9.993093e-06;
    r[186][X]= 4.664148e-01; r[186][Y]= 2.274950e-01; r[186][Z]= 7.835506e-01; q[186]= 9.993093e-06;
    r[187][X]= 1.935032e-01; r[187][Y]= 8.027849e-01; r[187][Z]= 9.587369e-01; q[187]= 9.993093e-06;
    r[188][X]= 8.518451e-01; r[188][Y]= 6.317273e-01; r[188][Z]= 2.702910e-01; q[188]= 9.993093e-06;
    r[189][X]= 7.939434e-01; r[189][Y]= 5.566740e-01; r[189][Z]= 4.511294e-01; q[189]= 9.993093e-06;
    r[190][X]= 3.503468e-01; r[190][Y]= 9.795095e-01; r[190][Z]= 5.918723e-01; q[190]= 9.993093e-06;
    r[191][X]= 6.554078e-01; r[191][Y]= 4.511187e-01; r[191][Z]= 9.388715e-01; q[191]= 9.993093e-06;
    r[192][X]= 4.909605e-01; r[192][Y]= 1.601743e-01; r[192][Z]= 7.312521e-01; q[192]= 9.993093e-06;
    r[193][X]= 3.067834e-01; r[193][Y]= 2.126336e-01; r[193][Z]= 9.278500e-03; q[193]= 9.993093e-06;
    r[194][X]= 9.672149e-01; r[194][Y]= 3.494468e-01; r[194][Z]= 6.679157e-01; q[194]= 9.993093e-06;
    r[195][X]= 7.734737e-01; r[195][Y]= 6.481692e-01; r[195][Z]= 2.909202e-01; q[195]= 9.993093e-06;
    r[196][X]= 2.815826e-01; r[196][Y]= 6.479193e-01; r[196][Z]= 3.046160e-02; q[196]= 9.993093e-06;
    r[197][X]= 7.853350e-02; r[197][Y]= 5.691187e-01; r[197][Z]= 9.141507e-01; q[197]= 9.993093e-06;
    r[198][X]= 2.358938e-01; r[198][Y]= 7.832672e-01; r[198][Z]= 4.072115e-01; q[198]= 9.993093e-06;
    r[199][X]= 6.127498e-01; r[199][Y]= 1.454210e-02; r[199][Z]= 1.320850e-02; q[199]= 9.993093e-06;
    r[200][X]= 7.869783e-01; r[200][Y]= 6.553801e-01; r[200][Z]= 2.416570e-01; q[200]= 9.993093e-06;
    r[201][X]= 8.674916e-01; r[201][Y]= 7.445623e-01; r[201][Z]= 6.741449e-01; q[201]= 9.993093e-06;
    r[202][X]= 4.021716e-01; r[202][Y]= 1.917640e-01; r[202][Z]= 9.382909e-01; q[202]= 9.993093e-06;
    r[203][X]= 1.717080e-01; r[203][Y]= 1.498474e-01; r[203][Z]= 8.959584e-01; q[203]= 9.993093e-06;
    r[204][X]= 8.413390e-02; r[204][Y]= 9.493770e-02; r[204][Z]= 3.750051e-01; q[204]= 9.993093e-06;
    r[205][X]= 4.836147e-01; r[205][Y]= 4.082118e-01; r[205][Z]= 9.349136e-01; q[205]= 9.993093e-06;
    r[206][X]= 6.718129e-01; r[206][Y]= 1.897454e-01; r[206][Z]= 8.861083e-01; q[206]= 9.993093e-06;
    r[207][X]= 2.945007e-01; r[207][Y]= 7.132446e-01; r[207][Z]= 9.952393e-01; q[207]= 9.993093e-06;
    r[208][X]= 6.630270e-01; r[208][Y]= 3.073893e-01; r[208][Z]= 7.020429e-01; q[208]= 9.993093e-06;
    r[209][X]= 5.789060e-01; r[209][Y]= 3.618169e-01; r[209][Z]= 6.938101e-01; q[209]= 9.993093e-06;
    r[210][X]= 6.409934e-01; r[210][Y]= 7.358271e-01; r[210][Z]= 2.540722e-01; q[210]= 9.993093e-06;
    r[211][X]= 6.739867e-01; r[211][Y]= 8.137075e-01; r[211][Z]= 2.472267e-01; q[211]= 9.993093e-06;
    r[212][X]= 2.212572e-01; r[212][Y]= 2.255680e-01; r[212][Z]= 1.003981e-01; q[212]= 9.993093e-06;
    r[213][X]= 2.540437e-01; r[213][Y]= 7.536200e-01; r[213][Z]= 8.445597e-01; q[213]= 9.993093e-06;
    r[214][X]= 4.375150e-02; r[214][Y]= 9.418944e-01; r[214][Z]= 8.099751e-01; q[214]= 9.993093e-06;
    r[215][X]= 2.543909e-01; r[215][Y]= 6.389982e-01; r[215][Z]= 9.772842e-01; q[215]= 9.993093e-06;
    r[216][X]= 1.663617e-01; r[216][Y]= 1.516279e-01; r[216][Z]= 8.938812e-01; q[216]= 9.993093e-06;
    r[217][X]= 4.358935e-01; r[217][Y]= 9.724421e-01; r[217][Z]= 6.164103e-01; q[217]= 9.993093e-06;
    r[218][X]= 8.969336e-01; r[218][Y]= 6.479022e-01; r[218][Z]= 6.601480e-02; q[218]= 9.993093e-06;
    r[219][X]= 3.433049e-01; r[219][Y]= 1.496022e-01; r[219][Z]= 7.816526e-01; q[219]= 9.993093e-06;
    r[220][X]= 1.285045e-01; r[220][Y]= 9.609555e-01; r[220][Z]= 6.623552e-01; q[220]= 9.993093e-06;
    r[221][X]= 6.370185e-01; r[221][Y]= 4.432269e-01; r[221][Z]= 4.122740e-02; q[221]= 9.993093e-06;
    r[222][X]= 6.331516e-01; r[222][Y]= 7.159289e-01; r[222][Z]= 2.654562e-01; q[222]= 9.993093e-06;
    r[223][X]= 6.699865e-01; r[223][Y]= 9.931264e-01; r[223][Z]= 5.042830e-02; q[223]= 9.993093e-06;
    r[224][X]= 6.459974e-01; r[224][Y]= 4.625895e-01; r[224][Z]= 5.478487e-01; q[224]= 9.993093e-06;
    r[225][X]= 7.286598e-01; r[225][Y]= 7.467237e-01; r[225][Z]= 6.556760e-02; q[225]= 9.993093e-06;
    r[226][X]= 6.170427e-01; r[226][Y]= 2.776112e-01; r[226][Z]= 1.660960e-02; q[226]= 9.993093e-06;
    r[227][X]= 2.529837e-01; r[227][Y]= 9.895960e-01; r[227][Z]= 9.463426e-01; q[227]= 9.993093e-06;
    r[228][X]= 7.728232e-01; r[228][Y]= 7.686914e-01; r[228][Z]= 9.938379e-01; q[228]= 9.993093e-06;
    r[229][X]= 4.723388e-01; r[229][Y]= 4.086991e-01; r[229][Z]= 9.589221e-01; q[229]= 9.993093e-06;
    r[230][X]= 3.096460e-01; r[230][Y]= 5.302121e-01; r[230][Z]= 8.437358e-01; q[230]= 9.993093e-06;
    r[231][X]= 3.586270e-02; r[231][Y]= 3.116630e-02; r[231][Z]= 7.758842e-01; q[231]= 9.993093e-06;
    r[232][X]= 4.605404e-01; r[232][Y]= 2.264686e-01; r[232][Z]= 8.159375e-01; q[232]= 9.993093e-06;
    r[233][X]= 3.902618e-01; r[233][Y]= 4.651819e-01; r[233][Z]= 9.060126e-01; q[233]= 9.993093e-06;
    r[234][X]= 8.839687e-01; r[234][Y]= 8.747633e-01; r[234][Z]= 2.379610e-02; q[234]= 9.993093e-06;
    r[235][X]= 1.205499e-01; r[235][Y]= 5.451384e-01; r[235][Z]= 9.054577e-01; q[235]= 9.993093e-06;
    r[236][X]= 7.792653e-01; r[236][Y]= 3.505072e-01; r[236][Z]= 6.447998e-01; q[236]= 9.993093e-06;
    r[237][X]= 4.813815e-01; r[237][Y]= 4.073483e-01; r[237][Z]= 9.361451e-01; q[237]= 9.993093e-06;
    r[238][X]= 4.003202e-01; r[238][Y]= 1.703649e-01; r[238][Z]= 2.753905e-01; q[238]= 9.993093e-06;
    r[239][X]= 4.858866e-01; r[239][Y]= 4.145169e-01; r[239][Z]= 9.367168e-01; q[239]= 9.993093e-06;
    r[240][X]= 3.810247e-01; r[240][Y]= 9.008560e-02; r[240][Z]= 2.952515e-01; q[240]= 9.993093e-06;
    r[241][X]= 2.449854e-01; r[241][Y]= 9.884974e-01; r[241][Z]= 9.687867e-01; q[241]= 9.993093e-06;
    r[242][X]= 3.180296e-01; r[242][Y]= 2.327368e-01; r[242][Z]= 1.531173e-01; q[242]= 9.993093e-06;
    r[243][X]= 6.175902e-01; r[243][Y]= 5.287101e-01; r[243][Z]= 5.575691e-01; q[243]= 9.993093e-06;
    r[244][X]= 7.902868e-01; r[244][Y]= 3.953765e-01; r[244][Z]= 4.809434e-01; q[244]= 9.993093e-06;
    r[245][X]= 2.539813e-01; r[245][Y]= 1.726044e-01; r[245][Z]= 4.742949e-01; q[245]= 9.993093e-06;
    r[246][X]= 2.544396e-01; r[246][Y]= 6.417211e-01; r[246][Z]= 9.697039e-01; q[246]= 9.993093e-06;
    r[247][X]= 3.747960e-01; r[247][Y]= 3.902130e-02; r[247][Z]= 1.192820e-01; q[247]= 9.993093e-06;
    r[248][X]= 2.873313e-01; r[248][Y]= 4.900452e-01; r[248][Z]= 3.974390e-02; q[248]= 9.993093e-06;
    r[249][X]= 4.868475e-01; r[249][Y]= 4.095014e-01; r[249][Z]= 9.354048e-01; q[249]= 9.993093e-06;
    r[250][X]= 9.713837e-01; r[250][Y]= 3.485178e-01; r[250][Z]= 5.016862e-01; q[250]= 9.993093e-06;
    r[251][X]= 3.373298e-01; r[251][Y]= 2.799238e-01; r[251][Z]= 6.232292e-01; q[251]= 9.993093e-06;
    r[252][X]= 8.650503e-01; r[252][Y]= 2.402334e-01; r[252][Z]= 3.410007e-01; q[252]= 9.993093e-06;
    r[253][X]= 5.270839e-01; r[253][Y]= 7.455215e-01; r[253][Z]= 4.059568e-01; q[253]= 9.993093e-06;
    r[254][X]= 4.755424e-01; r[254][Y]= 4.111226e-01; r[254][Z]= 9.261002e-01; q[254]= 9.993093e-06;
    r[255][X]= 6.194311e-01; r[255][Y]= 4.428910e-01; r[255][Z]= 6.098368e-01; q[255]= 9.993093e-06;
    r[256][X]= 1.690718e-01; r[256][Y]= 9.717770e-01; r[256][Z]= 3.738027e-01; q[256]= 9.993093e-06;
    r[257][X]= 4.227319e-01; r[257][Y]= 5.465376e-01; r[257][Z]= 3.537850e-02; q[257]= 9.993093e-06;
    r[258][X]= 3.557144e-01; r[258][Y]= 2.001384e-01; r[258][Z]= 6.401714e-01; q[258]= 9.993093e-06;
    r[259][X]= 4.931200e-01; r[259][Y]= 4.095510e-01; r[259][Z]= 9.231425e-01; q[259]= 9.993093e-06;
    r[260][X]= 7.117772e-01; r[260][Y]= 6.902739e-01; r[260][Z]= 1.829436e-01; q[260]= 9.993093e-06;
    r[261][X]= 6.463240e-02; r[261][Y]= 8.198007e-01; r[261][Z]= 9.310350e-01; q[261]= 9.993093e-06;
    r[262][X]= 2.463758e-01; r[262][Y]= 7.901024e-01; r[262][Z]= 3.987949e-01; q[262]= 9.993093e-06;
    r[263][X]= 2.615353e-01; r[263][Y]= 7.731661e-01; r[263][Z]= 8.504039e-01; q[263]= 9.993093e-06;
    r[264][X]= 2.842541e-01; r[264][Y]= 1.261420e-02; r[264][Z]= 9.924714e-01; q[264]= 9.993093e-06;
    r[265][X]= 2.350849e-01; r[265][Y]= 9.900961e-01; r[265][Z]= 9.616180e-01; q[265]= 9.993093e-06;
    r[266][X]= 9.690383e-01; r[266][Y]= 5.680671e-01; r[266][Z]= 4.175588e-01; q[266]= 9.993093e-06;
    r[267][X]= 2.923423e-01; r[267][Y]= 1.735015e-01; r[267][Z]= 8.860862e-01; q[267]= 9.993093e-06;
    r[268][X]= 1.841845e-01; r[268][Y]= 8.276337e-01; r[268][Z]= 8.863673e-01; q[268]= 9.993093e-06;
    r[269][X]= 2.210020e-02; r[269][Y]= 5.813012e-01; r[269][Z]= 4.295164e-01; q[269]= 9.993093e-06;
    r[270][X]= 2.593593e-01; r[270][Y]= 1.774219e-01; r[270][Z]= 4.251579e-01; q[270]= 9.993093e-06;
    r[271][X]= 4.949677e-01; r[271][Y]= 7.805596e-01; r[271][Z]= 4.595995e-01; q[271]= 9.993093e-06;
    r[272][X]= 4.628237e-01; r[272][Y]= 5.806021e-01; r[272][Z]= 9.052190e-02; q[272]= 9.993093e-06;
    r[273][X]= 4.388238e-01; r[273][Y]= 2.154238e-01; r[273][Z]= 9.277748e-01; q[273]= 9.993093e-06;
    r[274][X]= 9.405030e-02; r[274][Y]= 9.468625e-01; r[274][Z]= 2.250895e-01; q[274]= 9.993093e-06;
    r[275][X]= 7.244050e-02; r[275][Y]= 9.327710e-01; r[275][Z]= 6.261965e-01; q[275]= 9.993093e-06;
    r[276][X]= 8.670530e-02; r[276][Y]= 8.479211e-01; r[276][Z]= 5.621278e-01; q[276]= 9.993093e-06;
    r[277][X]= 1.982148e-01; r[277][Y]= 6.697485e-01; r[277][Z]= 3.651273e-01; q[277]= 9.993093e-06;
    r[278][X]= 3.046903e-01; r[278][Y]= 8.629267e-01; r[278][Z]= 1.867001e-01; q[278]= 9.993093e-06;
    r[279][X]= 2.364737e-01; r[279][Y]= 9.915632e-01; r[279][Z]= 9.625573e-01; q[279]= 9.993093e-06;
    r[280][X]= 9.534796e-01; r[280][Y]= 4.448758e-01; r[280][Z]= 3.793409e-01; q[280]= 9.993093e-06;
    r[281][X]= 9.928797e-01; r[281][Y]= 3.381901e-01; r[281][Z]= 4.699331e-01; q[281]= 9.993093e-06;
    r[282][X]= 4.454979e-01; r[282][Y]= 2.847480e-01; r[282][Z]= 8.955446e-01; q[282]= 9.993093e-06;
    r[283][X]= 4.838629e-01; r[283][Y]= 4.195193e-01; r[283][Z]= 9.392855e-01; q[283]= 9.993093e-06;
    r[284][X]= 7.423236e-01; r[284][Y]= 1.932730e-02; r[284][Z]= 7.294261e-01; q[284]= 9.993093e-06;
    r[285][X]= 5.832959e-01; r[285][Y]= 6.133283e-01; r[285][Z]= 3.026630e-01; q[285]= 9.993093e-06;
    r[286][X]= 4.465433e-01; r[286][Y]= 8.376560e-02; r[286][Z]= 6.797591e-01; q[286]= 9.993093e-06;
    r[287][X]= 1.740989e-01; r[287][Y]= 1.410923e-01; r[287][Z]= 9.036877e-01; q[287]= 9.993093e-06;
    r[288][X]= 2.032621e-01; r[288][Y]= 8.479360e-02; r[288][Z]= 9.130916e-01; q[288]= 9.993093e-06;
    r[289][X]= 1.658714e-01; r[289][Y]= 5.684340e-01; r[289][Z]= 5.107333e-01; q[289]= 9.993093e-06;
    r[290][X]= 8.784999e-01; r[290][Y]= 1.025457e-01; r[290][Z]= 8.769812e-01; q[290]= 9.993093e-06;
    r[291][X]= 3.162612e-01; r[291][Y]= 9.933904e-01; r[291][Z]= 5.268250e-02; q[291]= 9.993093e-06;
    r[292][X]= 6.024554e-01; r[292][Y]= 2.934860e-02; r[292][Z]= 9.962364e-01; q[292]= 9.993093e-06;
    r[293][X]= 5.516308e-01; r[293][Y]= 7.345093e-01; r[293][Z]= 3.645935e-01; q[293]= 9.993093e-06;
    r[294][X]= 7.612352e-01; r[294][Y]= 4.031502e-01; r[294][Z]= 7.059930e-02; q[294]= 9.993093e-06;
    r[295][X]= 7.431577e-01; r[295][Y]= 6.929313e-01; r[295][Z]= 2.153398e-01; q[295]= 9.993093e-06;
    r[296][X]= 4.132090e-02; r[296][Y]= 9.769734e-01; r[296][Z]= 8.306987e-01; q[296]= 9.993093e-06;
    r[297][X]= 1.970274e-01; r[297][Y]= 6.682970e-01; r[297][Z]= 3.378919e-01; q[297]= 9.993093e-06;
    r[298][X]= 6.312134e-01; r[298][Y]= 4.286448e-01; r[298][Z]= 6.180144e-01; q[298]= 9.993093e-06;
    r[299][X]= 8.315464e-01; r[299][Y]= 2.977569e-01; r[299][Z]= 7.069453e-01; q[299]= 9.993093e-06;
    r[300][X]= 6.409240e-01; r[300][Y]= 6.853335e-01; r[300][Z]= 9.733689e-01; q[300]= 9.993093e-06;
    r[301][X]= 5.933629e-01; r[301][Y]= 4.167516e-01; r[301][Z]= 8.435750e-02; q[301]= 9.993093e-06;
    r[302][X]= 5.221529e-01; r[302][Y]= 5.584600e-03; r[302][Z]= 3.664298e-01; q[302]= 9.993093e-06;
    r[303][X]= 4.510212e-01; r[303][Y]= 6.970100e-02; r[303][Z]= 6.745322e-01; q[303]= 9.993093e-06;
    r[304][X]= 1.661761e-01; r[304][Y]= 1.542092e-01; r[304][Z]= 8.935586e-01; q[304]= 9.993093e-06;
    r[305][X]= 2.834428e-01; r[305][Y]= 1.203330e-02; r[305][Z]= 9.888049e-01; q[305]= 9.993093e-06;
    r[306][X]= 3.970287e-01; r[306][Y]= 1.834573e-01; r[306][Z]= 9.420498e-01; q[306]= 9.993093e-06;
    r[307][X]= 1.141202e-01; r[307][Y]= 7.862020e-01; r[307][Z]= 9.664140e-01; q[307]= 9.993093e-06;
    r[308][X]= 1.088802e-01; r[308][Y]= 5.506056e-01; r[308][Z]= 9.258125e-01; q[308]= 9.993093e-06;
    r[309][X]= 2.391168e-01; r[309][Y]= 7.350962e-01; r[309][Z]= 9.432735e-01; q[309]= 9.993093e-06;
    r[310][X]= 5.762252e-01; r[310][Y]= 3.048869e-01; r[310][Z]= 7.293723e-01; q[310]= 9.993093e-06;
    r[311][X]= 2.140300e-03; r[311][Y]= 5.819753e-01; r[311][Z]= 4.128774e-01; q[311]= 9.993093e-06;
    r[312][X]= 4.519406e-01; r[312][Y]= 2.136002e-01; r[312][Z]= 7.547529e-01; q[312]= 9.993093e-06;
    r[313][X]= 3.770624e-01; r[313][Y]= 6.475763e-01; r[313][Z]= 8.181164e-01; q[313]= 9.993093e-06;
    r[314][X]= 4.565669e-01; r[314][Y]= 9.421645e-01; r[314][Z]= 6.105129e-01; q[314]= 9.993093e-06;
    r[315][X]= 5.347706e-01; r[315][Y]= 1.299227e-01; r[315][Z]= 6.313245e-01; q[315]= 9.993093e-06;
    r[316][X]= 4.936303e-01; r[316][Y]= 7.808295e-01; r[316][Z]= 4.502500e-01; q[316]= 9.993093e-06;
    r[317][X]= 4.697134e-01; r[317][Y]= 9.521661e-01; r[317][Z]= 6.200533e-01; q[317]= 9.993093e-06;
    r[318][X]= 1.700800e-01; r[318][Y]= 9.768800e-01; r[318][Z]= 3.831212e-01; q[318]= 9.993093e-06;
    r[319][X]= 4.484549e-01; r[319][Y]= 6.365104e-01; r[319][Z]= 6.767560e-02; q[319]= 9.993093e-06;
    r[320][X]= 9.181723e-01; r[320][Y]= 7.733048e-01; r[320][Z]= 6.288640e-01; q[320]= 9.993093e-06;
    r[321][X]= 2.351270e-01; r[321][Y]= 6.982184e-01; r[321][Z]= 8.515429e-01; q[321]= 9.993093e-06;
    r[322][X]= 4.787458e-01; r[322][Y]= 4.013712e-01; r[322][Z]= 9.342588e-01; q[322]= 9.993093e-06;
    r[323][X]= 4.001097e-01; r[323][Y]= 1.882590e-02; r[323][Z]= 2.168754e-01; q[323]= 9.993093e-06;
    r[324][X]= 2.902373e-01; r[324][Y]= 6.300303e-01; r[324][Z]= 2.046990e-02; q[324]= 9.993093e-06;
    r[325][X]= 9.801521e-01; r[325][Y]= 5.677516e-01; r[325][Z]= 4.035305e-01; q[325]= 9.993093e-06;
    r[326][X]= 6.496070e-02; r[326][Y]= 2.662055e-01; r[326][Z]= 4.990308e-01; q[326]= 9.993093e-06;
    r[327][X]= 8.969371e-01; r[327][Y]= 5.255881e-01; r[327][Z]= 4.382167e-01; q[327]= 9.993093e-06;
    r[328][X]= 6.369236e-01; r[328][Y]= 9.591399e-01; r[328][Z]= 6.910638e-01; q[328]= 9.993093e-06;
    r[329][X]= 6.093814e-01; r[329][Y]= 4.725995e-01; r[329][Z]= 5.937109e-01; q[329]= 9.993093e-06;
    r[330][X]= 2.262190e-01; r[330][Y]= 1.381464e-01; r[330][Z]= 9.003524e-01; q[330]= 9.993093e-06;
    r[331][X]= 2.431488e-01; r[331][Y]= 6.749640e-01; r[331][Z]= 8.608714e-01; q[331]= 9.993093e-06;
    r[332][X]= 9.703014e-01; r[332][Y]= 2.977927e-01; r[332][Z]= 7.086155e-01; q[332]= 9.993093e-06;
    r[333][X]= 5.450343e-01; r[333][Y]= 9.134600e-01; r[333][Z]= 3.616570e-01; q[333]= 9.993093e-06;
    r[334][X]= 6.050059e-01; r[334][Y]= 4.656065e-01; r[334][Z]= 5.915897e-01; q[334]= 9.993093e-06;
    r[335][X]= 3.094643e-01; r[335][Y]= 9.922765e-01; r[335][Z]= 5.262790e-02; q[335]= 9.993093e-06;
    r[336][X]= 4.400231e-01; r[336][Y]= 1.574770e-02; r[336][Z]= 2.763208e-01; q[336]= 9.993093e-06;
    r[337][X]= 7.147642e-01; r[337][Y]= 4.638974e-01; r[337][Z]= 9.711105e-01; q[337]= 9.993093e-06;
    r[338][X]= 5.499625e-01; r[338][Y]= 7.355512e-01; r[338][Z]= 3.650767e-01; q[338]= 9.993093e-06;
    r[339][X]= 1.187146e-01; r[339][Y]= 2.162690e-02; r[339][Z]= 4.467074e-01; q[339]= 9.993093e-06;
    r[340][X]= 1.812806e-01; r[340][Y]= 5.645619e-01; r[340][Z]= 9.957389e-01; q[340]= 9.993093e-06;
    r[341][X]= 5.308051e-01; r[341][Y]= 5.130930e-02; r[341][Z]= 6.983950e-01; q[341]= 9.993093e-06;
    r[342][X]= 6.631892e-01; r[342][Y]= 6.952277e-01; r[342][Z]= 2.432669e-01; q[342]= 9.993093e-06;
    r[343][X]= 1.765489e-01; r[343][Y]= 9.053465e-01; r[343][Z]= 6.471056e-01; q[343]= 9.993093e-06;
    r[344][X]= 3.338318e-01; r[344][Y]= 9.835421e-01; r[344][Z]= 8.878057e-01; q[344]= 9.993093e-06;
    r[345][X]= 6.349998e-01; r[345][Y]= 7.347705e-01; r[345][Z]= 2.519209e-01; q[345]= 9.993093e-06;
    r[346][X]= 8.478030e-01; r[346][Y]= 9.851480e-02; r[346][Z]= 9.050156e-01; q[346]= 9.993093e-06;
    r[347][X]= 3.139796e-01; r[347][Y]= 9.754692e-01; r[347][Z]= 9.223179e-01; q[347]= 9.993093e-06;
    r[348][X]= 6.369967e-01; r[348][Y]= 7.167122e-01; r[348][Z]= 2.577933e-01; q[348]= 9.993093e-06;
    r[349][X]= 4.612772e-01; r[349][Y]= 4.161525e-01; r[349][Z]= 8.877579e-01; q[349]= 9.993093e-06;
    r[350][X]= 1.220926e-01; r[350][Y]= 3.283005e-01; r[350][Z]= 2.373917e-01; q[350]= 9.993093e-06;
    r[351][X]= 1.794536e-01; r[351][Y]= 1.444035e-01; r[351][Z]= 9.003375e-01; q[351]= 9.993093e-06;
    r[352][X]= 2.680507e-01; r[352][Y]= 7.696674e-01; r[352][Z]= 3.009475e-01; q[352]= 9.993093e-06;
    r[353][X]= 4.238382e-01; r[353][Y]= 2.128448e-01; r[353][Z]= 9.408362e-01; q[353]= 9.993093e-06;
    r[354][X]= 2.896751e-01; r[354][Y]= 6.348923e-01; r[354][Z]= 1.792360e-02; q[354]= 9.993093e-06;
    r[355][X]= 5.529032e-01; r[355][Y]= 8.584417e-01; r[355][Z]= 1.717870e-01; q[355]= 9.993093e-06;
    r[356][X]= 2.782094e-01; r[356][Y]= 7.770963e-01; r[356][Z]= 7.759170e-01; q[356]= 9.993093e-06;
    r[357][X]= 9.603123e-01; r[357][Y]= 5.467005e-01; r[357][Z]= 4.129080e-01; q[357]= 9.993093e-06;
    r[358][X]= 1.558380e-02; r[358][Y]= 6.946238e-01; r[358][Z]= 1.569400e-03; q[358]= 9.993093e-06;
    r[359][X]= 8.593835e-01; r[359][Y]= 6.524087e-01; r[359][Z]= 5.888396e-01; q[359]= 9.993093e-06;
    r[360][X]= 9.200705e-01; r[360][Y]= 3.777825e-01; r[360][Z]= 2.156604e-01; q[360]= 9.993093e-06;
    r[361][X]= 3.809146e-01; r[361][Y]= 3.532410e-01; r[361][Z]= 6.581032e-01; q[361]= 9.993093e-06;
    r[362][X]= 8.875729e-01; r[362][Y]= 1.983107e-01; r[362][Z]= 7.508000e-04; q[362]= 9.993093e-06;
    r[363][X]= 3.754472e-01; r[363][Y]= 2.966760e-02; r[363][Z]= 2.008778e-01; q[363]= 9.993093e-06;
    r[364][X]= 4.653650e-01; r[364][Y]= 1.009549e-01; r[364][Z]= 9.640924e-01; q[364]= 9.993093e-06;
    r[365][X]= 3.614776e-01; r[365][Y]= 5.220592e-01; r[365][Z]= 8.147396e-01; q[365]= 9.993093e-06;
    r[366][X]= 6.594881e-01; r[366][Y]= 9.441490e-02; r[366][Z]= 6.614670e-01; q[366]= 9.993093e-06;
    r[367][X]= 3.830157e-01; r[367][Y]= 4.822315e-01; r[367][Z]= 9.161129e-01; q[367]= 9.993093e-06;
    r[368][X]= 4.128051e-01; r[368][Y]= 9.026970e-02; r[368][Z]= 3.177154e-01; q[368]= 9.993093e-06;
    r[369][X]= 3.762979e-01; r[369][Y]= 7.810635e-01; r[369][Z]= 1.197400e-01; q[369]= 9.993093e-06;
    r[370][X]= 1.837093e-01; r[370][Y]= 9.190315e-01; r[370][Z]= 4.212947e-01; q[370]= 9.993093e-06;
    r[371][X]= 2.354486e-01; r[371][Y]= 6.181627e-01; r[371][Z]= 1.579579e-01; q[371]= 9.993093e-06;
    r[372][X]= 5.651614e-01; r[372][Y]= 9.760027e-01; r[372][Z]= 2.391389e-01; q[372]= 9.993093e-06;
    r[373][X]= 4.253473e-01; r[373][Y]= 9.198382e-01; r[373][Z]= 2.840915e-01; q[373]= 9.993093e-06;
    r[374][X]= 8.292498e-01; r[374][Y]= 5.974528e-01; r[374][Z]= 3.467540e-01; q[374]= 9.993093e-06;
    r[375][X]= 3.999800e-01; r[375][Y]= 2.906330e-02; r[375][Z]= 2.296183e-01; q[375]= 9.993093e-06;
    r[376][X]= 9.074204e-01; r[376][Y]= 1.712184e-01; r[376][Z]= 8.466388e-01; q[376]= 9.993093e-06;
    r[377][X]= 7.776713e-01; r[377][Y]= 7.199529e-01; r[377][Z]= 4.833861e-01; q[377]= 9.993093e-06;
    r[378][X]= 3.942914e-01; r[378][Y]= 5.145888e-01; r[378][Z]= 9.962784e-01; q[378]= 9.993093e-06;
    r[379][X]= 4.685802e-01; r[379][Y]= 5.682405e-01; r[379][Z]= 7.395732e-01; q[379]= 9.993093e-06;
    r[380][X]= 8.906200e-03; r[380][Y]= 2.773881e-01; r[380][Z]= 8.253802e-01; q[380]= 9.993093e-06;
    r[381][X]= 2.668028e-01; r[381][Y]= 2.127092e-01; r[381][Z]= 3.170343e-01; q[381]= 9.993093e-06;
    r[382][X]= 5.993136e-01; r[382][Y]= 7.487454e-01; r[382][Z]= 2.767052e-01; q[382]= 9.993093e-06;
    r[383][X]= 6.525746e-01; r[383][Y]= 6.778380e-02; r[383][Z]= 9.739673e-01; q[383]= 9.993093e-06;
    r[384][X]= 3.691860e-02; r[384][Y]= 9.881868e-01; r[384][Z]= 8.241665e-01; q[384]= 9.993093e-06;
    r[385][X]= 2.535318e-01; r[385][Y]= 6.745943e-01; r[385][Z]= 7.564900e-01; q[385]= 9.993093e-06;
    r[386][X]= 6.536385e-01; r[386][Y]= 9.965687e-01; r[386][Z]= 8.599800e-03; q[386]= 9.993093e-06;
    r[387][X]= 6.631698e-01; r[387][Y]= 1.319420e-02; r[387][Z]= 3.616641e-01; q[387]= 9.993093e-06;
    r[388][X]= 3.938669e-01; r[388][Y]= 2.254012e-01; r[388][Z]= 6.834748e-01; q[388]= 9.993093e-06;
    r[389][X]= 6.356969e-01; r[389][Y]= 9.737976e-01; r[389][Z]= 1.205844e-01; q[389]= 9.993093e-06;
    r[390][X]= 2.436132e-01; r[390][Y]= 9.646862e-01; r[390][Z]= 4.171784e-01; q[390]= 9.993093e-06;
    r[391][X]= 9.276844e-01; r[391][Y]= 8.683637e-01; r[391][Z]= 3.482830e-02; q[391]= 9.993093e-06;
    r[392][X]= 9.657989e-01; r[392][Y]= 3.537790e-01; r[392][Z]= 6.709154e-01; q[392]= 9.993093e-06;
    r[393][X]= 6.397859e-01; r[393][Y]= 6.636408e-01; r[393][Z]= 3.930808e-01; q[393]= 9.993093e-06;
    r[394][X]= 2.841950e-01; r[394][Y]= 1.109820e-02; r[394][Z]= 9.898307e-01; q[394]= 9.993093e-06;
    r[395][X]= 3.338357e-01; r[395][Y]= 1.681920e-01; r[395][Z]= 8.913395e-01; q[395]= 9.993093e-06;
    r[396][X]= 5.116090e-02; r[396][Y]= 2.939200e-03; r[396][Z]= 8.867509e-01; q[396]= 9.993093e-06;
    r[397][X]= 4.708371e-01; r[397][Y]= 9.547396e-01; r[397][Z]= 6.182435e-01; q[397]= 9.993093e-06;
    r[398][X]= 4.517283e-01; r[398][Y]= 6.968310e-02; r[398][Z]= 6.857889e-01; q[398]= 9.993093e-06;
    r[399][X]= 3.916615e-01; r[399][Y]= 2.249520e-01; r[399][Z]= 6.914474e-01; q[399]= 9.993093e-06;
    r[400][X]= 4.469028e-01; r[400][Y]= 9.287069e-01; r[400][Z]= 5.754095e-01; q[400]= 9.993093e-06;
    r[401][X]= 8.441570e-02; r[401][Y]= 9.546136e-01; r[401][Z]= 7.938960e-01; q[401]= 9.993093e-06;
    r[402][X]= 9.154634e-01; r[402][Y]= 7.059277e-01; r[402][Z]= 7.829363e-01; q[402]= 9.993093e-06;
    r[403][X]= 1.341534e-01; r[403][Y]= 1.381927e-01; r[403][Z]= 5.573066e-01; q[403]= 9.993093e-06;
    r[404][X]= 2.935473e-01; r[404][Y]= 2.351480e-02; r[404][Z]= 9.920232e-01; q[404]= 9.993093e-06;
    r[405][X]= 4.900458e-01; r[405][Y]= 5.826397e-01; r[405][Z]= 7.041600e-01; q[405]= 9.993093e-06;
    r[406][X]= 3.940933e-01; r[406][Y]= 2.253583e-01; r[406][Z]= 6.891149e-01; q[406]= 9.993093e-06;
    r[407][X]= 1.733536e-01; r[407][Y]= 1.491691e-01; r[407][Z]= 9.140937e-01; q[407]= 9.993093e-06;
    r[408][X]= 1.387469e-01; r[408][Y]= 3.991855e-01; r[408][Z]= 6.656632e-01; q[408]= 9.993093e-06;
    r[409][X]= 6.329752e-01; r[409][Y]= 7.395154e-01; r[409][Z]= 2.462009e-01; q[409]= 9.993093e-06;
    r[410][X]= 6.126582e-01; r[410][Y]= 3.688115e-01; r[410][Z]= 6.613964e-01; q[410]= 9.993093e-06;
    r[411][X]= 1.026275e-01; r[411][Y]= 2.636577e-01; r[411][Z]= 9.824040e-01; q[411]= 9.993093e-06;
    r[412][X]= 5.979140e-01; r[412][Y]= 3.688913e-01; r[412][Z]= 9.458058e-01; q[412]= 9.993093e-06;
    r[413][X]= 1.632913e-01; r[413][Y]= 5.431284e-01; r[413][Z]= 6.897200e-03; q[413]= 9.993093e-06;
    r[414][X]= 7.047238e-01; r[414][Y]= 3.997731e-01; r[414][Z]= 6.034167e-01; q[414]= 9.993093e-06;
    r[415][X]= 8.239799e-01; r[415][Y]= 6.708920e-02; r[415][Z]= 9.102806e-01; q[415]= 9.993093e-06;
    r[416][X]= 5.516765e-01; r[416][Y]= 7.312488e-01; r[416][Z]= 3.635425e-01; q[416]= 9.993093e-06;
    r[417][X]= 9.481257e-01; r[417][Y]= 7.516960e-01; r[417][Z]= 7.206651e-01; q[417]= 9.993093e-06;
    r[418][X]= 3.276180e-01; r[418][Y]= 2.366602e-01; r[418][Z]= 1.614410e-01; q[418]= 9.993093e-06;
    r[419][X]= 1.260523e-01; r[419][Y]= 2.321554e-01; r[419][Z]= 7.824625e-01; q[419]= 9.993093e-06;
    r[420][X]= 6.096222e-01; r[420][Y]= 7.228681e-01; r[420][Z]= 2.872594e-01; q[420]= 9.993093e-06;
    r[421][X]= 5.680047e-01; r[421][Y]= 7.389191e-01; r[421][Z]= 3.177889e-01; q[421]= 9.993093e-06;
    r[422][X]= 6.787075e-01; r[422][Y]= 4.189831e-01; r[422][Z]= 5.583010e-01; q[422]= 9.993093e-06;
    r[423][X]= 3.557127e-01; r[423][Y]= 1.011897e-01; r[423][Z]= 9.039196e-01; q[423]= 9.993093e-06;
    r[424][X]= 7.897092e-01; r[424][Y]= 9.386740e-01; r[424][Z]= 2.204521e-01; q[424]= 9.993093e-06;
    r[425][X]= 4.733527e-01; r[425][Y]= 4.114872e-01; r[425][Z]= 9.169298e-01; q[425]= 9.993093e-06;
    r[426][X]= 4.507226e-01; r[426][Y]= 2.807021e-01; r[426][Z]= 8.755396e-01; q[426]= 9.993093e-06;
    r[427][X]= 5.424132e-01; r[427][Y]= 1.066629e-01; r[427][Z]= 5.022072e-01; q[427]= 9.993093e-06;
    r[428][X]= 3.715144e-01; r[428][Y]= 5.311499e-01; r[428][Z]= 3.619070e-02; q[428]= 9.993093e-06;
    r[429][X]= 9.694591e-01; r[429][Y]= 5.739996e-01; r[429][Z]= 9.890744e-01; q[429]= 9.993093e-06;
    r[430][X]= 4.256594e-01; r[430][Y]= 3.499372e-01; r[430][Z]= 8.813680e-01; q[430]= 9.993093e-06;
    r[431][X]= 4.788224e-01; r[431][Y]= 4.483767e-01; r[431][Z]= 8.338234e-01; q[431]= 9.993093e-06;
    r[432][X]= 8.703940e-01; r[432][Y]= 1.113770e-01; r[432][Z]= 8.845622e-01; q[432]= 9.993093e-06;
    r[433][X]= 4.870350e-01; r[433][Y]= 4.123304e-01; r[433][Z]= 9.389641e-01; q[433]= 9.993093e-06;
    r[434][X]= 7.471701e-01; r[434][Y]= 8.818106e-01; r[434][Z]= 3.135555e-01; q[434]= 9.993093e-06;
    r[435][X]= 1.115623e-01; r[435][Y]= 9.757223e-01; r[435][Z]= 3.773268e-01; q[435]= 9.993093e-06;
    r[436][X]= 7.903241e-01; r[436][Y]= 9.028322e-01; r[436][Z]= 6.709140e-01; q[436]= 9.993093e-06;
    r[437][X]= 4.758433e-01; r[437][Y]= 4.556865e-01; r[437][Z]= 8.574179e-01; q[437]= 9.993093e-06;
    r[438][X]= 1.191198e-01; r[438][Y]= 2.062275e-01; r[438][Z]= 6.788041e-01; q[438]= 9.993093e-06;
    r[439][X]= 4.621878e-01; r[439][Y]= 2.191444e-01; r[439][Z]= 7.624047e-01; q[439]= 9.993093e-06;
    r[440][X]= 1.445780e-01; r[440][Y]= 3.144762e-01; r[440][Z]= 4.304484e-01; q[440]= 9.993093e-06;
    r[441][X]= 2.712719e-01; r[441][Y]= 4.783002e-01; r[441][Z]= 2.624562e-01; q[441]= 9.993093e-06;
    r[442][X]= 4.462046e-01; r[442][Y]= 8.660921e-01; r[442][Z]= 4.303944e-01; q[442]= 9.993093e-06;
    r[443][X]= 8.952862e-01; r[443][Y]= 2.106014e-01; r[443][Z]= 6.124900e-03; q[443]= 9.993093e-06;
    r[444][X]= 9.606068e-01; r[444][Y]= 3.807855e-01; r[444][Z]= 2.200633e-01; q[444]= 9.993093e-06;
    r[445][X]= 5.852413e-01; r[445][Y]= 1.264235e-01; r[445][Z]= 2.852339e-01; q[445]= 9.993093e-06;
    r[446][X]= 6.478202e-01; r[446][Y]= 7.411938e-01; r[446][Z]= 2.632262e-01; q[446]= 9.993093e-06;
    r[447][X]= 4.848639e-01; r[447][Y]= 4.089852e-01; r[447][Z]= 9.354079e-01; q[447]= 9.993093e-06;
    r[448][X]= 2.507121e-01; r[448][Y]= 9.892093e-01; r[448][Z]= 9.684419e-01; q[448]= 9.993093e-06;
    r[449][X]= 2.692138e-01; r[449][Y]= 6.389329e-01; r[449][Z]= 9.694403e-01; q[449]= 9.993093e-06;
    r[450][X]= 2.377985e-01; r[450][Y]= 7.406127e-01; r[450][Z]= 9.129351e-01; q[450]= 9.993093e-06;
    r[451][X]= 3.386372e-01; r[451][Y]= 3.085162e-01; r[451][Z]= 6.216590e-02; q[451]= 9.993093e-06;
    r[452][X]= 7.402349e-01; r[452][Y]= 1.776560e-02; r[452][Z]= 9.479945e-01; q[452]= 9.993093e-06;
    r[453][X]= 1.737346e-01; r[453][Y]= 9.696180e-01; r[453][Z]= 3.843305e-01; q[453]= 9.993093e-06;
    r[454][X]= 9.024968e-01; r[454][Y]= 1.365204e-01; r[454][Z]= 8.475435e-01; q[454]= 9.993093e-06;
    r[455][X]= 3.954231e-01; r[455][Y]= 5.172741e-01; r[455][Z]= 2.190200e-03; q[455]= 9.993093e-06;
    r[456][X]= 3.673120e-01; r[456][Y]= 9.115603e-01; r[456][Z]= 4.667580e-02; q[456]= 9.993093e-06;
    r[457][X]= 2.363769e-01; r[457][Y]= 9.973970e-01; r[457][Z]= 9.784382e-01; q[457]= 9.993093e-06;
    r[458][X]= 6.892704e-01; r[458][Y]= 2.344152e-01; r[458][Z]= 5.649641e-01; q[458]= 9.993093e-06;
    r[459][X]= 1.977969e-01; r[459][Y]= 3.711294e-01; r[459][Z]= 2.042685e-01; q[459]= 9.993093e-06;
    r[460][X]= 3.207200e-01; r[460][Y]= 2.456476e-01; r[460][Z]= 3.150970e-02; q[460]= 9.993093e-06;
    r[461][X]= 4.682278e-01; r[461][Y]= 2.322716e-01; r[461][Z]= 7.868326e-01; q[461]= 9.993093e-06;
    r[462][X]= 5.293811e-01; r[462][Y]= 5.257190e-02; r[462][Z]= 1.281882e-01; q[462]= 9.993093e-06;
    r[463][X]= 4.814052e-01; r[463][Y]= 4.064474e-01; r[463][Z]= 9.365027e-01; q[463]= 9.993093e-06;
    r[464][X]= 3.916906e-01; r[464][Y]= 3.384480e-02; r[464][Z]= 2.363765e-01; q[464]= 9.993093e-06;
    r[465][X]= 4.690833e-01; r[465][Y]= 3.196250e-01; r[465][Z]= 8.874106e-01; q[465]= 9.993093e-06;
    r[466][X]= 4.815483e-01; r[466][Y]= 4.071391e-01; r[466][Z]= 9.273102e-01; q[466]= 9.993093e-06;
    r[467][X]= 2.313120e-01; r[467][Y]= 7.977766e-01; r[467][Z]= 4.078127e-01; q[467]= 9.993093e-06;
    r[468][X]= 4.518375e-01; r[468][Y]= 2.446049e-01; r[468][Z]= 8.532871e-01; q[468]= 9.993093e-06;
    r[469][X]= 8.653740e-01; r[469][Y]= 4.085620e-01; r[469][Z]= 1.732070e-01; q[469]= 9.993093e-06;
    r[470][X]= 9.807751e-01; r[470][Y]= 5.710414e-01; r[470][Z]= 9.844820e-01; q[470]= 9.993093e-06;
    r[471][X]= 1.671888e-01; r[471][Y]= 9.722523e-01; r[471][Z]= 3.741441e-01; q[471]= 9.993093e-06;
    r[472][X]= 2.928569e-01; r[472][Y]= 5.937928e-01; r[472][Z]= 7.312085e-01; q[472]= 9.993093e-06;
    r[473][X]= 3.936396e-01; r[473][Y]= 1.834283e-01; r[473][Z]= 5.512000e-03; q[473]= 9.993093e-06;
    r[474][X]= 1.240540e-01; r[474][Y]= 9.215065e-01; r[474][Z]= 9.428267e-01; q[474]= 9.993093e-06;
    r[475][X]= 2.856527e-01; r[475][Y]= 6.308636e-01; r[475][Z]= 1.678970e-02; q[475]= 9.993093e-06;
    r[476][X]= 9.779246e-01; r[476][Y]= 9.754175e-01; r[476][Z]= 7.265989e-01; q[476]= 9.993093e-06;
    r[477][X]= 3.108851e-01; r[477][Y]= 9.201230e-02; r[477][Z]= 7.928353e-01; q[477]= 9.993093e-06;
    r[478][X]= 2.088936e-01; r[478][Y]= 6.762553e-01; r[478][Z]= 7.177462e-01; q[478]= 9.993093e-06;
    r[479][X]= 3.225024e-01; r[479][Y]= 2.287382e-01; r[479][Z]= 1.602944e-01; q[479]= 9.993093e-06;
    r[480][X]= 3.112072e-01; r[480][Y]= 8.521408e-01; r[480][Z]= 2.205493e-01; q[480]= 9.993093e-06;
    r[481][X]= 1.626118e-01; r[481][Y]= 8.202666e-01; r[481][Z]= 5.318151e-01; q[481]= 9.993093e-06;
    r[482][X]= 2.124493e-01; r[482][Y]= 7.967520e-01; r[482][Z]= 5.908972e-01; q[482]= 9.993093e-06;
    r[483][X]= 6.889952e-01; r[483][Y]= 3.249358e-01; r[483][Z]= 5.322943e-01; q[483]= 9.993093e-06;
    r[484][X]= 2.174740e-02; r[484][Y]= 4.281430e-02; r[484][Z]= 8.776071e-01; q[484]= 9.993093e-06;
    r[485][X]= 1.692043e-01; r[485][Y]= 1.458984e-01; r[485][Z]= 8.903457e-01; q[485]= 9.993093e-06;
    r[486][X]= 5.485250e-02; r[486][Y]= 4.413439e-01; r[486][Z]= 3.369118e-01; q[486]= 9.993093e-06;
    r[487][X]= 7.530378e-01; r[487][Y]= 3.307324e-01; r[487][Z]= 5.512974e-01; q[487]= 9.993093e-06;
    r[488][X]= 4.278710e-02; r[488][Y]= 1.313136e-01; r[488][Z]= 8.952999e-01; q[488]= 9.993093e-06;
    r[489][X]= 4.030789e-01; r[489][Y]= 5.240499e-01; r[489][Z]= 9.669881e-01; q[489]= 9.993093e-06;
    r[490][X]= 7.308682e-01; r[490][Y]= 9.290232e-01; r[490][Z]= 5.782360e-02; q[490]= 9.993093e-06;
    r[491][X]= 4.528262e-01; r[491][Y]= 2.204809e-01; r[491][Z]= 7.665425e-01; q[491]= 9.993093e-06;
    r[492][X]= 7.351880e-02; r[492][Y]= 5.337377e-01; r[492][Z]= 5.381479e-01; q[492]= 9.993093e-06;
    r[493][X]= 4.829276e-01; r[493][Y]= 2.663587e-01; r[493][Z]= 8.132023e-01; q[493]= 9.993093e-06;
    r[494][X]= 3.521615e-01; r[494][Y]= 5.654891e-01; r[494][Z]= 3.139011e-01; q[494]= 9.993093e-06;
    r[495][X]= 4.975477e-01; r[495][Y]= 7.544468e-01; r[495][Z]= 5.271600e-01; q[495]= 9.993093e-06;
    r[496][X]= 4.658774e-01; r[496][Y]= 2.290643e-01; r[496][Z]= 7.838834e-01; q[496]= 9.993093e-06;
    r[497][X]= 2.807409e-01; r[497][Y]= 6.284957e-01; r[497][Z]= 2.153700e-02; q[497]= 9.993093e-06;
    r[498][X]= 4.682633e-01; r[498][Y]= 8.296577e-01; r[498][Z]= 7.368002e-01; q[498]= 9.993093e-06;
    r[499][X]= 4.733331e-01; r[499][Y]= 3.019186e-01; r[499][Z]= 8.890947e-01; q[499]= 9.993093e-06;
    r[500][X]= 7.012500e-02; r[500][Y]= 4.189419e-01; r[500][Z]= 8.550628e-01; q[500]= 9.993093e-06;
    r[501][X]= 1.441829e-01; r[501][Y]= 7.913139e-01; r[501][Z]= 6.132985e-01; q[501]= 9.993093e-06;
    r[502][X]= 1.870435e-01; r[502][Y]= 5.026524e-01; r[502][Z]= 5.999738e-01; q[502]= 9.993093e-06;
    r[503][X]= 2.450834e-01; r[503][Y]= 6.876646e-01; r[503][Z]= 8.384747e-01; q[503]= 9.993093e-06;
    r[504][X]= 5.283172e-01; r[504][Y]= 8.463487e-01; r[504][Z]= 1.576466e-01; q[504]= 9.993093e-06;
    r[505][X]= 9.691857e-01; r[505][Y]= 2.890469e-01; r[505][Z]= 7.954612e-01; q[505]= 9.993093e-06;
    r[506][X]= 3.090546e-01; r[506][Y]= 9.971512e-01; r[506][Z]= 9.429721e-01; q[506]= 9.993093e-06;
    r[507][X]= 6.280935e-01; r[507][Y]= 9.503080e-02; r[507][Z]= 3.811045e-01; q[507]= 9.993093e-06;
    r[508][X]= 6.811633e-01; r[508][Y]= 3.817675e-01; r[508][Z]= 6.175555e-01; q[508]= 9.993093e-06;
    r[509][X]= 6.721672e-01; r[509][Y]= 1.791827e-01; r[509][Z]= 2.731031e-01; q[509]= 9.993093e-06;
    r[510][X]= 9.604940e-01; r[510][Y]= 5.306052e-01; r[510][Z]= 9.522110e-01; q[510]= 9.993093e-06;
    r[511][X]= 2.019551e-01; r[511][Y]= 3.950060e-02; r[511][Z]= 3.419927e-01; q[511]= 9.993093e-06;
}