#include "msmLibrary.h"

double mylog2(double x) {
    return log(x) / log(2.0);
}

int data_read(char *filename, double **q, double **r, int *N, double *Lx, double *Ly, double *Lz)
{
    int nn, i;
    double *rr, *qq;
    FILE *fp;
    
    fp = fopen(filename, "r");
    fscanf(fp, "%lf %lf %lf\n", Lx, Ly, Lz);
    fscanf(fp, "%d\n", &nn);
    
    *N = nn;
    *q = (double *)calloc(nn, sizeof(double));
    *r = (double *)calloc(nn * 3, sizeof(double));
    
    rr = *r;
    qq = *q;
    for (i = 0; i < nn; i++)
    {
        fscanf(fp, "%lf %lf %lf %lf", &(rr[i * 3]), &(rr[i * 3 + 1]), &(rr[i * 3 + 2]), &(qq[i]));
    }
    fclose(fp);
    return 0;
}

double choosebetaNew(double e,double a,double h){
    double x = e;
    int iter = 1;
    int maxiter = 200;
    while ( fabs( (erfc(x*a)/a - e * h)/(-2*exp(-x*x*a*a)/sqrt(PI)))/fabs(x) > e  && iter < maxiter)
    {
        x = x - (erfc(x*a)/a - e * h)/(-2*exp(-x*x*a*a)/sqrt(PI));
        iter++;
    }
    if (iter >= maxiter) {
        fprintf(stderr,"beta didn't convergenge\n");
    }
    return x;
}

int chooseL(int Mx,int My,int Mz,int N)
{
    int L, Lwanted, Lallowed, M;
    int Fx = 0 , Fy = 0, Fz = 0; /* Number of factors of 2 in Mx, My and Mz */
    M = Mx*My*Mz; /* # grid points in the first level */
    
    while ((Mx % 2) == 0) { Fx++; Mx /= 2; }
    while ((My % 2) == 0) { Fy++; My /= 2; }
    while ((Mz % 2) == 0) { Fz++; Mz /= 2; }
    Lallowed = minOf3Integers(Fx,Fy,Fz) + 1;
    Lwanted = ceil(1.0 + (1./3)*mylog2(M/sqrt(N)));
    
    if (Lwanted > Lallowed) {
        //fprintf(stderr,"There is not enough factor of 2s. Dropping L from %d to %d\n",Lwanted,Lallowed);
        L = Lallowed;
    } else
        L = Lwanted;
    /* TODO: fprintf(stderr,"L = 1 manually !!!\n"); */
    if (L != 1) {
        L = 1;
    }  
    return L;
}
int chooseM(double htilde,double Ax)
{
    int i,M = 0 ;
    /* h/l * 5/6 <= M <= h/l 5/4 */
    double lalpha =  Ax;
    
    double c5over8 = (5.0/8.0) * lalpha / htilde;
    double c5over6 = (5.0/6.0) * lalpha / htilde;
    double c5over4 = (5.0/4.0) * lalpha / htilde;
    
    int iupper = ceil(mylog2(c5over4));
    
    /* TODO: I can also loop starting from top to bottom. */
    int kestimate = 0;
    
    for ( i = 0 ; i <= iupper; i++) {
        if (c5over8 < pow(2,i)  && pow(2,i) <= c5over4) {
            M = pow(2,i);
            kestimate = i;
            break;
        }
    }
    
    if (pow(2,kestimate) < c5over6) {
        int iupper = ceil(mylog2(c5over6));
        M = 0;
        for ( i = 0 ; i <= iupper; i++) {
            if (c5over8 < pow(2,i)  && pow(2,i) <= c5over6) {
                M = 3*pow(2,i-1);
                break;
            }
        }
    }
    
    if (M == 0) {
        fprintf(stderr,"M is found to be zero\n");
    }
    
    return M;
}

double choosekmax(double e,double beta,double h) {
    double x=0.0;
    
    int iter = 1;
    int maxiter = 200;
    while ( fabs( (erfc(PI*x/beta)*2*beta/sqrt(PI)-e*h)  /(-4*exp(-x*x*PI*PI/(beta*beta))))/fabs(x) > e  && iter < maxiter)
    {
        x = x - (erfc(PI*x/beta)*2*beta/sqrt(PI)-e*h) /(-4*exp(-x*x*PI*PI/(beta*beta)));
        iter++;
    }
    if (iter >= maxiter) {
        fprintf(stderr,"kmax didn't convergenge\n");
    }
    return x;
}
double short_range_periodic(int i,int j, double *r,double *q,int s,double a,double *fshort,double Ax,double Ay,double Az)
{
    int px,py,pz;
    double ushort = 0.0;
    double a2 ;
    double rix = r[MAXDIM*i + DIMX];
    double riy = r[MAXDIM*i + DIMY];
    double riz = r[MAXDIM*i + DIMZ];
    
    double rjx  = r[MAXDIM*j + DIMX];
    double rjy  = r[MAXDIM*j + DIMY];
    double rjz  = r[MAXDIM*j + DIMZ];
    
    double rijx = rix - rjx;
    double rijy = riy - rjy;
    double rijz = riz - rjz;
    int pxmin = ( rijx - a ) / Ax ;
    int pxmax = ( rijx + a ) / Ax ;
    int pymin = ( rijy - a ) / Ay ;
    int pymax = ( rijy + a ) / Ay ;
    int pzmin = ( rijz - a ) / Az ;
    int pzmax = ( rijz + a ) / Az ;
    
    a2 = a * a;
    for ( px = pxmin ; px <= pxmax ; px++)
    {
        for ( py = pymin ; py <= pymax ; py++)
        {
            for ( pz = pzmin ; pz <= pzmax ; pz++)
            {
                if (px == 0 && py == 0 && pz == 0 && i==j)
                    continue;
                double Aprx = Ax * px - rijx;
                double Apry = Ay * py - rijy;
                double Aprz = Az * pz - rijz;
                double r2   = Aprx * Aprx + Apry * Apry + Aprz * Aprz ;
                if (r2 < a2)
                {
                    double rlen = sqrt(r2);
                    /* double gvalue = gammarho(rlen/a, s);
                     double gvalued = gammarhod(rlen/a, s); */
                    double g =  1.0/rlen - (1./a) * TAUP4(rlen/a);
                    ushort +=  g;
	                    double forcex = 0.5*q[i]*(1./r2 + 1./(a*a)*TAUP4D(rlen/a))*Aprx/rlen;
	                    double forcey = 0.5*q[i]*(1./r2 + 1./(a*a)*TAUP4D(rlen/a))*Apry/rlen;
	                    double forcez = 0.5*q[i]*(1./r2 + 1./(a*a)*TAUP4D(rlen/a))*Aprz/rlen;
	                    fshort[MAXDIM*i+DIMX] +=  forcex;
	                    fshort[MAXDIM*i+DIMY] +=  forcey;
	                    fshort[MAXDIM*i+DIMZ] +=  forcez;
	                    fshort[MAXDIM*j+DIMX] += -forcex;
	                    fshort[MAXDIM*j+DIMY] += -forcey;
	                    fshort[MAXDIM*j+DIMZ] += -forcez;
                }
            }
        }
    }
    return ushort;
}
void msm4g_linkedlist_add(LinkedList *list,int data)
{
    LinkedListElement *item;
    
    item = (LinkedListElement *)malloc(sizeof(LinkedListElement));
    item->data = data;
    item->next = NULL;
    item->prev = NULL;
    
    /* the list is empty */
    if ( list->tail == NULL)
    {
        list->head = item;
        list->tail = item;
    } else
    {
        item->prev = list->tail;
        list->tail->next = item;
        list->tail = item;
    }
}

void msm4g_linkedlist_addvalue(LinkedList *list,int data,double value)
{
    LinkedListElement *item;
    
    item = (LinkedListElement *)malloc(sizeof(LinkedListElement));
    item->data = data;
    item->value = value;
    item->next = NULL;
    item->prev = NULL;
    
    /* the list is empty */
    if ( list->tail == NULL)
    {
        list->head = item;
        list->tail = item;
    } else
    {
        item->prev = list->tail;
        list->tail->next = item;
        list->tail = item;
    }
}

void msm4g_linkedlist_destroy(LinkedList *list)
{
    LinkedListElement *current, *next;
    current = list->head;
    while (current != NULL)
    {
        next=current->next;
        free(current);
        current=next;
    }
    free(list);
    list = NULL;
}

double msm4g_linkedlist_getvalue(LinkedList *list,int index)
{
    int order = 0;
    LinkedListElement *curr;
    curr = list->head;
    while (curr != NULL)
    {
        if (order == index)
            return curr->value;
        curr=curr->next;
        order++;
    }
    
    return -1;
}

LinkedList *msm4g_linkedlist_new()
{
    LinkedList *newlist ;
    newlist = (LinkedList *)malloc(sizeof(LinkedList));
    newlist->head = NULL;
    newlist->tail = NULL;
    return newlist;
}

int msm4g_linkedlist_search(LinkedList *list,int data)
{
    int index = 0;
    LinkedListElement *curr;
    curr = list->head;
    while (curr != NULL)
    {
        if (curr->data == data)
            return index;
        curr=curr->next;
        index++;
    }
    
    return -1;
}

int minOf3Integers(int a, int b, int c)
{
    int min = a <= b ? a : b;
    min  = min <= c ? min : c;
    return min;
}

void print_stats(msm_stats stats)
{
    int i,l;
    printf("%30s %d\n","N",stats.N);
    printf("%30s %d %d %d\n","Ax Ay Az",stats.Ax,stats.Ay,stats.Az);
    printf("%30s %d\n","L",stats.L);
    printf("%30s %d\n","mu",stats.mu);
    printf("%30s %25.16e\n","h (estimated)",stats.h);
    printf("%30s %f %f %f\n","h (calculated)",stats.hx,stats.hy,stats.hz);
    printf("%30s %25.16e\n","a",stats.a);
    printf("%30s %25.16e\n","aL",stats.aL);
    printf("%30s %25.16e\n","alpha",stats.alpha);
    printf("%30s %25.16e\n","beta",stats.beta);
    printf("%30s %d\n","kmax",stats.kmax);
    printf("%30s %d\n","p",stats.p);
    printf("%30s %d\n","s",stats.s);
    printf("%30s %f\n","csr",stats.csr);
    printf("%30s %f\n","chisum",stats.chisum);
    printf("%30s %d %d %d\n","chinx chiny chinz",stats.chinx,stats.chiny,stats.chinz);
    printf("%30s %d\n","chisize",stats.chisize);
    printf("%30s ","chi");
    for (i = 0 ; i < stats.chisize ; i++)
        printf("%6.2e ",stats.chi[i]);
    printf("\n");
    printf("%30s %f %f %f\n","rminx rminy rminz",stats.rminx,stats.rminy,stats.rminz);
    printf("%30s %f %f %f\n","rmaxx rmaxy rmaxz",stats.rmaxx,stats.rmaxy,stats.rmaxz);
    printf("%30s %d %d %d\n","Mx My Mz",stats.Mx,stats.My,stats.My);
    printf("%30s %d %d %d\n","MLx MLy MLz",stats.MLx,stats.MLy,stats.MLz);
    printf("%30s %d %d %d\n","ML2x ML2y ML2z",stats.ML2x,stats.ML2y,stats.ML2z);
    for (l = 0 ; l <= stats.L; l++)
        printf("%30s (l=%d) %d %d %d\n","Mlx Mly Mlz",l,stats.Mlx[l],stats.Mly[l],stats.Mlz[l]);
    printf("%30s %d %d %d\n","clenx cleny clenz",stats.clenx,stats.cleny,stats.clenz);
    printf("%30s ","cvecx");
    for (i = 0 ; i < stats.clenx ; i++)
        printf("%6.2e ",stats.cvecx[i]);
    printf("\n");
    printf("%30s ","cvecy");
    for (i = 0 ; i < stats.cleny ; i++)
        printf("%6.2e ",stats.cvecy[i]);
    printf("\n");
    printf("%30s ","cvecz");
    for (i = 0 ; i < stats.clenz ; i++)
        printf("%6.2e ",stats.cvecz[i]);
    printf("\n");
    printf("%30s %d x %d x %d = %d\n","sLnx x sLny x sLnz = sLsize",stats.sLnx,stats.sLny,stats.sLnz,stats.sLsize);
    for (l = 0 ; l < stats.L; l++)
        printf("%30s (l=%d) %d %d %d\n","slsizex slsizey slsizez",l,stats.slsizex[l],stats.slsizey[l],stats.slsizez[l]);
    printf("%30s %f\n","time total",stats.cputime_total);
    printf("%30s %f\n","time ulong",stats.cputime_potential_long);
    printf("%30s %f\n","time stencilL",stats.cputime_stencilL);
    printf("%30s %f\n","time stencill",stats.cputime_stencill);
    printf("%30s %f\n","time top_level",stats.cputime_top_level);
    printf("%30s %f\n","time restriction",stats.cputime_restriction);
    printf("%30s %f\n","time short_range",stats.cputime_short_range);
    printf("%30s %f\n","time self",stats.cputime_self);
    printf("%30s %f\n","time prolongation",stats.cputime_prolongation);
    printf("%30s %f\n","time anterpolation",stats.cputime_anterpolation);
    printf("%30s %f\n","time interpolation",stats.cputime_interpolation);
    printf("%30s %f\n","time init",stats.cputime_initialization);
    printf("%30s %f\n","time grid to grid",stats.cputime_grid_to_grid_mapping);
    
    printf("%30s %f\n","netcharge",stats.netcharge);
    for (l = 0 ; l < stats.L; l++)
        printf("%30s (l=%d) %f\n","netgridcharge",l,stats.netgridcharge[l]);
    
    printf("%-20s : %25.16e\n","ushort_real(MSM)",stats.potential_short_real);
    printf("%-20s : %25.16e\n","ushort_self(MSM)",stats.potential_short_self);
    printf("%-20s : %25.16e\n","ushort_corr(MSM)",stats.potential_short_csr);
    printf("%-20s : %25.16e\n","ushort_tota(MSM)",stats.potential_short_total);
    
    printf("%-20s : %25.16e\n","ulong_real(MSM)",stats.potential_long_real);
    printf("%-20s : %25.16e\n","ulong_self(MSM)",stats.potential_long_self);
    printf("%-20s : %25.16e\n","ulong_four(MSM)",stats.potential_long_fourier);
    printf("%-20s : %25.16e\n","ulong_tota(MSM)",stats.potential_long_total);
    printf("%-20s : %25.16e\n","utotal",stats.potential_total);

    printf("%-20s : %25.16e\n","(A) ushort_real",stats.potential_short_real);
    printf("%-20s : %25.16e\n","(B) ushort_self",stats.potential_short_self);
    printf("%-20s : %25.16e\n","(C) ushort_csr" ,stats.potential_short_csr);
    printf("%-20s : %25.16e\n","(D) ulong_real" ,stats.potential_long_real);
    printf("%-20s : %25.16e\n","(E) ulong_four" ,stats.potential_long_fourier);
    printf("%-20s : %25.16e\n","(F) ulong_self" ,stats.potential_long_self);
    printf("%-20s : %25.16e\n","    utotal"     ,stats.potential_total);
    
    
    
}
