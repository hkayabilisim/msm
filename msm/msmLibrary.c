#include "msmLibrary.h"

float myerfc(float x)
{
    float t,z,ans;
    z=fabs(x);
    t=1.0/(1.0+0.5*z);
    ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
    t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
    t*(-0.82215223+t*0.17087277)))))))));
    return x >= 0.0 ? ans : 2.0-ans;
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
    while ( fabs( (myerfc(x*a)/a - e * h)/(-2*exp(-x*x*a*a)/sqrt(PI)))/fabs(x) > e  && iter < maxiter)
    {
        x = x - (myerfc(x*a)/a - e * h)/(-2*exp(-x*x*a*a)/sqrt(PI));
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
    int Fx = 0 , Fy = 0, Fz = 0; // Number of factors of 2 in Mx, My and Mz
    M = Mx*My*Mz; // # grid points in the first level
    
    while ((Mx % 2) == 0) { Fx++; Mx /= 2; }
    while ((My % 2) == 0) { Fy++; My /= 2; }
    while ((Mz % 2) == 0) { Fz++; Mz /= 2; }
    Lallowed = minOf3Integers(Fx,Fy,Fz) + 1;
    Lwanted = ceil(1.0 + (1./3)*log2(M/sqrt(N)));
    
    if (Lwanted > Lallowed) {
        fprintf(stderr,"There is not enough factor of 2s. Dropping L from %d to %d\n",Lwanted,Lallowed);
        L = Lallowed;
    } else
        L = Lwanted;
    if (L != 1) {
        // TODO: fprintf(stderr,"L = 1 manually !!!\n");
        L = 1;
    }
    return L;
}
int chooseM(double htilde,double Ax)
{
    int M = 0 ;
    // h/l * 5/6 <= M <= h/l 5/4
    double lalpha =  Ax;
    
    double c5over8 = (5.0/8.0) * lalpha / htilde;
    double c5over6 = (5.0/6.0) * lalpha / htilde;
    double c5over4 = (5.0/4.0) * lalpha / htilde;
    
    int iupper = ceil(log2(c5over4));
    
    // TODO: I can also loop starting from top to bottom.
    int kestimate = 0;
    
    for (int i = 0 ; i <= iupper; i++) {
        if (c5over8 < pow(2,i)  && pow(2,i) <= c5over4) {
            M = pow(2,i);
            kestimate = i;
            break;
        }
    }
    
    if (pow(2,kestimate) < c5over6) {
        int iupper = ceil(log2(c5over6));
        M = 0;
        for (int i = 0 ; i <= iupper; i++) {
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
    while ( fabs( (myerfc(PI*x/beta)*2*beta/sqrt(PI)-e*h)  /(-4*exp(-x*x*PI*PI/(beta*beta))))/fabs(x) > e  && iter < maxiter)
    {
        x = x - (myerfc(PI*x/beta)*2*beta/sqrt(PI)-e*h) /(-4*exp(-x*x*PI*PI/(beta*beta)));
        iter++;
    }
    if (iter >= maxiter) {
        fprintf(stderr,"kmax didn't convergenge\n");
    }
    return x;
}
double short_range_periodic(int i,int j, double *r,double *q,int s,double a,double *fshort,double Ax,double Ay,double Az)
{
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
    for (int px = pxmin ; px <= pxmax ; px++)
    {
        for (int py = pymin ; py <= pymax ; py++)
        {
            for (int pz = pzmin ; pz <= pzmax ; pz++)
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
                    //double gvalue = gammarho(rlen/a, s);
                    //double gvalued = gammarhod(rlen/a, s);
                    double g =  1.0/rlen - (1./a) * TAUP4(rlen/a);
                    double forcex = 0; //q[i]*q[j]*(1./r2 + 1/(a*a)*gvalued)*(rix-rjx)/rlen;
                    double forcey = 0; //q[i]*q[j]*(1./r2 + 1/(a*a)*gvalued)*(riy-rjy)/rlen;
                    double forcez = 0; //q[i]*q[j]*(1./r2 + 1/(a*a)*gvalued)*(riz-rjz)/rlen;
                    ushort +=  g;
                    fshort[MAXDIM*i+DIMX] += -forcex;
                    fshort[MAXDIM*i+DIMY] += -forcey;
                    fshort[MAXDIM*i+DIMZ] += -forcez;
                    fshort[MAXDIM*j+DIMX] +=  forcex;
                    fshort[MAXDIM*j+DIMY] +=  forcey;
                    fshort[MAXDIM*j+DIMZ] +=  forcez;
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
