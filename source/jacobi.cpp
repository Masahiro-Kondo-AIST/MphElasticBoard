#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

#define DIM 3

static int converge(double A[3][3]){
    return A[0][1]*A[0][1]+A[1][2]*A[1][2]+A[2][0]*A[2][0] <= 1.0e-16*(A[0][0]*A[0][0]+A[1][1]*A[1][1]+A[2][2]*A[2][2]);
}

static double sign(double d){
    if(d<0)return -1;
    else return 1;
}

int jacobi(double R[3][3], double A[3][3])
{
    int iter=0;
    double theta, t, tau, c, s;
    R[0][0]=1;R[0][1]=0;R[0][2]=0;
    R[1][0]=0;R[1][1]=1;R[1][2]=0;
    R[2][0]=0;R[2][1]=0;R[2][2]=1;
    
    while(1){
        for(int p=0;p<DIM;++p){
            for(int q=p+1;q<DIM;++q){
                if(converge(A))goto converged;
                if(fabs(A[p][q])<=1.0e-100)continue;
                theta = (A[q][q]-A[p][p])/2.0/A[p][q];
                t = sign(theta)/(fabs(theta)+sqrt(theta*theta+1));
                c=1/sqrt(t*t+1);
                s=t*c;
                tau=s/(1+c);
                double Apq=A[p][q],App=A[p][p],Aqq=A[q][q];
                A[p][q]=0;
                A[q][p]=0;
                A[p][p]=App-t*Apq;
                A[q][q]=Aqq+t*Apq;
                for(int r=0;r<DIM;++r){
                    double Rrp=R[r][p],Rrq=R[r][q];
                    R[r][p]=Rrp-s*(Rrq+tau*Rrp);
                    R[r][q]=Rrq+s*(Rrp-tau*Rrq);
                    if(r==p || r==q)continue;
                    double Arp=A[r][p],Arq=A[r][q];
                    A[p][r]=A[r][p]=Arp-s*(Arq+tau*Arp);
                    A[q][r]=A[r][q]=Arq+s*(Arp-tau*Arq);
                }
                iter++;
            }
        }
    }
 converged:
    // fprintf(stderr, "iter: %d\n", iter);
    return 0;
}
