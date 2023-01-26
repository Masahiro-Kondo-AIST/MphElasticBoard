//================================================================================================//
//------------------------------------------------------------------------------------------------//
//    MPH-I : Moving Particle Hydrodynamics for Elastic Board                                     //
//------------------------------------------------------------------------------------------------//
//    Developed by    : Masahiro Kondo                                                            //
//    Distributed in  : 2023                                                                      //
//    Lisence         : GPLv3                                                                     //
//    For instruction : see README                                                                //
//    For theory      : see the following references                                              //
//     [1] JSCES Paper No.20100016,   https://doi.org/10.11421/jsces.2010.20100016                //
//    Copyright (c) 2010   Masahiro Kondo                                                         //
//================================================================================================//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <assert.h>

#include "errorfunc.h"
#include "log.h"
#include "jacobi.h"

const double DOUBLE_ZERO[32]={0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0};

using namespace std;

#define DIM 3

#define FREE_PARTICLE  0
#define FIXED_PARTICLE 1


#define MAX_STORED_COUNT 256
typedef struct{
    int count;
    int stored[MAX_STORED_COUNT];
}t_plist;
const int INT_MINUS[MAX_STORED_COUNT]={ -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                        -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1,
                                       -1,-1,-1,-1,-1,-1,-1,-1, -1,-1,-1,-1,-1,-1,-1,-1};

#define  DEFAULT_LOG  "sample.log"
#define  DEFAULT_DATA "sample.data"
#define  DEFAULT_GRID "sample.grid"
#define  DEFAULT_PROF "sample%03d.prof"
#define  DEFAULT_ENE  "sample.ene"

// Calculation and Output
static double DomainMin[DIM];
static double DomainMax[DIM];
static double OutputInterval=0.0;
static double OutputNext=0.0;
static double EndTime=0.0;
static double Time=0.0;
static double Dt=1.0e100;
static double Radius=0.0;


// Phisics
static double ParticleSpacing;
static double Density;
static double Young;
static double Poiss;
static double Thickness;
static double CoefficientTensor[DIM][DIM][DIM][DIM];  //係数行列(coeffieint tensor)

// Particle
static int ParticleCount;
static int *Property;              // 粒子種類(particle type)
static double (*InitialPosition)[DIM];           // 初期座標(initial corrdinate)
static double (*Position)[DIM];           // 座標(coordinate)
static double (*Velocity)[DIM];           // 運動量(momentum)
static double (*Acceleration)[DIM];       // 加速度(acceleration)
static double (*DeformationGradient)[DIM][DIM];      // 変形勾配テンソル(deformation gradient tensor)
static double (*CurvatureTensor)[DIM][DIM][DIM]; // 曲率テンソル(curvature tensor)
static double (*NormalizeTensor)[DIM][DIM];     // 規格化テンソル(normalizer)
static t_plist (*InitialNeighbor); // 近傍粒子(終了記号:"-1")(neighbor particles (end with "-1")


// BackGroundCells
#define CellId(iCX,iCY,iCZ)  ((iCX)*CellCount[1]*CellCount[2]+(iCY)*CellCount[2]+(iCZ))
static int CellCount[DIM];
static double MaxRadius = 0.0;
static double CellWidth = 0.0;
static int *CellParticleCount;  // セルに入る粒子数(number of particles in the cell)
static int *CellParticleBegin;  // セルの始まり(beginning of particles int the cell)
static int *CellParticle;       // セルに入る粒子番号の配列(array of particle id in the cells)

// Energy
static double Kinetic;
static double BendingPotential;

static void readDataFile(char *filename);
static void readGridFile(char *filename);
static void writeProfFile(char *filename);
static void initializeDomain( void );
static void calculateNeighbor( t_plist *neighbor, const double (*position)[DIM] );
static void calculateCoefficientTensor();
static void calculateNormalizer();
static void resetAcceleration();
static void calculateDeformationGradient();
static void calculateCurvatureTensor();
static void calculateMotionEquation();
static void calculateExternalForce();
static void calculateConvection();

static void writeEnergyFile(char *filename);

int main(int argc, char *argv[])
{
    char logfilename[1024]  = DEFAULT_LOG;
    char datafilename[1024] = DEFAULT_DATA;
    char gridfilename[1024] = DEFAULT_GRID;
    char proffilename[1024] = DEFAULT_PROF;
    char energyfilename[1024]=DEFAULT_ENE;
    {
        if(argc>1)strcpy(datafilename,argv[1]);
        if(argc>2)strcpy(gridfilename,argv[2]);
        if(argc>3)strcpy(proffilename,argv[3]);
        if(argc>4)strcpy(energyfilename,argv[4]);
        if(argc>5)strcpy( logfilename,argv[5]);
    }
    log_open(logfilename);
    {
        time_t t=time(NULL);
        log_printf("start reading files at %s\n",ctime(&t));
    }
    readDataFile(datafilename);
    readGridFile(gridfilename);
    {
        time_t t=time(NULL);
        log_printf("start initialization at %s\n",ctime(&t));
    }
    initializeDomain();
    calculateNeighbor( InitialNeighbor, InitialPosition );
    calculateCoefficientTensor();
    calculateNormalizer();
    {
        time_t t=time(NULL);
        log_printf("start main roop at %s\n",ctime(&t));
    }
    int iStep=(int)(Time/Dt);
    while(Time < EndTime + 1.0e-5*Dt){
        resetAcceleration();
        calculateDeformationGradient();
        calculateCurvatureTensor();
        calculateMotionEquation();
        calculateExternalForce();
        if( Time + 1.0e-5*Dt >= OutputNext ){
            char filename[256];
            sprintf(filename,proffilename,iStep);
            writeProfFile(filename);
            log_printf("@ Prof Output Time : %e\n", Time );
            OutputNext += OutputInterval;
        }
        writeEnergyFile( energyfilename );
        calculateConvection();

        Time += Dt;
        iStep++;
    }
    {
        time_t t=time(NULL);
        log_printf("end main roop at %s\n",ctime(&t));
    }
    
    return 0;
    
}

static void readDataFile(char *filename)
{
    FILE * fp;
    char buf[1024];
    // char command[1024];
    // char sval[1024];
    // int ival;
    // double dval;
    const int reading_global=0;
    // jconst int reading_particle_property=1;
    int mode=reading_global;
    // int iProperty;
    

    fp=fopen(filename,"r");
    mode=reading_global;
    while(fp!=NULL && !feof(fp) && !ferror(fp)){
        fgets(buf,sizeof(buf),fp);
        if(buf[0]=='#'){}
        else if(sscanf(buf," Dt %lf",&Dt)==1){mode=reading_global;}
        else if(sscanf(buf," OutputInterval %lf",&OutputInterval)==1){mode=reading_global;}
        else if(sscanf(buf," EndTime %lf",&EndTime)==1){mode=reading_global;}
        else if(sscanf(buf," Radius %lf",&Radius)==1){mode=reading_global;}
        else if(sscanf(buf," Density %lf", &Density)==1){mode=reading_global;}
        else if(sscanf(buf," Young %lf", &Young)==1){mode=reading_global;}
        else if(sscanf(buf," Poiss %lf", &Poiss)==1){mode=reading_global;}
        else if(sscanf(buf," Thickness %lf", &Thickness)==1){mode=reading_global;}
        else{
            log_printf("Invalid line in data file \"%s\"\n", buf);
        }
    }
    fclose(fp);
    return;
}

static void readGridFile(char *filename)
{
    FILE *fp=fopen(filename,"r");
    char buf[1024];   

    fgets(buf,sizeof(buf),fp);
    sscanf(buf,"%lf",&Time);
    fgets(buf,sizeof(buf),fp);
    sscanf(buf,"%d  %lf  %lf %lf %lf  %lf %lf %lf",
           &ParticleCount,
           &ParticleSpacing,
           &DomainMin[0], &DomainMax[0],
           &DomainMin[1], &DomainMax[1],
           &DomainMin[2], &DomainMax[2]);

    Property = (int *)malloc(ParticleCount*sizeof(int));
    InitialPosition = (double (*)[DIM])malloc(ParticleCount*sizeof(double [DIM]));
    Position = (double (*)[DIM])malloc(ParticleCount*sizeof(double [DIM]));
    Velocity = (double (*)[DIM])malloc(ParticleCount*sizeof(double [DIM]));
    Acceleration = (double (*)[DIM])malloc(ParticleCount*sizeof(double [DIM]));
    DeformationGradient = (double (*)[DIM][DIM])malloc(ParticleCount*sizeof(double [DIM][DIM]));
    CurvatureTensor = (double (*)[DIM][DIM][DIM])malloc(ParticleCount*sizeof(double [DIM][DIM][DIM]));
    NormalizeTensor= (double (*)[DIM][DIM])malloc(ParticleCount*sizeof(double [DIM][DIM]));
    InitialNeighbor = (t_plist (*))malloc(ParticleCount*sizeof(t_plist));

    
    double (*a)[DIM] = InitialPosition;
    double (*q)[DIM] = Position;
    double (*v)[DIM] = Velocity;
    
    for(int iP=0;iP<ParticleCount;++iP){
        fgets(buf,sizeof(buf),fp);
        sscanf(buf,"%d  %lf %lf %lf  %lf %lf %lf  %lf %lf %lf",
               &Property[iP],
               &q[iP][0],&q[iP][1],&q[iP][2],
               &v[iP][0],&v[iP][1],&v[iP][2],
               &a[iP][0],&a[iP][1],&a[iP][2]
               );
    }
    fclose(fp);
    return;
}

static void writeProfFile(char *filename)
{
    FILE *fp=fopen(filename,"w");

    fprintf(fp,"%e\n",Time);
    fprintf(fp,"%d %e %e %e %e %e %e %e\n",
            ParticleCount,
            ParticleSpacing,
            DomainMin[0], DomainMax[0],
            DomainMin[1], DomainMax[1],
            DomainMin[2], DomainMax[2]);

    const double (*a)[DIM] = InitialPosition;
    const double (*q)[DIM] = Position;
    const double (*v)[DIM] = Velocity;
    
    for(int iP=0;iP<ParticleCount;++iP){
            fprintf(fp,"%d %e %e %e %e %e %e %\n",
                    Property[iP],
                    q[iP][0], q[iP][1], q[iP][2],
                    v[iP][0], v[iP][1], v[iP][2],
                    a[iP][0], a[iP][1], a[iP][2]);
    }
    fflush(fp);
    fclose(fp);
}

static void initializeDomain( void )
{
    CellWidth = ParticleSpacing;
    MaxRadius = Radius;
    
    int range = (int)(MaxRadius/CellWidth)+1;
    CellCount[0] = (int)((DomainMax[0] - DomainMin[0])/CellWidth)+1 +2*range;
    CellCount[1] = (int)((DomainMax[1] - DomainMin[1])/CellWidth)+1 +2*range;
    CellCount[2] = (int)((DomainMax[2] - DomainMin[2])/CellWidth)+1 +2*range;

    CellParticleCount = (int *)malloc( CellCount[0]*CellCount[1]*CellCount[2]*sizeof(int) );
    CellParticleBegin = (int *)malloc( CellCount[0]*CellCount[1]*CellCount[2]*sizeof(int) );

    CellParticle = (int *)malloc( ParticleCount * sizeof(int) );
}

static void calculateNeighbor( t_plist (*neighbor), const double (*position)[DIM] )
{
    int range = (int)(MaxRadius/CellWidth) + 1;
    
    const double (*q)[DIM] = position;
    for(int iC=0;iC<CellCount[0]*CellCount[1]*CellCount[2];++iC){
        CellParticleCount[iC]=0;
    }
    for(int iP=0; iP<ParticleCount; ++iP){
        const int iCX=(int)((q[iP][0]-DomainMin[0])/CellWidth)+range;
        const int iCY=(int)((q[iP][1]-DomainMin[1])/CellWidth)+range;
        const int iCZ=(int)((q[iP][2]-DomainMin[2])/CellWidth)+range;
        const int iC=CellId(iCX,iCY,iCZ);
        CellParticleCount[iC]++;
    }
    int iCellParticle = 0;
    for(int iC=0;iC<CellCount[0]*CellCount[1]*CellCount[2];++iC){
        CellParticleBegin[iC] = iCellParticle;
        iCellParticle+=CellParticleCount[iC];
    }
    for(int iC=0;iC<CellCount[0]*CellCount[1]*CellCount[2];++iC){
        CellParticleCount[iC]=0;
    }
    for(int iP=0; iP<ParticleCount; ++iP){
        const int iCX=(int)((q[iP][0]-DomainMin[0])/CellWidth)+range;
        const int iCY=(int)((q[iP][1]-DomainMin[1])/CellWidth)+range;
        const int iCZ=(int)((q[iP][2]-DomainMin[2])/CellWidth)+range;
        const int iC=CellId(iCX,iCY,iCZ);
        CellParticle[ CellParticleBegin[iC] + CellParticleCount[iC] ] = iP;
        CellParticleCount[iC]++;
    }
    
    // calculate neighbor
    for(int iP=0;iP<ParticleCount;++iP){
        neighbor[iP].count=0;
        memcpy(neighbor[iP].stored,INT_MINUS,sizeof(int)*MAX_STORED_COUNT);
    }
    for(int iP=0;iP<ParticleCount;++iP){
        const int iCX=(int)((q[iP][0]-DomainMin[0])/CellWidth)+range;
        const int iCY=(int)((q[iP][1]-DomainMin[1])/CellWidth)+range;
        const int iCZ=(int)((q[iP][2]-DomainMin[2])/CellWidth)+range;
        
        for(int jCX=iCX-range;jCX<=iCX+range;++jCX){
            for(int jCY=iCY-range;jCY<=iCY+range;++jCY){
                for(int jCZ=iCZ-range;jCZ<=iCZ+range;++jCZ){
                    const int jC=CellId(jCX,jCY,jCZ);
                    for(int jCP=CellParticleBegin[jC];jCP<CellParticleBegin[jC]+CellParticleCount[jC];++jCP){
                        int jP=CellParticle[jCP];
                        const double qij[DIM] = {q[jP][0]-q[iP][0],q[jP][1]-q[iP][1],q[jP][2]-q[iP][2]};;
                        const double qij2= qij[0]*qij[0]+qij[1]*qij[1]+qij[2]*qij[2];
                        if(qij2 <= Radius*Radius){
                            if(neighbor[iP].count>=MAX_STORED_COUNT){
                                log_printf("Too many neighbors\nNeighbors must be less than %d", MAX_STORED_COUNT);
                                exit(0);
                            }
                            neighbor[iP].stored[neighbor[iP].count] = jP;
                            neighbor[iP].count++;
                        }
                    }
                }
            }
        }
    }
}

static void calculateCoefficientTensor()
{
    double (*K)[DIM][DIM][DIM] = CoefficientTensor;
    double E = Young;
    double v = Poiss;
    double h = Thickness;
    double D = E*h*h*h/(12.0*(1-v*v));

    for(int sD=0;sD<DIM;++sD){
        for(int tD=0;tD<DIM;++tD){
            for(int vD=0;vD<DIM;++vD){
                for(int wD=0;wD<DIM;++wD){
                    K[sD][tD][vD][wD]=0.0;
                    if(sD==tD && vD==wD){
                        K[sD][tD][vD][wD] += D * v;
                    }
                    if(sD==vD && tD==wD){
                        K[sD][tD][vD][wD] += D * (1-v);
                    }
                }
            }
        }
    }
}

static inline double weight(double r[DIM]){
    return 1 - (r[0]*r[0]+r[1]*r[1]+r[2]*r[2])/Radius/Radius;
}

static void calculateNormalizer()
{
    double (*a)[DIM] = InitialPosition;
    double (*_A)[DIM][DIM] = NormalizeTensor;
    
    for(int iP=0;iP<ParticleCount;++iP){
        double Ai[DIM][DIM] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
        for(int iN=0;iN<InitialNeighbor[iP].count;++iN){
            const int jP=InitialNeighbor[iP].stored[iN];
            double aij[DIM];
            for(int iD=0;iD<DIM;++iD){
                aij[iD] = a[jP][iD] - a[iP][iD];
            }
            double wij=weight(aij);
            for(int sD=0;sD<DIM;++sD){
                for(int tD=0;tD<DIM;++tD){
                    Ai[sD][tD]+=aij[sD]*aij[tD]*wij;
                }
            }
        }

        // rotate Ai to diagnal matrix with R
        double R[DIM][DIM]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
        jacobi(R,Ai);

        // calculate _Ai (inverse of Ai)
        int minD=0;
        for(int iD=0;iD<DIM;++iD){
            if(Ai[minD][minD]>Ai[iD][iD])minD=iD;
        }
        double _Ai[DIM][DIM]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
        for(int iD=0;iD<DIM;++iD){
            if(iD==minD){_Ai[iD][iD]=0.0;}
            else{        _Ai[iD][iD]=1.0/Ai[iD][iD];}
        }

        // re-rotate _Ai with transpose of R
        memcpy(_A[iP],&DOUBLE_ZERO,sizeof(double [DIM][DIM]));
        for(int sD=0;sD<DIM;++sD){
            for(int tD=0;tD<DIM;++tD){
                for(int uD=0;uD<DIM;++uD){
                    for(int vD=0;vD<DIM;++vD){
                        _A[iP][sD][vD]+=R[sD][tD]*_Ai[tD][uD]*R[vD][uD];
                    }
                }
            }
        }

        _A[iP][0][1]=_A[iP][1][0]=0.5*_A[iP][0][1]+0.5*_A[iP][1][0];
        _A[iP][1][2]=_A[iP][2][1]=0.5*_A[iP][1][2]+0.5*_A[iP][2][1];
        _A[iP][2][0]=_A[iP][0][2]=0.5*_A[iP][2][0]+0.5*_A[iP][0][2];
    }
}

static void resetAcceleration()
{
    for(int iP=0;iP<ParticleCount;++iP){
        for(int iD=0;iD<DIM;++iD){
            Acceleration[iP][iD]=0.0;
        }
    }
}

static void calculateDeformationGradient()
{
    double (*a)[DIM] = InitialPosition;
    double (*q)[DIM] = Position;
    double (*F)[DIM][DIM] = DeformationGradient;
    double (*_A)[DIM][DIM] = NormalizeTensor;

    for(int iP=0;iP<ParticleCount;++iP){
        double Gi[DIM][DIM];// Sum( qij aij wij ), F=G/A
        memcpy(Gi,DOUBLE_ZERO,sizeof(double [DIM][DIM]));
        for(int iN=0;iN<InitialNeighbor[iP].count;++iN){
            const int jP=InitialNeighbor[iP].stored[iN];
            double qij[DIM],aij[DIM],wij;
            for(int iD=0;iD<DIM;++iD){
                qij[iD]=q[jP][iD]-q[iP][iD];
                aij[iD]=a[jP][iD]-a[iP][iD];
            }
            wij=weight(aij);
            for(int sD=0;sD<DIM;++sD){
                for(int aD=0;aD<DIM;++aD){
                    Gi[sD][aD] += qij[sD]*aij[aD]*wij;
                }
            }
        }
        memcpy(F[iP],DOUBLE_ZERO,sizeof(double [DIM][DIM]));
        for(int sD=0;sD<DIM;++sD){
            for(int aD=0;aD<DIM;++aD){
                for(int tD=0;tD<DIM;++tD){
                    F[iP][sD][tD]+=Gi[sD][aD]*_A[iP][aD][tD];
                }
            }
        }
    }
}

static void calculateCurvatureTensor()
{
    double (*a)[DIM] = InitialPosition;
    double (*F)[DIM][DIM] = DeformationGradient;
    double (*C)[DIM][DIM][DIM] = CurvatureTensor;
    double (*_A)[DIM][DIM] = NormalizeTensor;

    for(int iP=0;iP<ParticleCount;++iP){
        double Di[DIM][DIM][DIM]; // Sum( Fij aij wij ), C=D/A
        memcpy(Di,DOUBLE_ZERO,sizeof(double [DIM][DIM][DIM]));
        for(int iN=0;iN<InitialNeighbor[iP].count;++iN){
            const int jP=InitialNeighbor[iP].stored[iN];
            double Fij[DIM][DIM],aij[DIM],wij;
            for(int iD=0;iD<DIM;++iD){
                for(int jD=0;jD<DIM;++jD){
                    Fij[iD][jD]=F[jP][iD][jD]-F[iP][iD][jD];
                }
                aij[iD]=a[jP][iD]-a[iP][iD];
            }
            wij=weight(aij);
            for(int rD=0;rD<DIM;++rD){
                for(int sD=0;sD<DIM;++sD){
                    for(int aD=0;aD<DIM;++aD){
                        Di[rD][sD][aD] += Fij[rD][sD]*aij[aD]*wij;
                    }
                }
            }
        }
        memcpy(C[iP],DOUBLE_ZERO,sizeof(double [DIM][DIM][DIM]));
        for(int rD=0;rD<DIM;++rD){
            for(int sD=0;sD<DIM;++sD){
                for(int aD=0;aD<DIM;++aD){
                    for(int tD=0;tD<DIM;++tD){
                        C[iP][rD][sD][tD]+=Di[rD][sD][aD]*_A[iP][aD][tD];
                    }
                }
            }
        }
    }
}

static void calculateMotionEquation()
{
    const double (*K)[DIM][DIM][DIM]=CoefficientTensor;
    const double (*C)[DIM][DIM][DIM]=CurvatureTensor;
    const double (*a)[DIM]=InitialPosition;
    const double (*_A)[DIM][DIM]=NormalizeTensor;
    // const double mass = Density*Thickness*ParticleSpacing*ParticleSpacing;
    // double (*v)[DIM]= Velocity;

    for(int iP=0;iP<ParticleCount;++iP){
        double Mi[DIM][DIM][DIM]; // Moment Tenssor, M=K*C
        memcpy(Mi,DOUBLE_ZERO,sizeof(double [DIM][DIM][DIM]));
        for(int rD=0;rD<DIM;++rD){
             for(int sD=0;sD<DIM;++sD){
                 for(int tD=0;tD<DIM;++tD){
                     for(int vD=0;vD<DIM;++vD){
                         for(int wD=0;wD<DIM;++wD){
                             Mi[rD][sD][tD] += K[sD][tD][vD][wD]*C[iP][rD][vD][wD];
                         }
                     }
                 }
             }
        }

        for(int jN=0;jN<InitialNeighbor[iP].count;++jN){
            const int jP=InitialNeighbor[iP].stored[jN];
            double aij[DIM];
            for(int cD=0;cD<DIM;++cD){
                aij[cD]=a[jP][cD]-a[iP][cD];
            }
            double wij=weight(aij);
            double aijwij_Ai[DIM]={0.0,0.0,0.0};
            for(int tD=0;tD<DIM;++tD){
                for(int cD=0;cD<DIM;++cD){
                    aijwij_Ai[tD] += aij[cD] * wij * _A[iP][cD][tD];
                }
            }

            // 粒子jと粒子lの相互作用 (interaction of "j" and "l")
            for(int lN=0;lN<InitialNeighbor[jP].count;++lN){
                const int lP=InitialNeighbor[jP].stored[lN];
                double ajl[DIM];
                for(int aD=0;aD<DIM;++aD){
                    ajl[aD]=a[lP][aD]-a[jP][aD];
                }
                double wjl=weight(ajl);
                double ajlwjl_Aj[DIM]={0.0,0.0,0.0};
                for(int sD=0;sD<DIM;++sD){
                    for(int aD=0;aD<DIM;++aD){
                        ajlwjl_Aj[sD] += ajl[aD] * wjl * _A[jP][aD][sD];
                    }
                }
                double dp[DIM]={0.0,0.0,0.0};
                for(int rD=0;rD<DIM;++rD){
                    for(int sD=0;sD<DIM;++sD){
                        for(int tD=0;tD<DIM;++tD){
                            dp[rD]+= -Mi[rD][sD][tD]*ajlwjl_Aj[sD]*aijwij_Ai[tD];
                        }
                    }
                }
                for(int rD=0;rD<DIM;++rD){
                    Acceleration[lP][rD]+= +dp[rD]/Density/Thickness;
                    Acceleration[jP][rD]+= -dp[rD]/Density/Thickness;
                }
            }
            
            // 粒子iと粒子kの相互作用 (interaction of "i" and "k")
            for(int kN=0;kN<InitialNeighbor[iP].count;++kN){
                const int kP=InitialNeighbor[iP].stored[kN];
                double aik[DIM];
                for(int bD=0;bD<DIM;++bD){
                    aik[bD]=a[kP][bD]-a[iP][bD];
                }
                double wik=weight(aik);
                double aikwik_Ai[DIM]={0.0,0.0,0.0};
                for(int sD=0;sD<DIM;++sD){
                    for(int bD=0;bD<DIM;++bD){
                        aikwik_Ai[sD] += aik[bD] * wik * _A[iP][bD][sD];
                    }
                }
                double dp[DIM]={0.0,0.0,0.0};
                for(int rD=0;rD<DIM;++rD){
                    for(int sD=0;sD<DIM;++sD){
                        for(int tD=0;tD<DIM;++tD){
                            dp[rD]+= -Mi[rD][sD][tD]*aikwik_Ai[sD]*aijwij_Ai[tD];
                        }
                    }
                }
                for(int rD=0;rD<DIM;++rD){
                    Acceleration[kP][rD]+= -dp[rD]/Density/Thickness;
                    Acceleration[iP][rD]+= +dp[rD]/Density/Thickness;
                }
            }
        }
    }
}

static void calculateExternalForce( void ){
    for(int iP=0;iP<ParticleCount;++iP){
        if(Property[iP]==FIXED_PARTICLE){
            for(int iD=0;iD<DIM;++iD){
                Acceleration[iP][iD]=0.0;
            }
        }
    }
}

static void calculateConvection()
{
    for(int iP=0;iP<ParticleCount;++iP){
        for(int iD=0;iD<DIM;++iD){
            Velocity[iP][iD] += Acceleration[iP][iD]*Dt;
            Position[iP][iD] += Velocity[iP][iD]*Dt;
        }
    }
}

static void writeEnergyFile(char *filename)
{
    static FILE *fp;
    static int init_flag=0;
    if(init_flag==0){
        fp = fopen(filename, "w");
        init_flag=1;
    }
    
    const double (*K)[DIM][DIM][DIM]=CoefficientTensor;
    const double (*C)[DIM][DIM][DIM]=CurvatureTensor;
    const double (*v)[DIM]=Velocity;
    
    BendingPotential =0.0;
    Kinetic = 0.0;
//    double AngularMomentum[DIM]={0.0,0.0,0.0};
    for(int iP=0;iP<ParticleCount;++iP){
        double Mi[DIM][DIM][DIM]; // Moment Tenssor, M=K*C
        memcpy(Mi,DOUBLE_ZERO,sizeof(double [DIM][DIM][DIM]));
        for(int rD=0;rD<DIM;++rD){
             for(int sD=0;sD<DIM;++sD){
                 for(int tD=0;tD<DIM;++tD){
                     for(int vD=0;vD<DIM;++vD){
                         for(int wD=0;wD<DIM;++wD){
                             Mi[rD][sD][tD] += K[sD][tD][vD][wD]*C[iP][rD][vD][wD];
                         }
                     }
                 }
             }
        }
        for(int rD=0;rD<DIM;++rD){
            for(int sD=0;sD<DIM;++sD){
                for(int tD=0;tD<DIM;++tD){
                    BendingPotential += 0.5 * Mi[rD][sD][tD] * C[iP][rD][sD][tD] * ParticleSpacing * ParticleSpacing;
                }
            }
        }

        for(int iD=0;iD<DIM;++iD){
            Kinetic += 0.5 * Density * v[iP][iD] * v[iP][iD] * Thickness * ParticleSpacing * ParticleSpacing;
        }

//        AngularMomentum[2] += Position[iP][0]*Velocity[iP][1] - Position[iP][1]*Velocity[iP][0];
//        AngularMomentum[0] += Position[iP][1]*Velocity[iP][2] - Position[iP][2]*Velocity[iP][1];
//        AngularMomentum[1] += Position[iP][2]*Velocity[iP][0] - Position[iP][0]*Velocity[iP][2];
        
    }
//    fprintf(fp, "%e  %e %e %e  %e %e %e\n",
//            Time,
//            BendingPotential, Kinetic, BendingPotential+Kinetic,
//            AngularMomentum[0], AngularMomentum[1], AngularMomentum[2]);

    fprintf(fp, "%e  %e %e %e\n",
            Time,
            BendingPotential, Kinetic, BendingPotential+Kinetic);
}



