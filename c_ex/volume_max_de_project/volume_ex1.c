// 必須在演算過程中, 設法限制各變數的上下限!
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <time.h>
   
// 最大族群數, NP
#define MAXPOP  5000
// 最大向量維度, D
#define MAXDIM  35
// MAXIMAPROBLEM =1 最大化 0 最小化
#define MAXIMAPROBLEM 1
// 最大化時 PENALITY 必須為負值, 否則為正值
#define PENALITY -1000
/*
#define MAXIMAPROBLEM 0
#define PENALITY 1000
*/
   
/*------Constants for rnd_uni()--------------------------------------------*/
   
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
   
/*------------------------Globals---------------------------------------*/
   
long  rnd_uni_init;                 /* serves as a seed for rnd_uni()   */
double c[MAXPOP][MAXDIM], d[MAXPOP][MAXDIM];
double (*pold)[MAXPOP][MAXDIM], (*pnew)[MAXPOP][MAXDIM], (*pswap)[MAXPOP][MAXDIM];
   
/*---------Function declarations----------------------------------------*/
   
void  assignd(int D, double a[], double b[]);
double rnd_uni(long *idum);    /* uniform pseudo random number generator */
double extern evaluate(int D, double tmp[], long *nfeval); /* obj. funct. */
   
/*---------Function definitions-----------------------------------------*/
// 指定向量 b 為 a
void  assignd(int D, double a[], double b[])
{
   int j;
   for (j=0; j<D; j++)
   {
      a[j] = b[j];
   }
}
   
// 產生 0 ~ 1 間的亂數
double rnd_uni(long *idum)
{
  long j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;
   
  if (*idum <= 0)
  {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--)
    {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
   
}/*------End of rnd_uni()--------------------------*/
   
// 將上下限轉為全域變數
double inibound_h;      /* upper parameter bound              */
double inibound_l;      /* lower parameter bound              */
// 與機構合成相關的全域變數
// 宣告一個座標結構
struct Coord {
    double x;
    double y;
  // 這裡保留 double z;
};
   
main(int argc, char *argv[])
{
   char  chr;             /* y/n choice variable                */
   char  *strat[] =       /* strategy-indicator                 */
   {
            "",
            "DE/best/1/exp",
            "DE/rand/1/exp",
            "DE/rand-to-best/1/exp",
            "DE/best/2/exp",
            "DE/rand/2/exp",
            "DE/best/1/bin",
            "DE/rand/1/bin",
            "DE/rand-to-best/1/bin",
            "DE/best/2/bin",
            "DE/rand/2/bin"
   };
   
   int   i, j, L, n;      /* counting variables                 */
   int   r1, r2, r3, r4;  /* placeholders for random indexes    */
   int   r5;              /* placeholders for random indexes    */
   int   D;               /* Dimension of parameter vector      */
   int   NP;              /* number of population members       */
   int   imin;            /* index to member with lowest energy */
   int   refresh;         /* refresh rate of screen output      */
   int   strategy;        /* choice parameter for screen output */
   int   gen, genmax, seed;   
   
   long  nfeval;          /* number of function evaluations     */
   
   double trial_cost;      /* buffer variable                    */
   // 將上下限轉為全域變數, 可能要根據各變數加以設定
   //double inibound_h;      /* upper parameter bound              */
   //double inibound_l;      /* lower parameter bound              */
   double tmp[MAXDIM], best[MAXDIM], bestit[MAXDIM]; /* members  */
   double cost[MAXPOP];    /* obj. funct. values                 */
   double cvar;            /* computes the cost variance         */
   double cmean;           /* mean cost                          */
   double F,CR;            /* control variables of DE            */
   double cmin;            /* help variables                     */
   
   FILE  *fpin_ptr;
   FILE  *fpout_ptr;
   
// 計算執行過程所需時間起點, 需要導入 time.h
  clock_t start = clock();
   
/*------Initializations----------------------------*/
   
// 將結果寫入 out.dat
 fpout_ptr = fopen("out.dat","w");          /* open output file for reading,    */
// 目前已經採用 strategy 3 可以得到最佳結果
  strategy = 3;
  genmax = 2000;
  refresh = 100;
  // 配合機構尺寸合成, 每一個體有 9 個機構尺寸值與 5 個通過點角度值
  D = 2;
  NP = 200;
  inibound_h = 50.;
  inibound_l = 0.;
/*得到最佳解
  F = 0.85;
CR 必須介於 0 to 1. 之間
  CR = 1.;
*/
  F = 0.85;
  CR = 1.;
  seed = 3;
   
 //fclose(fpin_ptr);
   
/*-----Checking input variables for proper range----------------------------*/
   
  if (D > MAXDIM)
  {
     printf("\nError! D=%d > MAXDIM=%d\n",D,MAXDIM);
     exit(1);
  }
  if (D <= 0)
  {
     printf("\nError! D=%d, should be > 0\n",D);
     exit(1);
  }
  if (NP > MAXPOP)
  {
     printf("\nError! NP=%d > MAXPOP=%d\n",NP,MAXPOP);
     exit(1);
  }
  if (NP <= 0)
  {
     printf("\nError! NP=%d, should be > 0\n",NP);
     exit(1);
  }
  if ((CR < 0) || (CR > 1.0))
  {
     printf("\nError! CR=%f, should be ex [0,1]\n",CR);
     exit(1);
  }
  if (seed <= 0)
  {
     printf("\nError! seed=%d, should be > 0\n",seed);
     exit(1);
  }
  if (refresh <= 0)
  {
     printf("\nError! refresh=%d, should be > 0\n",refresh);
     exit(1);
  }
  if (genmax <= 0)
  {
     printf("\nError! genmax=%d, should be > 0\n",genmax);
     exit(1);
  }
  if ((strategy < 0) || (strategy > 10))
  {
     printf("\nError! strategy=%d, should be ex {1,2,3,4,5,6,7,8,9,10}\n",strategy);
     exit(1);
  }
  if (inibound_h < inibound_l)
  {
     printf("\nError! inibound_h=%f < inibound_l=%f\n",inibound_h, inibound_l);
     exit(1);
  }
   
/*-----Initialize random number generator-----------------------------*/
   
 rnd_uni_init = -(long)seed;  /* initialization of rnd_uni() */
 nfeval       =  0;  /* reset number of function evaluations */
   
/*------Initialization------------------------------------------------*/
/*------Right now this part is kept fairly simple and just generates--*/
/*------random numbers in the range [-initfac, +initfac]. You might---*/
/*------want to extend the init part such that you can initialize-----*/
/*------each parameter separately.------------------------------------*/
   
   for (i=0; i<NP; i++)
   {
      for (j=0; j<D; j++) /* spread initial population members */
      {
        c[i][j] = inibound_l + rnd_uni(&rnd_uni_init)*(inibound_h - inibound_l);
      }
      cost[i] = evaluate(D,c[i],&nfeval); /* obj. funct. value */
   }
   cmin = cost[0];
   imin = 0;
   for (i=1; i<NP; i++)
   {
     if(MAXIMAPROBLEM == 1)
     {
       // 改為最大化
        if (cost[i]>cmin)
        {
          cmin = cost[i];
          imin = i;
        }
      }
      else
      {
        // 最小化問題
        if (cost[i]<cmin)
        {
          cmin = cost[i];
          imin = i;
        }
      }
   }
   
   assignd(D,best,c[imin]);            /* save best member ever          */
   assignd(D,bestit,c[imin]);          /* save best member of generation */
   
   pold = &c; /* old population (generation G)   */
   pnew = &d; /* new population (generation G+1) */
   
/*=======================================================================*/
/*=========Iteration loop================================================*/
/*=======================================================================*/
   
   gen = 0;                          /* generation counter reset */
   while ((gen < genmax) /*&& (kbhit() == 0)*/) /* remove comments if conio.h */
   {                                            /* is accepted by compiler    */
      gen++;
      imin = 0;
   
      for (i=0; i<NP; i++)         /* Start of loop through ensemble  */
      {
     do                        /* Pick a random population member */
     {                         /* Endless loop for NP < 2 !!!     */
       r1 = (int)(rnd_uni(&rnd_uni_init)*NP);
     }while(r1==i);            
   
     do                        /* Pick a random population member */
     {                         /* Endless loop for NP < 3 !!!     */
       r2 = (int)(rnd_uni(&rnd_uni_init)*NP);
     }while((r2==i) || (r2==r1));
   
     do                        /* Pick a random population member */
     {                         /* Endless loop for NP < 4 !!!     */
       r3 = (int)(rnd_uni(&rnd_uni_init)*NP);
     }while((r3==i) || (r3==r1) || (r3==r2));
   
     do                        /* Pick a random population member */
     {                         /* Endless loop for NP < 5 !!!     */
       r4 = (int)(rnd_uni(&rnd_uni_init)*NP);
     }while((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));
   
     do                        /* Pick a random population member */
     {                         /* Endless loop for NP < 6 !!!     */
       r5 = (int)(rnd_uni(&rnd_uni_init)*NP);
     }while((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));
   
/*=======EXPONENTIAL CROSSOVER============================================================*/
   
/*-------DE/best/1/exp--------------------------------------------------------------------*/
/*-------Our oldest strategy but still not bad. However, we have found several------------*/
/*-------optimization problems where misconvergence occurs.-------------------------------*/
     if (strategy == 1) /* strategy DE0 (not in our paper) */
     {
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D);
       L = 0;
       do
       {                       
         tmp[n] = bestit[n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
         n = (n+1)%D;
         L++;
       }while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
     }
/*-------DE/rand/1/exp-------------------------------------------------------------------*/
/*-------This is one of my favourite strategies. It works especially well when the-------*/
/*-------"bestit[]"-schemes experience misconvergence. Try e.g. F=0.7 and CR=0.5---------*/
/*-------as a first guess.---------------------------------------------------------------*/
     else if (strategy == 2) /* strategy DE1 in the techreport */
     {
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D);
       L = 0;
       do
       {                       
         tmp[n] = (*pold)[r1][n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
         n = (n+1)%D;
         L++;
       }while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
     }
/*-------DE/rand-to-best/1/exp-----------------------------------------------------------*/
/*-------This strategy seems to be one of the best strategies. Try F=0.85 and CR=1.------*/
/*-------If you get misconvergence try to increase NP. If this doesn't help you----------*/
/*-------should play around with all three control variables.----------------------------*/
     else if (strategy == 3) /* similiar to DE2 but generally better */
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
       L = 0;
       do
       {                       
         tmp[n] = tmp[n] + F*(bestit[n] - tmp[n]) + F*((*pold)[r1][n]-(*pold)[r2][n]);
         n = (n+1)%D;
         L++;
       }while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
     }
/*-------DE/best/2/exp is another powerful strategy worth trying--------------------------*/
     else if (strategy == 4)
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
       L = 0;
       do
       {                           
         tmp[n] = bestit[n] + 
              ((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
         n = (n+1)%D;
         L++;
       }while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
     }
/*-------DE/rand/2/exp seems to be a robust optimizer for many functions-------------------*/
     else if (strategy == 5)
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
       L = 0;
       do
       {                           
         tmp[n] = (*pold)[r5][n] + 
              ((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
         n = (n+1)%D;
         L++;
       }while((rnd_uni(&rnd_uni_init) < CR) && (L < D));
     }
   
/*=======Essentially same strategies but BINOMIAL CROSSOVER===============================*/
   
/*-------DE/best/1/bin--------------------------------------------------------------------*/
     else if (strategy == 6) 
     {
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
           for (L=0; L<D; L++) /* perform D binomial trials */
           {
         if ((rnd_uni(&rnd_uni_init) < CR) || L == (D-1)) /* change at least one parameter */
         {                       
           tmp[n] = bestit[n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
         }
         n = (n+1)%D;
           }
     }
/*-------DE/rand/1/bin-------------------------------------------------------------------*/
     else if (strategy == 7) 
     {
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
           for (L=0; L<D; L++) /* perform D binomial trials */
           {
         if ((rnd_uni(&rnd_uni_init) < CR) || L == (D-1)) /* change at least one parameter */
         {                       
           tmp[n] = (*pold)[r1][n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
         }
         n = (n+1)%D;
           }
     }
/*-------DE/rand-to-best/1/bin-----------------------------------------------------------*/
     else if (strategy == 8) 
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
           for (L=0; L<D; L++) /* perform D binomial trials */
           {
         if ((rnd_uni(&rnd_uni_init) < CR) || L == (D-1)) /* change at least one parameter */
         {                       
           tmp[n] = tmp[n] + F*(bestit[n] - tmp[n]) + F*((*pold)[r1][n]-(*pold)[r2][n]);
         }
         n = (n+1)%D;
           }
     }
/*-------DE/best/2/bin--------------------------------------------------------------------*/
     else if (strategy == 9)
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
           for (L=0; L<D; L++) /* perform D binomial trials */
           {
         if ((rnd_uni(&rnd_uni_init) < CR) || L == (D-1)) /* change at least one parameter */
         {                       
           tmp[n] = bestit[n] + 
              ((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
         }
         n = (n+1)%D;
           }
     }
/*-------DE/rand/2/bin--------------------------------------------------------------------*/
     else
     { 
       assignd(D,tmp,(*pold)[i]);
       n = (int)(rnd_uni(&rnd_uni_init)*D); 
           for (L=0; L<D; L++) /* perform D binomial trials */
           {
         if ((rnd_uni(&rnd_uni_init) < CR) || L == (D-1)) /* change at least one parameter */
         {                       
           tmp[n] = (*pold)[r5][n] + 
              ((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
         }
         n = (n+1)%D;
           }
     }
   
   
/*=======Trial mutation now in tmp[]. Test how good this choice really was.==================*/
   
     trial_cost = evaluate(D,tmp,&nfeval);  /* Evaluate new vector in tmp[] */
   if(MAXIMAPROBLEM == 1)
   {
    // 改為最大化
       if (trial_cost >= cost[i])   /* improved objective function value ? */
       {                                  
          cost[i]=trial_cost;         
          assignd(D,(*pnew)[i],tmp);
          if (trial_cost>cmin)          /* Was this a new minimum? */
          {                               /* if so...*/
             cmin=trial_cost;           /* reset cmin to new low...*/
             imin=i;
             assignd(D,best,tmp);           
          }                           
       }                            
       else
       {
          assignd(D,(*pnew)[i],(*pold)[i]); /* replace target with old value */
       }
    }
    else
    {
          // 最小化問題
       if (trial_cost <= cost[i])   /* improved objective function value ? */
       {                                  
          cost[i]=trial_cost;         
          assignd(D,(*pnew)[i],tmp);
          if (trial_cost<cmin)          /* Was this a new minimum? */
          {                               /* if so...*/
             cmin=trial_cost;           /* reset cmin to new low...*/
             imin=i;
             assignd(D,best,tmp);           
          }                           
       }                            
       else
       {
          assignd(D,(*pnew)[i],(*pold)[i]); /* replace target with old value */
       }
    }
   
      }   /* End mutation loop through pop. */
   
      assignd(D,bestit,best);  /* Save best population member of current iteration */
   
      /* swap population arrays. New generation becomes old one */
   
      pswap = pold;
      pold  = pnew;
      pnew  = pswap;
   
/*----Compute the energy variance (just for monitoring purposes)-----------*/
   
      cmean = 0.;          /* compute the mean value first */
      for (j=0; j<NP; j++)
      {
         cmean += cost[j];
      }
      cmean = cmean/NP;
   
      cvar = 0.;           /* now the variance              */
      for (j=0; j<NP; j++)
      {
         cvar += (cost[j] - cmean)*(cost[j] - cmean);
      }
      cvar = cvar/(NP-1);
   
   
/*----Output part----------------------------------------------------------*/
   
      if (gen%refresh==1)   /* display after every refresh generations */
      { /* ABORT works only if conio.h is accepted by your compiler */
    printf("\n\n                         PRESS ANY KEY TO ABORT"); 
    printf("\n\n\n Best-so-far cost funct. value=%-15.10g\n",cmin);
   
    for (j=0;j<D;j++)
    {
      printf("\n best[%d]=%-15.10g",j,best[j]);
    }
    printf("\n\n Generation=%d  NFEs=%ld   Strategy: %s    ",gen,nfeval,strat[strategy]);
    printf("\n NP=%d    F=%-4.2g    CR=%-4.2g   cost-variance=%-10.5g\n",
               NP,F,CR,cvar);
      }
   
      fprintf(fpout_ptr,"%ld   %-15.10g\n",nfeval,cmin);
   }
/*=======================================================================*/
/*=========End of iteration loop=========================================*/
/*=======================================================================*/
   
/*-------Final output in file-------------------------------------------*/
   
   
   fprintf(fpout_ptr,"\n\n\n Best-so-far obj. funct. value = %-15.10g\n",cmin);
   
   for (j=0;j<D;j++)
   {
     fprintf(fpout_ptr,"\n best[%d]=%-15.10g",j,best[j]);
   }
   fprintf(fpout_ptr,"\n\n Generation=%d  NFEs=%ld   Strategy: %s    ",gen,nfeval,strat[strategy]);
   fprintf(fpout_ptr,"\n NP=%d    F=%-4.2g    CR=%-4.2g    cost-variance=%-10.5g\n",
           NP,F,CR,cvar); 
   
  fclose(fpout_ptr);
   
  /* Code you want timed here */
  printf("Time elapsed: %f\n", ((double)clock() - start) / CLOCKS_PER_SEC);
   return(0);
}
   
/*-----------End of main()------------------------------------------*/
   
// 適應函式 fittness function (cost function)
double evaluate(int D, double tmp[], long *nfeval)
{
   double result=0, surface = 192.0, z, volume, penality;
   (*nfeval)++;
   z = (surface-tmp[0]*tmp[1])/(2.0*(tmp[0]+tmp[1]));
   volume = tmp[0]*tmp[1]*z;
   
  if(volume <= 0){
    return PENALITY;
  }
// 只限制長度與寬度必須大於 0
  if(tmp[0] <= inibound_l){
    return PENALITY;
  }
   
  if(tmp[1] <= inibound_l){
    return PENALITY;
  }
/*
  if((tmp[0] <= inibound_l)|| (tmp[0] >inibound_h)){
    return PENALITY;
  }
   
  if((tmp[1] <= inibound_l) || (tmp[1] >inibound_h)){
    return PENALITY;
  }
  */
  // volume must >0 and max volume
  // 目前為最小化問題
   //return 1+1/(volume*volume);
  return volume;
}