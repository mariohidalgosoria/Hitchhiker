/*The program gives the final sizes for an ensemble of different 
systems  of interacting random walkers on a 1D lattice, 
and for different diffusion rate models
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#define GNUPLOT "gnuplot -persist"
#define N 6000 //number of cells in the Lattice
#define P 6000//number of particles;  
float ran2(long *idum);
int ProbSampleReplace(int n, double *p);
void PolyDiff(int k, double *p);//Diffusion of a RW with mass M


int LA[N]={0};//array of search in each step

int x=0;//initial position of one random walker

int main(void) 
{

  int i,j,l,k,s,id,T;
  double r,u,prob,b;
  long iseed;
  double t, w,d,wb;
  double pib[2]={0,0};//prob of breaking=pib[0]; prob of diffusion=pib[1]
  double piw[3]={0,0,0};//prob of breaking=pib[0]; prob of diffusion=pib[1]
  int ini=1;
  int fin=100;//number of different systems simulated
  FILE *dat=fopen("HK3XIA.dat","w");//data file of the initial condition in the lattice, just for checking
  FILE *datos=fopen("HK3XTA.dat","w");//data file for the final sizes in the lattice


  srand(time(NULL));
  prob=0.5;//prob of walking
  u=0; 
  s=0;  
  t=0;
  w=0;
  b=0;
  r=0;
  
 
  wb=0.005;//rate of breaking
  //the total monte carlo steps achieved are =(T/0.01)/M 
  //with M=1000, for Rouse; M=500, for Zimm; M=2000 for Arrhenius
  //i.e. every M updates a monte carlo step is given
  //T=10000;//Rouse
  //T=5000;//Zimm
  T=20000;//Arrhenius
  //loop on the different systems of interacting RW's 
  for(id=ini;id<=fin;id++) {
  
   
    for(i=0;i<N;i++){//initialize in each case the lattice with monomers
     LA[i]=1;
    }
    x=2999;//initial position for a tracked particle, the middle of the lattice
    if(id ==fin){//in the final system save the initial data of N
      for(k=0;k<N;k++){ 
	if(LA[k]!=0) fprintf(dat,"%d\n",LA[k]);    
      }
    }
    t=0;//restart the time 
    while(t<T){//iteration over the discrete time steps Jumps
      /*select a place on the lattice*/
      iseed=-(long) rand();
      j= (int) (ran2(&iseed)*N);
      while(LA[j]==0){
	j= (int) (ran2(&iseed)*N);
      }
      //d=1.0/pow(LA[j],1); //Rouse
      //d=1.0/pow(LA[j],0.6); //Zimm
      d=exp(-LA[j]); //arrhenius
      iseed=-(long) rand();
      w=ran2(&iseed);//U r.v. for time t 
      
      pib[0]=wb/(wb+d);//probability of breaking
      pib[1]=d/(wb+d);//probability of diffusion
      
      s=ProbSampleReplace(2,pib);//select breaking or diffusion
      if(s==1){//do breaking
	if(LA[j]>1){//mass at last 2
	  if((LA[j]%2)==0){//if LA[j] is even
	    if(j==0){//BC there is breaking but j=0
	      LA[j+1]=LA[j+1]+(LA[j]/2);
	      LA[j]=LA[j]/2;
	    }//end else BC j=0
	    else if(j==(N-1)){//BC breaking j=N-1
	      LA[j-1]=LA[j-1]+(LA[j]/2);
	      LA[j]=LA[j]/2;
	    }//end else BC j=N-1
	    else {//j!=0 || j!=(N-=1
	      LA[j+1]=LA[j+1]+(LA[j]/2);
	      LA[j-1]=LA[j-1]+(LA[j]/2);
	      LA[j]=0;
	      if(j==x){
		iseed=-(long) rand();
		r=ran2(&iseed);
		//genre u Unif
		if(r>prob) x+=1;
		else x-=1;
	      }
	    }//end else j!=0
	  }//end if LA[j] is even
	
	  else{//else LA[j] is odd
	    iseed=-(long) rand();
	    b=ran2(&iseed);
	    if(b<=0.5){//if place heavier to the right 
		
	  
	      if(j==0){//BC there is breaking but j=0
		LA[j+1]=LA[j+1]+((LA[j]-1)/2);
		LA[j]=LA[j]-((LA[j]-1)/2);
	      }//end else BC j=0
	      else if(j==(N-1)){//BC breaking j=N-1
		LA[j-1]=LA[j-1]+((LA[j]-1)/2);
		LA[j]=LA[j]-((LA[j]-1)/2);
	      }//end else BC j=N-1
	      else {//j!=0 || j!=(N-=1
		LA[j+1]=LA[j+1]+((LA[j]-1)/2);
		LA[j-1]=LA[j-1]+(LA[j]-((LA[j]-1)/2));
		LA[j]=0;
		if(j==x) {
		  x+=1;
		}
	      }//end else j!=0
	      
	    }// end if  place heavier to the right
	      
	    else{//place heavier to the left
	      if(j==0){//BC there is breaking but j=0
		LA[j+1]=LA[j+1]+(LA[j]-((LA[j]-1)/2));
		LA[j]=((LA[j]-1)/2);
		if(j==x){
		  x+=1;
		}
	      }//end else BC j=0
	      else if(j==(N-1)){//BC breaking j=N-1
	
		LA[j-1]=LA[j-1]+(LA[j] -  ((LA[j]-1)/2));
		LA[j]=((LA[j]-1)/2);
		if(j==x){
		  x-=1;
		}
	      }//end else BC j=N-1
	      else {//j!=0 || j!=(N-=1
		LA[j+1]=LA[j+1] + (LA[j]-((LA[j]-1)/2));
		LA[j-1]=LA[j-1]+((LA[j]-1)/2);
		LA[j]=0;
		if(j==x) {
		  x-=1;
		}
	      }//end else j!=0
	      
	    }//end else heavier to the right
	    
	  
	  }//end else LA[j] is odd
	  
	}//end if LA[j]>1
	
      }//end if s==1, breaking
    
      else if(s==2){//do diffusion      
	if(LA[j]==1) {
	  iseed=-(long) rand();
	  u=ran2(&iseed);
	  if(u<=prob){//move to the right
	    if(j==(N-1)){//periodic B.C. 
	      LA[0]=LA[j]+LA[0];//agrregation
	      LA[j]=0;
	      if(j==x){ 
		x+=1;
	      }
	    }//end if Periodic BC
	    else{//do the step to the right
	      LA[j+1]=LA[j]+LA[j+1];//aggregation
	      LA[j]=0;
	      if(j==x) {
		x+=1;
	      }
	    }//end else step to the right
	  }//end if u<=prob
	  else {//move to the left
	    if(j==0){//periodic B.C.
	      LA[N-1]=LA[j]+LA[N-1];//agrregation
	      LA[j]=0;
	      if(j==x){ 
		x-=1;
	      }
	    }//end if Periodic B.C.
	    else{//do the step to the left
	      LA[j-1]=LA[j]+LA[j-1];//aggregation
	      LA[j]=0;
	      if(j==x){ 
		x-=1;
	      }
	    }//end else step to the left
	  }//end else move to the left
	  
	}//end  if(LA[j]==1)
	
	else{//LA[] is not a monomer
	 
	  //piw[0]=piw[2]=(double) 1.0/(2*LA[j]);//Rouse 
	  //piw[1]=(double) (LA[j]-1)/LA[j];//Rouse
	  //piw[0]=piw[2]=(double) 1.0/(2*pow(LA[j],0.6));//Zimm
	  //piw[1]=(double) (pow(LA[j],0.6)-1)/pow(LA[j],0.6);//Zimm
	  piw[0]=piw[2]=(double) exp(-LA[j])/2;//Arrhenius
	  piw[1]=(double) 1-exp(-LA[j]);//Arrhenius
	  PolyDiff(j,piw);
	 
	}//end else is not a monomer
	
      
      }//end else if s==2
      
      t+=0.01;//advance in inner time
     
    }//end while
    
    for(k=0;k<N;k++){ //for each system save the final sizes
      if(LA[k]!=0) fprintf(datos,"%d\n",LA[k]);
    }
  }//end for id 
 
  fclose(dat);
  fclose(datos);
  
  return 0;
}  




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

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
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
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

int ProbSampleReplace(int n, double *p)
{
  float rU;
  int i, j,ans;
  int nm1 = n - 1;
  int perm[N]={0};
  double q[N]={0};
  long iseed;
  /* record element identities */
  for (i = 0; i < n; i++) {
    perm[i] = i + 1;
    q[i]=p[i];
  }
  /* compute cumulative probabilities */
  for (i = 1 ; i < n; i++)
    q[i] += q[i - 1];
  /* compute the sample */
  iseed=-(long) rand();
  rU = ran2(&iseed);
  for (j = 0; j < nm1; j++) {
    if (rU <= q[j])
      break;
  }
  //en C despues del for j=nm1
  ans = perm[j];
  return ans;
}



void PolyDiff(  int k, double *p)
{

  int s;
 
  s=ProbSampleReplace(3,p);
  if(s==1) {//jump to the left
    // printf("%d jump to the left %d\n",LA[k],s);
    if(k==0){//periodic B.C.
      LA[N-1]=LA[k]+LA[N-1];
      LA[k]=0;
      if(k==x){ 
	x-=1;
	//printf("L[%d]=%d\n",x,LA[x]);
      }   
    }//end if Reflc B.C.
    else{//do the step to the left
      LA[k-1]=LA[k-1]+LA[k];
      LA[k]=0;
      if(k==x){
	x-=1;
	//printf("L[%d]=%d\n",x,LA[x]);
      }
    }//end else step to the left
  }//end if step to the left s=1 
  else if(s==2){//no jump
    //printf("%d don't jump %d\n",LA[k],s);
    LA[k]=LA[k];
    if(k==x) {
      //x=x;
      //printf("L[%d]=%d\n",x,LA[x]);   
    } 
  }//end else if no jump s=2
  else if(s==3){//jump to the right
    //printf("%d jump to the right %d\n",LA[k],s);
    if(k==(N-1)){//periodic B.C. 
      LA[0]=LA[0]+LA[k];
      LA[k]=0;
      if(k==x){ 
	x+=1;
	//printf("L[%d]=%d\n",x,LA[x]);
      }
    }//end if Refl BC
    else{//do the step to the right
      LA[k+1]=LA[k+1]+LA[k];
      LA[k]=0;
      if(k==x){ 
	x+=1;
	//printf("L[%d]=%d\n",x,LA[x]);
      }
    }//end else step to the right
  }//end else if step to the right s=3
  
}



