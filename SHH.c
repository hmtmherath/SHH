/*Copyright [2015] [Neranjan Suranga Edirisinghe]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.*/

#include "Headers.h"
void main(int argc ,char** argv)
//int argc;  /* contains number of arguments */
////char **argv; /* contains the arguments themselves */
{
bool condition;
double GK2,GNa,mK2S,ENa,Eh,mhS,mK2,El,Ek,Pol,H,Gl,z1,Gh,t0,sqH;
double *restrict Y,*restrict u1,*restrict u2,*restrict u3,*restrict u4,*restrict results;
double Final_Time,Last,Sigma;
int Steps,i,j,k,ProbSize;
char *tokenB,*tokenA;
char *search = "=";
FILE *fp;
const gsl_rng_type * T;
gsl_rng * r;
gsl_rng_env_setup();
T = gsl_rng_default;
r = gsl_rng_alloc (T);
gsl_rng_set(r,(unsigned int) atoi(argv[1]));
char FileName[128];
ProbSize=4;
Y=malloc(ProbSize*sizeof(double)); 
u1=malloc(ProbSize*sizeof(double));
u2=malloc(ProbSize*sizeof(double));
u3=malloc(ProbSize*sizeof(double));
u4=malloc(ProbSize*sizeof(double));
for (i=0;i<ProbSize;i++)
{
Y[i]=0.0;u1[i]=0.0;u2[i]=0.0;u3[i]=0.0;u4[0]=0.0;
}
GNa = 105.0;GK2 = 30.0;Gh = 4.0;Gl = 8.0;ENa = 0.045;Ek = -0.07;Eh = -0.021;El = -0.046;Pol = 0.006;mK2S=  -0.0075;mhS = 0.038;z1=0.0,t0=0;
H = 1e-6;Final_Time = 1e-3;Steps=2000000;
Last = Final_Time/H + 1;
Y[0]=-0.034735512748588;
Y[1]=0.454706120407293;
Y[2]=0.032834744275582;
Y[3]=0.205951547789834;
Sigma=0.00001;
typedef struct{
double *Voltage;
double *HNa;
double *MH;
double *MK2;
}DSystem;

DSystem Results;

Results.Voltage=(double *) malloc(Steps*sizeof(double));
Results.HNa=(double *) malloc(Steps*sizeof(double));
Results.MH=(double *) malloc(Steps*sizeof(double));
Results.MK2=(double *) malloc(Steps*sizeof(double));

FILE *file = fopen ( "InitialConditions.txt", "r" );
if ( file != NULL )
{
  char line [ 128 ]; /* or other suitable maximum line size */
  while ( fgets ( line, sizeof line, file ) != NULL ) /* read a line */
  {
	     // Token will point to the part before the =.
              tokenA = strtok(line, search);
              // Token will point to the part after the =.
              tokenB = strtok(NULL, search);
    
              if (strcmp("V",tokenA)==0 )
              {
                 Y[0]=(double) atof(tokenB);
              }
	      else if (strcmp("hNa",tokenA)==0 )
              {
                 Y[1]=(double) atof(tokenB);
              }
              else if (strcmp("mh",tokenA)==0 )
              {
                 Y[2]=(double) atof(tokenB);
              }
              else if (strcmp("mK2",tokenA)==0 )
              {
                 Y[3]=(double) atof(tokenB);
              }
	 }
  fclose ( file );
	printf("Initial Conditions ::\nV %.15f\nhNA %.15f\nmh %.15f\nmK2 %.15f\n\n\n",Y[0],Y[1],Y[2],Y[3]);
}

file = fopen ( "RunSetup.txt", "r" );
if ( file != NULL )
{
  char line [ 128 ]; /* or other suitable maximum line size */
  while ( fgets ( line, sizeof line, file ) != NULL ) /* read a line */
  {
	tokenA = strtok(line, search);
	tokenB = strtok(NULL, search);

              if (strcmp("Time_Step",tokenA)==0 )
              {  
                 H=(double) atof(tokenB);
              }
              else if (strcmp("Num_Steps",tokenA)==0 )
              {  
                 Steps=(int) atof(tokenB);
              }
              else if (strcmp("Step_Time",tokenA)==0 )
              {  
                 Final_Time=(double) atof(tokenB);
              }
              else if (strcmp("Noise_STD",tokenA)==0 )
              {  
                 Sigma=(double) atof(tokenB);
              }
         } 
  fclose ( file );
        printf("Run Parameters ::\nTime_Step %f\nNum_Steps %d\nStep_Time %f\nNoise_STD %f\n\n\n",H,Steps,Final_Time,Sigma);
}

    

file = fopen ( "parameters.txt", "r" );
if ( file != NULL )
{
  char line [ 128 ]; /* or other suitable maximum line size */
  while ( fgets ( line, sizeof line, file ) != NULL ) /* read a line */
  {
    // Token will point to the part before the =.
         tokenA = strtok(line, search);
    // Token will point to the part after the =.
         tokenB = strtok(NULL, search);	
	
	if (strcmp("GNa",tokenA)==0 )
	{
		GNa=(double) atof(tokenB);
	}
	else if (strcmp("GK2",tokenA)==0 )
	{
		GK2=(double) atof(tokenB);
	}
	else if (strcmp("Gl",tokenA)==0 )
        {
                Gl=(double) atof(tokenB);
        }
	else if (strcmp("Gh",tokenA)==0 )
        {
                Gh=(double) atof(tokenB);
        }
	else if (strcmp("Eh",tokenA)==0 )
        {
                Eh=(double) atof(tokenB);
        }
	else if (strcmp("Ek",tokenA)==0 )
        {
                Ek=(double) atof(tokenB);
        }
	else if (strcmp("El",tokenA)==0 )
        {
                El=(double) atof(tokenB);
        }
	else if (strcmp("Pol",tokenA)==0 )
        {
                Pol=(double) atof(tokenB);
        }
	else if (strcmp("mK2S",tokenA)==0 )
        {
                mK2S=(double) atof(tokenB);
        }
	else if (strcmp("mhS",tokenA)==0 )
        {
                mhS=(double) atof(tokenB);
      	}
  }
  fclose ( file );
printf("System Parameters :: \nGNa %f\nmhS %f\nmK2S %f\nPol %f\nEl %f\nEk %f\nEh %f\nGh %f\nGl %f\nGK2 %f\n\nStarting Computation ::\n",GNa,mhS,mK2S,Pol,El,Ek,Eh,Gh,Gl,GK2);
}
sprintf(FileName,"Results-%d.txt",atoi(argv[1]));    
fp=fopen(FileName, "w");   

sqH=sqrt(H);
for (i=0;i<Steps;i++)
{
	for (j=0;j<Last;j++)
	{
		condition = true;
		while(condition) {
		z1=gsl_ran_gaussian (r,Sigma/sqH);
		rk4vec (ProbSize,t0,H,u1,u2,u3,Y,u4,GK2,GNa,ENa,Eh,mhS,El,Ek,Pol,H,Gl,z1,Gh,mK2S,f);
		if ( u4[2]<0.0 || u4[2]<0.0 || u4[2]<0.0 )
		{
			condition = true;
			printf("Value Too Large %d\n",z1);
		}
		else if ( u4[2]>1.0 || u4[2]>1.0 || u4[2]>1.0 )
		{
			condition = true;
			printf("Value Too Large %d\n",z1);
		}  
		else
		{
			condition = false;
		}	
		}//while

		Y=u4;
	}//j for loop
	for (k=0;k<ProbSize;k++)
	{
		Results.Voltage[i]=Y[0];
		Results.HNa[i]=Y[1];
		Results.MH[i]=Y[2];
		Results.MK2[i]=Y[3];
	}
fprintf(fp,"%d %.17g %.17g %.17g %.17g \n",i,Results.Voltage[i],Results.HNa[i],Results.MH[i],Results.MK2[i]);
//fprintf(fp,"%d %f %f %f %f \n",i,Results.Voltage[i],Results.HNa[i],Results.MH[i],Results.MK2[i]);
}//i for loop
fclose(fp);
free ( Y );
Y=NULL;
free ( u1 );
u1=NULL;
free ( u2 );
u2=NULL;
free ( u3 );
u3=NULL;
free (u4);
u4=NULL;
gsl_rng_free (r);
}//main function

