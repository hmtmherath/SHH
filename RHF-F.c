double *f(double Y[],double GK2, double GNa, double ENa, double Eh, double mhS, double El, double Ek,double Pol,double H,double Gl,double z1,double Gh,double mK2S)
{
double Y1,Y2,Y3,Y4,mNa,mK2;
double *ydot;
int i;
ydot = ( double * ) malloc ( 4 * sizeof ( double ) );
for (i=0;i<4;i++)
ydot[i]=0.0;

Y1=Y[0];Y2=Y[1];Y3=Y[2];Y4=Y[3];
mNa=1.0/(1.0+exp(-150.0*(Y1+0.0305)));
mK2=1.0/(1.0+exp(-83.0*(Y1+mK2S)));
ydot[0]=-2.0*(GK2*Y4*Y4*(Y1-Ek)+GNa*mNa*mNa*mNa*Y2*(Y1-ENa)+Gh*Y3*Y3*(Y1-Eh)+Gl*(Y1-El)+Pol);
ydot[1]=(((1.0/(1.0+exp(500.0*(Y1+0.0325))))-Y2)+(z1/sqrt(H)))/0.0405 ;
ydot[2]=((1.0/(1.0+2.0*exp(180.0*(Y1+mhS))+exp(500.0*(Y1+mhS))))-Y3)/0.1;
ydot[3]=(mK2-Y4)/2.0 ;
return ydot;
}
