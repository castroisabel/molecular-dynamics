#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define npart 50
#define tmax 200000
#define densi 0.2
#define temp 0.5
#define caixa pow(npart/densi,0.5)
#define dt 0.01
#define alfa 6.0
#define kapa 2.0
#define Ao 0.0
#define rcri caixa/2

void posicaoinicial(double *xx,double *yy);
void velocidadeinicial(double *vxx,double *vyy);
double force(double *fxx,double*fyy,double *xx,double *yy,double energia);
void integrate(double *fxx,double *fyy,double *xx,double *yy,double *vxx,double *vyy,double *vvxx,double *vvyy);
double integraVEL(double *fxx,double *fyy,double *vxx,double *vyy,double *vvxx,double *vvyy, double cinetica);
double termostato(double *vxx,double *vyy, double cinetica);
void CPC(double *xx,double *yy);

int main(void){

	clock_t t1,t2;
	t1 = clock();

	// criando a variável ponteiro para os arquivos
	FILE *file1, *file2, *file3, *file4;

	file1 = fopen("posi.dat", "w");	//abrindo os arquivos
	file2 = fopen("blender.xyz", "w");
	file3 = fopen("posf.dat", "w");
	file4 = fopen("Energy.dat", "w");


	// variáveis
	double x[npart], y[npart], vx[npart], vy[npart], fx[npart], fy[npart], vvx[npart], vvy[npart];
	double Ep=0.0,Ec=0.0;
	for(int i=0;i<npart;i++){
	x[i]=0; y[i]=0; vx[i]=0; vy[i]=0; fx[i]=0; fy[i]=0; vvx[i]=0; vvy[i]=0;  // zerando os vetores
		}


	int tf = 100; //quantidade de frames
	int tmod=tmax/tf;


	// início
	posicaoinicial(x,y);
	velocidadeinicial(vx,vy);
	force(fx,fy,x,y,Ep);

	fprintf(file2,"%d \n",npart);
    fprintf(file2,"Centro da partícula \n");
	for(int j=0;j<npart;j++){
		fprintf(file2,"Default      %lf      %lf      0 \n",x[j],y[j]);	// blender
		fprintf(file1,"%lf      %lf\n",x[j],y[j]);	// posição inicial
		}

	for(int t=1;t<=tmax;t++){
		integrate(fx,fy,x,y,vx,vy,vvx,vvy);
		CPC(x,y);
		Ep = force(fx,fy,x,y,Ep);
		Ec = integraVEL(fx,fy,vx,vy,vvx,vvy,Ec);
		Ec = termostato(vx,vy,Ec);
		printf("%d %lf %lf %lf\n",t,Ep,Ec,Ep+Ec);
		if(t%tmod==0){
			fprintf(file2,"%d \n",npart);
			fprintf(file2,"Centro da partícula \n");
			for(int j=0;j<npart;j++){
				fprintf(file2,"Default      %lf      %lf      0 \n",x[j],y[j]);    // blender
			}
		}
		fprintf(file4,"%d      %lf      %lf \n",t,Ep,Ep+Ec);  // energias

	}


	for(int j=0;j<npart;j++){
		fprintf(file3,"%lf      %lf\n",x[j],y[j]);	// posição final
	}

	t2 = clock();
	printf("%f",(((double)t2 - (double)t1) / 1000000.0F ));
	return 0;
}


// rede quadrada
void posicaoinicial(double *x,double *y){
	double a=(caixa)/(pow((double)npart,0.5)+1);
	double ax=a, ay=a;
	for(int i=0;i<npart;i++){
		x[i]=(-caixa*0.50)+ax; y[i]=(-caixa*0.50)+ay;
		ax=ax+a;
	if (ax>=caixa){ ax=a; ay=ay+a; }
	}
 }

 // velocidade aleatória
 void velocidadeinicial(double *vx,double *vy){
	double sumvx,sumvy,sumv2,fs;
	sumvx=sumvy=sumv2=fs=0.0;

	srand(time(NULL));
	for(int i=0;i<npart;i++){
		vx[i]=2.*rand()/RAND_MAX-1;
		vy[i]=2.*rand()/RAND_MAX-1;
		sumvx=sumvx+vx[i]; sumvy=sumvy+vy[i];
		sumv2=sumv2+(vx[i]*vx[i]+vy[i]*vy[i])*0.50;
	}

	sumv2=sumv2/(double)npart; sumvx=sumvx/(double)npart; sumvy=sumvy/(double)npart;
	fs=pow(temp/sumv2,0.5);

	// centro de massa estático (termostato alpha1)
	for(int i=0;i<npart;i++){
		vx[i]=(vx[i]-sumvx)*fs; vx[i]=(vy[i]-sumvy)*fs;
	}
}

// força de interação entre as partículas
double force(double *fx, double *fy, double *x, double *y, double Ep){
	double ene=0.0;
	for(int i=0;i<npart;i++){
	fx[i]=fy[i]=0.0;
	}
	for(int i=0;i<npart-1;i++){
	for(int j=i+1;j<npart;j++){
	    double xr,yr,r,rinv,ffx,ffy;
	    xr=x[i]-x[j]; yr=y[i]-y[j];
	    if (xr>=0.5*caixa){ xr=xr-caixa; }
	    if (xr<=-0.5*caixa){ xr=xr+caixa; }
	    if (yr>=0.5*caixa){ yr=yr-caixa; }
	    if (yr<=-0.5*caixa){ yr=yr+caixa; }
	    r=pow((xr*xr+yr*yr),0.5);
	    if(r<=rcri){
		    rinv=1.0/r;
		    ffx=4.0*xr*(alfa)*(pow(rinv,2.0))*(2.0*(pow(rinv,(2.0*alfa)))-(pow(rinv,alfa)));
		    ffy=4.0*yr*(alfa)*(pow(rinv,2.0))*(2.0*(pow(rinv,(2.0*alfa)))-(pow(rinv,alfa)));
		    ene=ene+4.0*(pow(rinv,(2.0*alfa))-pow(rinv,(alfa)));
		    fx[i]=fx[i]+ffx; fy[i]=fy[i]+ffy; fx[j]=fx[j]-ffx; fy[j]=fy[j]-ffy;
	    }
	}
	}

	Ep=ene/(double)npart;
	return Ep;
}

// integração das posições
void integrate(double *fx, double *fy, double *x, double *y, double *vx, double *vy, double *vvx, double *vvy){
	for(int i=0;i<npart;i++){
		vvx[i]=vx[i]+0.50*dt*fx[i]; vvy[i]=vy[i]+0.50*dt*fy[i];
		x[i]=x[i]+dt*vvx[i]; y[i]=y[i]+dt*vvy[i];
	}
}

// integração das velocidades
double integraVEL(double *fx, double *fy,double *vx, double *vy, double *vvx, double *vvy, double Ec){
	double Ecin=0;
	for(int i=0;i<npart;i++){
		vx[i]=vvx[i]+0.50*dt*fx[i]; vy[i]=vvy[i]+0.50*dt*fy[i];
		Ecin=Ecin+(vx[i]*vx[i]+vy[i]*vy[i])*0.50;
	}

	Ec = Ecin/(double)npart;
	return Ec;
}

// termostato
double termostato(double *vx, double *vy, double Ec){
	double tal=dt,vfac,Ecin;
	Ecin= Ec;
	vfac=pow(1.0+(dt/tal)*((temp/Ecin)-1.0),0.5);
	Ecin=0.0;
	for(int i=0;i<npart;i++){
		vx[i]=vx[i]*vfac; vy[i]=vy[i]*vfac;
		Ecin=Ecin+(vx[i]*vx[i]+vy[i]*vy[i])*0.5;
	}

	Ec=Ecin/(double)npart;
	return Ec;
}

// condição periódica de contorno
void CPC(double *x,double *y){
	for(int i=0;i<npart;i++){
		if (x[i]>=0.5*caixa){ x[i]=x[i]-caixa; }
		if (x[i]<=-0.5*caixa){ x[i]=x[i]+caixa; }
		if (y[i]>=0.5*caixa){ y[i]=y[i]-caixa; }
		if (y[i]<=-0.5*caixa){ y[i]=y[i]+caixa; }
	}
}

