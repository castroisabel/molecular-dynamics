# bibliotecas
import random
import math
import numpy as np
import matplotlib.pyplot as plt
from array import *
import time as tm

antes = tm.time()

# condições de contorno
npart=50; densi=0.2; caixa=(float(npart)/densi)**0.5; temp=0.5; ene=0.0; dt=0.01;
alfa=6.0; kapa=2.0; Ao=0.0; rcri=caixa/2; Ecin=0.0;
tmax=200000;Ecin=0;

# arquivos
posicaoi=open("posi.dat","w") # posições iniciais
blender=open("blender.xyz","w") # blender
posicaof=open("posf.dat","w") # posições finais
energia=open("Energia.dat","w") # energia potencial e mecânica


# vetores para as partículas
x=array('f',range(0,npart)); y=array('f',range(0,npart));
vx=array('f',range(0,npart)); vy=array('f',range(0,npart));
fx=array('f',range(0,npart)); fy=array('f',range(0,npart));
vvx=array('f',range(0,npart)); vvy=array('f',range(0,npart));


# posições iniciais rede quadrada        
def posicaoinitial():
	a=(caixa)/((pow(float(npart),0.5)+1)); ax=a;ay=a;
	for i in range(0,npart,1):
		x[i]=(-caixa*0.50)+ax; y[i]=(-caixa*0.50)+ay;
		ax=ax+a
		if ax>=caixa:
		    ax=a; ay=ay+a
	return

    
# velocidades iniciais aleatórias
def velocidadeinicial():
	sumvx=sumvy=sumv2=0.0
	for i in range(0,npart,1):
		vx[i]=random.uniform(-1, 1)
		vy[i]=random.uniform(-1, 1)
		sumvx=sumvx+vx[i]; sumvy=sumvy+vy[i];
		sumv2=sumv2+(vx[i]*vx[i]+vy[i]*vy[i])*0.50
	    
	sumv2=sumv2/float(npart); sumvx=sumvx/float(npart); sumvy=sumvy/float(npart);
	fs=(temp/sumv2)**0.50
	# termostato alpha1
	for i in range(0,npart,1):                               
		vx[i]=(vx[i]-sumvx)*fs; vx[i]=vay=(vy[i]-sumvy)*fs
	      
	return
       
    
# forças de interação  
def force():
	ene=0
	for i in range(0,npart,1):
		fx[i]=fy[i]=0.0
	for i in range(0,npart-1,1):       
		for j in range(i+1,npart,1):
			xr=x[i]-x[j]; yr=y[i]-y[j];
			if (xr>=0.50*caixa):
				xr=xr-caixa
			if (xr<=-0.50*caixa):
				xr=xr+caixa
			if (yr>=0.50*caixa):
				yr=yr-caixa
			if (yr<=-0.50*caixa):
				yr=yr+caixa
			r=(xr*xr+yr*yr)**0.5
			if (r<=rcri):
				rinv=1.0/r;
				ffx=4.0*xr*(alfa)*(pow(rinv,2.0))*(2.0*(pow(rinv,(2.0*alfa)))-(pow(rinv,alfa)))
				ffy=4.0*yr*(alfa)*(pow(rinv,2.0))*(2.0*(pow(rinv,(2.0*alfa)))-(pow(rinv,alfa)))
				ene=ene+4.0*(pow(rinv,(2.0*alfa))-pow(rinv,(alfa)))
				fx[i]=fx[i]+ffx; fy[i]=fy[i]+ffy; fx[j]=fx[j]-ffx; fy[j]=fy[j]-ffy;

	ene=ene/float(npart)
	return ene


# integração das posições
def integrate():
	for i in range(0,npart,1):
		vvx[i]=vx[i]+0.50*dt*fx[i]; vvy[i]=vy[i]+0.50*dt*fy[i];
		x[i]=x[i]+dt*vvx[i]; y[i]=y[i]+dt*vvy[i];
	return


# integração das velocidades
def integraVEL():
	Ecin=0
	for i in range(0,npart,1):
		vx[i]=vvx[i]+0.50*dt*fx[i]; vy[i]=vvy[i]+0.50*dt*fy[i]; 
		Ecin=Ecin+(vx[i]*vx[i]+vy[i]*vy[i])*0.50;
		
	Ecin=Ecin/float(npart)
	return Ecin
	
	
# termostato
def termostato(Ecin):
	tal=dt
	vfac=(1.0+(dt/tal)*((temp/Ecin)-1.0))**0.5
	Ecin=0.0;
	for i in range(0,npart,1):
	      vx[i]=vx[i]*vfac;vy[i]=vy[i]*vfac;
	      Ecin=Ecin+(vx[i]*vx[i]+vy[i]*vy[i])*0.50;
	      
	Ecin=Ecin/float(npart)
	return Ecin


# condição periódica de contorno
def CPC():
	for i in range(0,npart,1):
		if (x[i]>=0.5*caixa):
			x[i]=x[i]-caixa
		if (x[i]<=-0.5*caixa):
			x[i]=x[i]+caixa      
		if (y[i]>=0.5*caixa):
			y[i]=y[i]-caixa
		if (y[i]<=-0.5*caixa):
			y[i]=y[i]+caixa

	return

P = []  # Energia Potencial
T = []	# Energia Cinética
M = []	# Energia Mecânica



# inicío
posicaoinitial() 
velocidadeinicial()

# plt.scatter(x,y)
# plt.axis([-caixa*0.5, caixa*0.5, -caixa*0.5,caixa*0.5])

# arquivo posi.dat
blender.write(str(npart) + "\n")
blender.write("Centro da Partícula \n")
for i in range(0,npart,1):
	posicaoi.write(str(x[i]) + "      " + str(y[i]) + "\n")
	blender.write("Default  " + str(x[i]) + "       " + str(y[i]) + "                0 \n")
       
force()

# variáveis auxiliares
tobs = 100 # número de frames
tmod=tmax/tobs
time = np.linspace(1, tmax, tobs)*dt # vetor tempo


for t in range(1,tmax+1,1):
	integrate()
	CPC()
	Ep=force()
	Ecin=integraVEL()
	Ecin=termostato(Ecin)
	Emec=Ep+Ecin
	energia.write(str(t) + "      " + str(Ep) + "      " + str(Emec) + "\n")
	if(t%tmod==0):
		P.append(Ep)
		T.append(Ecin)
		M.append(Emec)
		blender.write(str(npart) + "\n")
		blender.write("Centro da Partícula \n")
		for i in range(0,npart,1):
			blender.write("Default  " + str(x[i]) + "       " + str(y[i]) + "                0 \n")
#		plt.scatter(x,y)
#		plt.pause(0.05)
	print(t,Ep,Ecin,Emec)


# arquivo posf.dat
for i in range(0,npart,1):
	posicaof.write(str(x[i]) + "      " + str(y[i]) + "\n")



plt.plot(time,P,'-r',time,T,'-b',time,M,'-m')
plt.xlabel('tempo')
plt.ylabel('Energias')
plt.legend(['Potencial', 'Cinética','Mecânica'])
plt.savefig('gráfico.png')

depois = tm.time()
intervalo = depois - antes
print(f'Tempo: {intervalo} segundos')
