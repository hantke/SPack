#Main Libs
import numpy as np

import SPack_Python.RW as RW
import SPack

def GR(NBin = 50,Xmin = -3,Xmax = 2,Lbox = 62.5,N_Prim=100,N_Sec=100):
	gr = []
	for i in range(NBin):
		Min = 10**(Xmin+i/float(NBin)*(Xmax-Xmin))
		Max = 10**(Xmin+(i+1)/float(NBin)*(Xmax-Xmin))
		Vol = 4./3.*3.1415926*((Max**3)-(Min**3))
		Den = N_Sec*Lbox**(-3)
		gr.append(Vol*Den*N_Prim)
	return np.array(gr)

######################
NBin  = 40
Xmin  = -2.5
Xmax  = 1.5
Lbox  = 62.5
iNBin = 10
LIMIT = 30
JN	= 27
x_axis = (np.array(range(NBin))+0.5)*(Xmax-Xmin)/NBin+Xmin
######################

####### READ #######
x,y,z = RW.Read_Advance('data.txt', Nprop = 3, Sep = ',')
x = x.astype(float)
y = y.astype(float)
z = z.astype(float)
####### READ #######

gg = np.zeros(JN*NBin,dtype=float)
JN_Random = ((x/Lbox*3).astype(int) + (y/Lbox*3).astype(int)*3 + (z/Lbox*3).astype(int)*9) #JN RANDOM SAMPLE

SPack.ACF_DD(x,y,z,JN_Random,gg,NBin,iNBin,LIMIT,Xmin,Xmax,Lbox,JN,1,1)

## JN ERROR CALCULATION

gg_Arr  = gg.reshape(JN, NBin)
gg_Full = np.sum(gg.reshape(JN, NBin),axis=0) 
gg_JN   = gg_Full-gg_Arr

Len = len(JN_Random)
gr_Full = GR(NBin,Xmin,Xmax,Lbox,N_Prim=Len,N_Sec=Len)
gr_JN   = np.zeros((JN,NBin))

for i in range(JN):	
	gr_JN[i] =  GR(NBin,Xmin,Xmax,Lbox,N_Prim=(Len-len(np.where(JN_Random == i)[0])),N_Sec=Len)

## JN ERROR CALCULATION
	
Xi = gg_Full/gr_Full-1

Xi_err = np.sqrt(np.sum((Xi - (gg_JN/gr_JN-1))**2,axis=0)*(JN-1)/JN)

RW.Print_Basic3([x_axis,Xi,Xi_err],'xi_1.txt')
