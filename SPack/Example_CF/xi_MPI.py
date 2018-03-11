#Main Libs
import numpy as np

import SPack_Python.RW as RW
import SPack

from mpi4py import MPI
from time import sleep
#mpirun -c 4 python xi_MPI.py

##
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
##

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
if rank == 0:
	Pos = RW.Read_Advance('data.txt', Nprop = 3, Sep = ',')
else: Pos = None
Pos = comm.bcast(Pos,root=0)
x,y,z = Pos
x = x.astype(float)
y = y.astype(float)
z = z.astype(float)

####### READ #######

gg     = np.zeros(JN*NBin,dtype=float)
gg_cpu = np.zeros(JN*NBin,dtype=float)

comm.Barrier()
if rank == 0:
	JN_Random = ((x/Lbox*3).astype(int) + (y/Lbox*3).astype(int)*3 + (z/Lbox*3).astype(int)*9) #JN RANDOM SAMPLE
else: JN_Random = None
JN_Random = comm.bcast(JN_Random,root=0)

SPack.ACF_DD(x,y,z,JN_Random,gg_cpu,NBin,iNBin,LIMIT,Xmin,Xmax,Lbox,JN,size,rank+1)
comm.Barrier()
comm.Reduce(gg_cpu, gg, MPI.SUM,0)


if rank == 0:
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
	
	RW.Print_Basic3([x_axis,Xi,Xi_err],'xi_2.txt')
