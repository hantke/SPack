import numpy as np
import SPack_Light as SL

def Line(x,X1,X2,Y1,Y2):
    return (Y2-Y1)/(X2-X1)*(x-X1)+Y1

def FuncSingle(x0,X,Y):
	N = len(X)
	if x0 < X[0]:	return Line(x0,X[0],X[1],Y[0],Y[1])
	for i in range(len(X)):	
		if (x0 < X[i]):	return Line(x0,X[i-1],X[i],Y[i-1],Y[i])
	return Line(x0,X[N-2],X[N-1],Y[N-2],Y[N-1])

N = 100000
X = np.linspace(0., 1., num=N)
Y = np.linspace(0., 10., num=N)

for i in range(N):
	#FuncSingle(np.random.random(),0.,1.,0.,10.)
	#SL.Line(1.,0.,10.,0.,30.)
	FuncSingle(np.random.random(),X,Y)
	#rand = np.random.random()
	#SL.FuncSingle(rand,X,Y)
	#print rand,SL.FuncSingle(rand,X,Y)
#print SL.FuncSingle(0.5,X,Y)

#Rand = np.random.random(N)

#A = SL.FuncArr(Rand,X,Y)

#for i in range(len(A)):
	#print Rand[i],A[i]
