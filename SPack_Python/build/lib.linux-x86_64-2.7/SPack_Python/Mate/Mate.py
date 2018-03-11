import numpy as np
import SPack_Light as SPL

def Between(x,Xmin,Xmax):
	if x >= Xmin and x < Xmax:	return True
	return False

def Between2(x,Xmin,Xmax):
	if x >= Xmin and x <= Xmax:	return True
	return False

def Between_Arr(X,Xmin,Xmax, take_min=True, take_max = False):
	if take_min: Cond1 = (X >= Xmin)
	else: Cond2 = (X > Xmin)
	if take_max: Cond2 = (X <= Xmax)
	else: Cond2 = (X < Xmax)
	return Cond1 & Cond2

def Line(x,X1,X2,Y1,Y2):
	m = (Y2-Y1)/(X2-X1)
	return Y1+m*(x-X1)

def is_Sort(A):
	for i in range(1,len(A)):
		if A[i-1] > A[i]:	return False
	return True

def Bars(Xmin,Xmax,NBin,X,Y,bar = 0.2):
	Axis = (np.array(range(NBin))+0.5)*(Xmax-Xmin)/float(NBin)+Xmin
	Data = np.zeros((NBin,3))
	ix = ((X-Xmin)/float(Xmax-Xmin)*NBin).astype(int)
	for i in range(NBin):
		tmp_id = np.where(i == ix)[0]
		tmp_Y = np.sort(Y[tmp_id])
		Len = len(tmp_Y)
		if Len > 1:
			Data[i] = np.nanmedian(tmp_Y),tmp_Y[int(Len*bar)],tmp_Y[int(Len*(1-bar))]
		elif Len == 1:
			Data[i] = np.nanmedian(tmp_Y),np.nan,np.nan
		else:
			Data[i] = np.nan,np.nan,np.nan
	return Axis, Data

def Histo(X,Xmin,Xmax,NBin):
	X_Arr = (np.array(range(NBin))+0.5)*(Xmax-Xmin)/NBin+Xmin
	Arr = np.zeros(NBin,dtype=int)
	SPL.Histo(X[~np.isnan(X)].astype(float),Arr,float(Xmin),float(Xmax),int(NBin))
	return X_Arr,Arr
