#Main Libs
import numpy as np
import SPack

def Find_ainA(a,A,a0,a1):#Look for a in Y-Table, initial a0,a1 = 0, len(A-1)
	if a1 - a0 < 3:
		t = a0
		while t <= a1 and t < len(A):
			if a == A[t]:	return t
			t += 1
		return -99	
	med = (a1+a0)/2
	if a < A[med]:	return Find_ainA(a,A,a0,med)
	if a > A[med]:	return Find_ainA(a,A,med,a1)
	return med	

def Find_AinB(A, B):# Only if A is complete in B
	b=max(0,Find_ainA(A[0],B,0,len(B)-1)-1)
	a = 0
	#b=0
	index = []
	while a<len(A):
		if (A[a] < B[b]) or (len(B)<b):
			print 'warning', a, b, A[a], B[b]
			return -1
		#if b%1000000==0: print a, b, A[a], B[b]
		if A[a]== B[b]:
			index.append(b)
			a+=1
		b+=1
	return np.array(index)

def Find_A1A2inB(A1, A2, B):# Only if A1[a] <= B[b] <= A2[a] NOT TESTED
	a = 0
	b = 0
	index = np.zeros(len(A1),dtype=int)-99
	while a<len(A1) and b <len(B):
				#if b%1000000==0: print a, b, A[a], B[b]
		if B[b] > A2[a]:
			a += 1
			b -= 1
		elif A2[a] >= B[b] and A1[a] <= B[b]:
			index[a] = b
			a+=1
		b+=1
	return index

def GR(NBin = 50,Xmin = -3,Xmax = 2,Lbox = 62.5,N_Prim=100,N_Sec=100):
	gr = []
	for i in range(NBin):
		Min = 10**(Xmin+i/float(NBin)*(Xmax-Xmin))
		Max = 10**(Xmin+(i+1)/float(NBin)*(Xmax-Xmin))
		Vol = 4./3.*3.1415926*((Max**3)-(Min**3))
		Den = N_Sec*Lbox**(-3)
		gr.append(Vol*Den*N_Prim)
	return np.array(gr)

def Xi(x, y, z, NBin  = 40, Xmin  = -2.5, Xmax  = 1.5, Lbox  = 62.5, iNBin = 10, LIMIT = 30, JN3	= 3):
	x_axis = (np.array(range(NBin))+0.5)*(Xmax-Xmin)/NBin+Xmin
	
	x = x.astype(float)
	y = y.astype(float)
	z = z.astype(float)
	JN = JN3**3
	
	gg = np.zeros(JN*NBin,dtype=float)
	JN_Random = ((x/Lbox*JN3).astype(int) + (y/Lbox*JN3).astype(int)*JN3 + (z/Lbox*JN3).astype(int)*JN3**2) #JN RANDOM SAMPLE
	
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
		
	Xi = gg_Full/gr_Full*2-1
	
	Xi_err = np.sqrt(np.sum((Xi - (gg_JN/gr_JN*2-1))**2,axis=0)*(JN-1)/JN)
	
	return np.array([x_axis,Xi,Xi_err])

def Shuffle(Pos,HaloID,HaloMass,Type,HaloProp = [],GalProp = [],Lbox = 500,Mmin = 8,Mmax = 17 ,MNBin = 90):
	
	##
	if any((HaloMass < Mmin) | (HaloMass >= Mmax)):
		print '\n\n\n #################################' 
		print      '\n            WANING\n\n'
		print 'HaloMass MUST be between Mmin and Mmax'	
		print  min(HaloMass),max(HaloMass),np.nanmedian(HaloMass)
		print '\n\n\n #################################' 
		return 0
	##
	if len(HaloProp) == 0: use_HaloProp = False
	else: use_HaloProp = True
	if len(GalProp) == 0: use_GalProp = False
	else: use_GalProp = True
	#SORTING
	ind = np.lexsort((Type,HaloID))
	Pos = Pos[ind]
	HaloID = HaloID[ind]
	Type = Type[ind]
	HaloMass = HaloMass[ind]
	if use_HaloProp: HaloProp = HaloProp[ind]
	if use_GalProp: GalProp = GalProp[ind]

	Len = len(Pos)
	
	Centrals = np.where(Type == 0)[0]
	NGr             = len(Centrals)
	Gr_Mass         = HaloMass[Centrals]
	Gr_Pos          = Pos[Centrals]
	Gr_ID           = HaloID[Centrals]
	Gr_ID2CentalGal = np.copy(Centrals)
	if use_HaloProp: Gr_Prop         = HaloProp[Centrals]

	lastHaloID = HaloID[0]
	for i in range(1,Len):
		if Type[i] == 0: lastHaloID = HaloID[i]
		elif HaloID[i] != lastHaloID: print 'WARNING, satelites without centrals ', i,HaloID[i]

	for i in range(Len):
		if Type[i] == 0: Pos_Central = np.copy(Pos[i])
		Pos[i] -= Pos_Central
	for i in range(len(Pos[0])):    
		Cond1 = np.where(Pos[:,i] > Lbox/2.)[0]
		Cond2 = np.where(Pos[:,i] < -Lbox/2.)[0]
		Pos[:,i][Cond1]  = Pos[:,i][Cond1] - Lbox
		Pos[:,i][Cond2]  = Pos[:,i][Cond2] + Lbox
		
	iGr_Mass = ((Gr_Mass-Mmin)/(Mmax-Mmin)*MNBin).astype(int)
	
	for i in range(MNBin):
		tmpID = np.where(iGr_Mass == i)[0]
		tmpID2 = np.copy(tmpID)
		np.random.shuffle(tmpID2)
		Gr_ID2CentalGal[tmpID] = Gr_ID2CentalGal[tmpID2]
		
	for i in range(NGr):
		Pos[Gr_ID2CentalGal[i]] = Gr_Pos[i]
		HaloID[Gr_ID2CentalGal[i]] = Gr_ID[i]
		HaloMass[Gr_ID2CentalGal[i]] = Gr_Mass[i]
		if use_HaloProp: HaloProp[Gr_ID2CentalGal[i]] = Gr_Prop[i]

	for i in range(Len):
		if Type[i] == 0: CentralID = i
		else:
			Pos[i] += Pos[CentralID]
			HaloID[i] += HaloID[CentralID]
			HaloMass[i] += HaloMass[CentralID]
			if use_HaloProp: HaloProp[i] += HaloProp[CentralID]
	for i in range(len(Pos[0])):    
		Cond1 = np.where(Pos[:,i] >= Lbox)[0]
		Cond2 = np.where(Pos[:,i] < 0)[0]
		Pos[:,i][Cond1]  = Pos[:,i][Cond1] - Lbox
		Pos[:,i][Cond2]  = Pos[:,i][Cond2] + Lbox
		
	if use_HaloProp and use_GalProp:
		return Pos,HaloID,HaloMass,Type,HaloProp,GalProp
	elif use_HaloProp:
		return Pos,HaloID,HaloMass,Type,HaloProp
	elif use_GalProp:
		return Pos,HaloID,HaloMass,Type,GalProp
	else:
		return Pos,HaloID,HaloMass,Type
	
