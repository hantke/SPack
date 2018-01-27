import numpy as np
################################################################

def Print_Basic(X,Y,FileName):
	File= open("Out/"+FileName,'w')
	for i in range(len(X)):
		print >>File,X[i],Y[i]
	File.close()
	print 'Write '+FileName
	
def Big_String(A,i):
	Str = ''
	for a in A:	Str += str(a[i])+'	'
	return Str
def Print_Basic2(A,FileName):
	File= open(FileName,'w')
	for i in range(len(A[0])):	print >>File,Big_String(A,i)
	File.close()
	print 'Write '+FileName
def Print_Basic3(A,FileName):
	File= open(FileName,'w')
	A = np.array(A)
	for i in range(len(A[0])):
		if all(A[:,i] == A[:,i]): print >>File,Big_String(A,i)
		else: print >>File,'#'+Big_String(A,i)
	File.close()
	print 'Write '+FileName

def Read_Basic(FileName):
	F = open(FileName, "r")
	X = []
	while True:
		linea = F.readline()
		L = linea.split()
		if not linea: break
		if L[0] != '#' and list(L[0])[0] != '#':
			X.append(float(L[0]))
	
	F.close()
	print 'Read ',FileName
	return X

def Read_Advance(FileName, Nprop = 1, SpecialPos = None, Arr = True, T = True,Sep = None):
	F = open(FileName, "r")
	X = []
	if SpecialPos == None:	SpecialPos = range(Nprop)
	if len(SpecialPos) != Nprop:	print 'FAIL! len(SpecialPos) != Nprop',len(SpecialPos),Nprop
	while True:
		linea = F.readline()
		L = linea.split(Sep)
		if not linea: break
		if L[0] != '#' and list(L[0])[0] != '#':
			Xtmp = []
			for i in SpecialPos:	Xtmp.append(float(L[i]))
			X.append(np.array(Xtmp))
	F.close()
	if Arr:	X = np.array(X)
	if T:	X = X.T
	
	print 'Read ',FileName
	return X
	