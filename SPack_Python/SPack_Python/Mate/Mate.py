def Between(x,Xmin,Xmax):
	if x >= Xmin and x < Xmax:	return True
	return False

def Between2(x,Xmin,Xmax):
	if x >= Xmin and x <= Xmax:	return True
	return False

def Line(x,X1,X2,Y1,Y2):
	m = (Y2-Y1)/(X2-X1)
	return Y1+m*(x-X1)

def is_Sort(A):
	for i in range(1,len(A)):
		if A[i-1] > A[i]:	return False
	return True
