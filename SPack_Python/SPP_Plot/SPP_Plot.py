from matplotlib.pyplot import rc
import matplotlib

def Load_RC():
	rc('axes', labelpad = 10, labelsize=7, linewidth=2)
	rc('xtick', labelsize=30,direction='in') 
	rc('ytick', labelsize=30,direction='in') 
	rc('xtick.minor',visible=True,size=6, width=1, pad = 4) # 1_dist
	rc('ytick.minor',visible=True,size=6, width=1, pad = 4) #
	rc('xtick.major',             size=12,width=2, pad = 3) # All_Dist
	rc('ytick.major',             size=12,width=2, pad = 3) #
	rc('figure.subplot',left = 0.05,bottom = 0.05,right = 0.95, top = 0.95 )
	rc('savefig',dpi=150,extension='pdf')
	rc('path',simplify=True)
	rc('figure',figsize= "10, 10")
	#rc('font',family='serif')
	#rc('mathtext',fontset='custom')
	rc('legend',numpoints=1,frameon=False,handletextpad=0.3)
	matplotlib.rcParams["pdf.fonttype"] = 3
	#ax.xaxis.labelpad = 20 #Dist_2
	#ax.tick_params(axis='x', pad=30) #Dist_1
	matplotlib.rcParams['pdf.fonttype'] = 42

