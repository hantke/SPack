from matplotlib.pyplot import rc
import matplotlib

def Load_RC():
	try:
		rc('axes', labelpad = 10, labelsize=7, linewidth=2)
		rc('xtick', labelsize=30,bottom=True,top=True,direction='in') 
		rc('ytick', labelsize=30,left=True,right=True,direction='in') 
		rc('xtick.minor',visible=True,size=6, width=1, pad = 4) # 1_dist
		rc('ytick.minor',visible=True,size=6, width=1, pad = 4) #
		rc('xtick.major',             size=12,width=2, pad = 3) # All_Dist
		rc('ytick.major',             size=12,width=2, pad = 3) #
		rc('figure.subplot',left = 0.05,bottom = 0.05,right = 0.95, top = 0.95 )
		rc('savefig',dpi=150)
		rc('figure',dpi=80)
		rc('path',simplify=True)
		rc('figure',figsize= "10, 10")
		rc('legend',numpoints=1,frameon=False,handletextpad=0.3)
		#rc('font',family='serif')
		matplotlib.rcParams['pdf.fonttype'] = 42
	except:
		print 'Using old Matplotlib version...'
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

    #"    rc('xtick', labelsize=30,bottom=True,top=True,direction='in') \n",
    #"    rc('ytick', labelsize=30,left=True,right=True,direction='in') \n",
    #"    rc('xtick.minor',visible=True,size=6, width=1, pad = 4) # 1_dist\n",
    #"    rc('ytick.minor',visible=True,size=6, width=1, pad = 4) #\n",
    #"    rc('xtick.major',             size=12,width=2, pad = 3) # All_Dist\n",
    #"    rc('ytick.major',             size=12,width=2, pad = 3) #\n",
    #"    rc('figure.subplot',left = 0.05,bottom = 0.05,right = 0.95, top = 0.95 )\n",
    #"    rc('savefig',dpi=150,extension='pdf')\n",
    #"    rc('path',simplify=True)\n",
    #"    rc('font',family='serif')\n",
    #"    rc('mathtext',fontset='custom')\n",
    #"    rc('legend',numpoints=1,frameon=False,handletextpad=0.3)\n",
    #"    matplotlib.rcParams[\"pdf.fonttype\"] = 3\n",
    #"    #ax.xaxis.labelpad = 20 #Dist_2\n",
    #"    #ax.tick_params(axis='x', pad=30) #Dist_1"
