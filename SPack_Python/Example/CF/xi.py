#Main Libs
import numpy as np

import SPack_Python.RW as RW
import SPack_Python.Advance as SPPA


NBin  = 40
Xmin  = -2.5
Xmax  = 1.5
Lbox  = 62.5
iNBin = 10
LIMIT = 30
JN3	= 3

x,y,z = RW.Read_Advance('data.txt', Nprop = 3, Sep = ',')

####### READ #######

x_axis,Xi,Xi_err = SPPA.Xi(x,y,z,NBin,Xmin,Xmax,Lbox,iNBin,LIMIT,JN3)
 
RW.Print_Basic3([x_axis,Xi,Xi_err],'xi_1.txt')
