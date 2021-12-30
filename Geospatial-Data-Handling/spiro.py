from math import cos,sin,pi
import numpy as np


R = 8
r = 1
a = 4
ti,xc,yc = 0,0,0
nRev = 16
tt = [34.0205777,-118.2854465]
f= open("spirograph.txt","w+")
	
for ti in np.arange(0,(pi*nRev),0.01):
		
	xc = tt[0] + ((R+r)*cos((r/R)*ti) - a*cos((1+r/R)*ti))/10000
	yc = tt[1] + ((R+r)*sin((r/R)*ti) - a*sin((1+r/R)*ti))/10000
	f.write(str(yc) + "," + str(xc) + "\n")
		
f.close()

