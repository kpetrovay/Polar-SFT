#!/usr/bin/env python 
# program parammap.py ; invoke with arguments flowtype tau [year]
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')     # don't require X11, so the code can be run with nohup
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf
import scipy.interpolate
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

flowtype=1
tau0=100.0
flowtype=int(sys.argv[1])
tau0=float(sys.argv[2])

vegtelen=float("inf")

# Observational constraints:
sigmafac=2

meanrevWSOB=4.33 
stdrevWSOB=0.36
meanrevdipmom=3.44 
stdrevdipmom=0.18
meanrelWSOB=0.94 
stdrelWSOB=0.07
meanreldipmom=0.80 
stdreldipmom=0.11
meanhmcycmin=70.0
stdhmcycmin=2.5

#infilename='profile' + str(sys.argv[1])+ '.txt'
infilename='params' + str(sys.argv[1])+ '.dat'
#outfilename_revWSO='maps_case' + str(sys.argv[1]) + '/revWSO_' + str(sys.argv[2]) +'.png'
#outfilename_revdipmom='maps_case' + str(sys.argv[1]) + '/revdipmom_' + str(sys.argv[2]) +'.png'
outfilename='maps_case' + str(sys.argv[1]) + '/' + 'maps_case' + str(sys.argv[1]) + '_' + str(sys.argv[2]) +'.png'

f=open(infilename,"r")
lines=f.readlines()
tau=[]
u0=[]
eta=[]
revBpol=[]
revWSOB=[]
revdipmom=[]
relWSOB=[]
reldipmom=[]
hmcycmin=[]
hmpb=[]
hmWS=[]
for x in lines:
	if (float(x.split('	')[2]))==tau0:
		u0.append(float(x.split('	')[0]))
		eta.append(float(x.split('	')[1]))
		tau.append(float(x.split('	')[2]))
		revBpol.append(float(x.split('	')[3]))
		revWSOB.append(float(x.split('	')[4]))
		revdipmom.append(float(x.split('	')[5]))
		relWSOB.append(float(x.split('	')[6]))
		reldipmom.append(float(x.split('	')[7]))
		hmcycmin.append(float(x.split('	')[8]))
		hmpb.append(float(x.split('	')[9]))
		hmWS.append(float(x.split('	')[10]))
f.close()

x = np.array(u0)
y = np.array(eta)

plt.ioff()        # set on/off interactive mode, so fig.is redrawn every draw() command or only at the end.
fig = plt.figure(1,figsize=(12,12))
fig.suptitle('Flow ' + str(flowtype)  + ', $\\tau=$' + str(sys.argv[2]))
fig.subplots_adjust(wspace=0.3)

plt.subplot(321)

loWSOB=meanrevWSOB-sigmafac*stdrevWSOB
hiWSOB=meanrevWSOB+sigmafac*stdrevWSOB

#fig=plt.figure()
#ax=fig.add_subplot(111)#, projection='3d')
xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)
zi = scipy.interpolate.griddata((x, y), revWSOB, (xi, yi), method='cubic')
plt.imshow(zi, vmin=-11.0, vmax=11.0,origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto",cmap='hsv')
plt.colorbar().set_label('WSO polar field reversal time')
E=plt.contour(xi, yi,zi,levels = [loWSOB,hiWSOB],colors=('k',),linestyles=('-',),linewidths=(0.5))
plt.clabel(E, fmt = '%.2f', colors = 'k', fontsize=14)
excl1a=plt.contourf(xi, yi,zi,levels = [-vegtelen,loWSOB],colors=('gray',),alpha=0.5)
excl1b=plt.contourf(xi, yi,zi,levels = [hiWSOB,vegtelen],colors=('gray',),alpha=0.5)
plt.xlabel('$u_0$')
plt.ylabel('$\eta$')

#plt.show()
#raw_input("Press [enter] to continue.")

#fig.savefig(outfilename_revWSO, format='png', dpi=300)

#sys.exit("Bye!")

plt.subplot(325)

lodipmom=meanrevdipmom-sigmafac*stdrevdipmom
hidipmom=meanrevdipmom+sigmafac*stdrevdipmom

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)
zi = scipy.interpolate.griddata((x, y), revdipmom, (xi, yi), method='cubic')
plt.imshow(zi, vmin=-11.0, vmax=11.0,origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto",cmap='hsv')
plt.colorbar().set_label('Dipolar moment reversal time')
E=plt.contour(xi, yi,zi,levels = [lodipmom,hidipmom],colors=('k',),linestyles=('-',),linewidths=(0.5))
plt.clabel(E, fmt = '%.2f', colors = 'k', fontsize=14)
plt.contourf(xi, yi,zi,levels = [-vegtelen,lodipmom],colors=('gray',),alpha=0.5)
plt.contourf(xi, yi,zi,levels = [hidipmom,vegtelen],colors=('gray',),alpha=0.5)
plt.xlabel('$u_0$')
plt.ylabel('$\eta$')

plt.subplot(322)

loWSOB2=meanrelWSOB-sigmafac*stdrelWSOB
hiWSOB2=meanrelWSOB+sigmafac*stdrelWSOB

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)
zi = scipy.interpolate.griddata((x, y), relWSOB, (xi, yi), method='cubic')
plt.imshow(zi, vmin=-1, vmax=1.0,origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto",cmap='YlGn')
plt.colorbar().set_label('WSO polar field at min./max.value')
E=plt.contour(xi, yi,zi,levels = [loWSOB2,hiWSOB2],colors=('k',),linestyles=('-',),linewidths=(0.5))
plt.clabel(E, fmt = '%.2f', colors = 'k', fontsize=14)
plt.contourf(xi, yi,zi,levels = [-vegtelen,loWSOB2],colors=('gray',),alpha=0.5)
plt.contourf(xi, yi,zi,levels = [hiWSOB2,vegtelen],colors=('gray',),alpha=0.5)
plt.xlabel('$u_0$')
plt.ylabel('$\eta$')


plt.subplot(326)

lodipmom2=meanreldipmom-sigmafac*stdreldipmom
hidipmom2=meanreldipmom+sigmafac*stdreldipmom

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)
zi = scipy.interpolate.griddata((x, y), reldipmom, (xi, yi), method='cubic')
plt.imshow(zi, vmin=-1, vmax=1.0,origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto",cmap='YlGn')
plt.colorbar().set_label('Dipolar moment at min./max value')
E=plt.contour(xi, yi,zi,levels = [lodipmom2,hidipmom2],colors=('k',),linestyles=('-',),linewidths=(0.5))
plt.clabel(E, fmt = '%.2f', colors = 'k', fontsize=14)
plt.contourf(xi, yi,zi,levels = [-vegtelen,lodipmom2],colors=('gray',),alpha=0.5)
plt.contourf(xi, yi,zi,levels = [hidipmom2,vegtelen],colors=('gray',),alpha=0.5)
plt.xlabel('$u_0$')
plt.ylabel('$\eta$')


plt.subplot(323)

lohmcycmin=meanhmcycmin-sigmafac*stdhmcycmin
hihmcycmin=meanhmcycmin+sigmafac*stdhmcycmin

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)
zi = scipy.interpolate.griddata((x, y), hmcycmin, (xi, yi), method='cubic')
plt.imshow(zi, vmin=45, vmax=90,origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto",cmap='copper')
plt.colorbar().set_label('Latitude of edge of polar topknot at cycle min.')
E=plt.contour(xi, yi,zi,levels = [lohmcycmin,hihmcycmin],colors=('k',),linestyles=('-',),linewidths=(0.5))
plt.clabel(E, fmt = '%.2f', colors = 'k', fontsize=14)
plt.contourf(xi, yi,zi,levels = [-vegtelen,lohmcycmin],colors=('gray',),alpha=0.5)
plt.contourf(xi, yi,zi,levels = [hihmcycmin,vegtelen],colors=('gray',),alpha=0.5)
plt.xlabel('$u_0$')
plt.ylabel('$\eta$')


ax=plt.subplot(324)

xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
xi, yi = np.meshgrid(xi, yi)
#zi = np.zeros((100,100))
#plt.imshow(zi, vmin=0, vmax=100,origin='lower', extent=[x.min(), x.max(), y.min(), y.max()], aspect="auto",cmap='Greys')
zi = scipy.interpolate.griddata((x, y), revWSOB, (xi, yi), method='cubic')
plt.contourf(xi, yi,zi,levels = [loWSOB,hiWSOB],colors=('gray',),alpha=0.5)
#zi = scipy.interpolate.griddata((x, y), revdipmom, (xi, yi), method='cubic')
#plt.contourf(xi, yi,zi,levels = [lodipmom,hidipmom],colors=('gray',),alpha=0.5)
zi = scipy.interpolate.griddata((x, y), relWSOB, (xi, yi), method='cubic')
plt.contourf(xi, yi,zi,levels = [loWSOB2,hiWSOB2],colors=('gray',),alpha=0.5)
#zi = scipy.interpolate.griddata((x, y), reldipmom, (xi, yi), method='cubic')
#plt.contourf(xi, yi,zi,levels = [lodipmom2,hidipmom2],colors=('gray',),alpha=0.5)
zi = scipy.interpolate.griddata((x, y), hmcycmin, (xi, yi), method='cubic')
plt.contourf(xi, yi,zi,levels = [lohmcycmin,hihmcycmin],colors=('gray',),alpha=0.5)
plt.xlabel('$u_0$')
plt.ylabel('$\eta$')














#plt.show()
#raw_input("Press [enter] to continue.")

fig.savefig(outfilename, format='png', dpi=200)



#sys.exit("Bye!")



