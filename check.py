#!/usr/bin/env python 
# program check.py ; invoke with arguments u0 eta tau flow
# in units m/s km^2/s year nondim
import sys
import numpy as np
#import scipy.integrate as integrate
#import scipy as sp
import matplotlib
matplotlib.use('Agg')     # don't require X11, so the code can be run with nohup
import matplotlib.pyplot as plt
#import sympy as sm
import scipy.special as special

# Units:  Megameter ; day ; Gauss

cycleper = 11.0                     # cycle period in years
Rsun = 695.7                        # solar radius 

N=181                               # no. of grid points

theta=np.linspace( 0, np.pi, N )    # x = theta: colatitude in radian
latitude=90.0-theta*180/np.pi       # latitude in degrees
hemi=np.sign(latitude)              # +1 in N, -1 in S, 0 on eq.

x=theta
dx=np.pi/(N-1)

#merid.flow amplitude m/s:
u0 = float(sys.argv[1])
#u0=10.0
u0*=8.64e-2

#diffusivity in km**2/s:
eta = float(sys.argv[2])
#eta = 300.0         # Cameron et al. 2007, reproducing Dikpati et al. 2006
eta*=8.64e-2        # eta ~ 20-30 in our units

#tau in years:
tau = float(sys.argv[3])
#tau = 10.0          # Cameron et al. 2007; term accounting for radial diff.
tau*=365.25

#flowtype=1
flowtype=int(sys.argv[4])

infilename = 'res_case'+ str(flowtype) + '/tr_' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3]) +'.idt'
#infilename = 'tr_' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3]) +'.idt'
plotfilename = 'plots_case' + str(flowtype) + '/plots_' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3]) +'.png'
bflyfilename = 'plots_case' + str(flowtype) + '/bfly_' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3]) +'.png'
datafilename='params' + str(flowtype) + '.dat'
#print "Output will be written to  ", datafilename, plotfilename

################ Define source ##############################

def source(latitude,t):
    global ampli

    # Time profile of source: set constant factor arbitrarily so B_max ~ 15 G
    # (a) simple sin^2: 
    #ampli = 10.0 * (np.sin(((t/cycleper)%1)*np.pi))**2  
    # (b) Hathaway et al. (1994):
    tc = 12.0 * ((((t/cycleper)%1)*cycleper/365.25)) 
    ahat = 0.00185   # "average" cycle, ~ Cycle 17
    bhat = 48.7  
    chat = 0.71  
    sourcescale = 0.003   # merid.flow (1): 0.003; (2): 0.015; (3): 0.0005
    ampli = sourcescale * ahat * tc**3 / (np.exp(tc**2/bhat**2) - chat)
    
    # Latitudinal profile as in Cameron et al. 2007:
    cycleno = int(t // cycleper) + 1
    evenodd = 1 - 2*(cycleno%2)             # 1 for even cycle, -1 for odd
    # Latitudinal profile as in Cameron et al. 2007:
    lambda0 = 35.0 - 30.0*((t/cycleper)%1)  # (t%cycleper) 
    # Latitudinal profile as in Jiang et al.2011:
    #lambda0 = 26.4 - 34.2*((t/cycleper)%1) + 16.1*((t/cycleper)%1)**2  # Jiang et al.2011
    fwhm = 6.0
    #bandn1 = evenodd * ampli * np.exp(-(latitude-lambda0-0.5)**2/2/fwhm**2)
    #bandn2 = -evenodd * ampli * np.exp(-(latitude-lambda0+0.5)**2/2/fwhm**2)
    #bands1 = evenodd * ampli * np.exp(-(latitude+lambda0-0.5)**2/2/fwhm**2)
    #bands2 = -evenodd * ampli * np.exp(-(latitude+lambda0+0.5)**2/2/fwhm**2)
    joynorm = 0.5/np.sin(20.0/180*np.pi)
    bandn1 = evenodd * ampli *np.exp(-(latitude-lambda0-joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    bandn2 = -evenodd * ampli *np.exp(-(latitude-lambda0+joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    bands1 = evenodd * ampli *np.exp(-(latitude+lambda0-joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    bands2 = -evenodd * ampli *np.exp(-(latitude+lambda0+joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    value=bandn1+bandn2+bands1+bands2
    #print tc, ampli, value[60]
    return value

###############################################################################

#print u0, eta, tau

# read file into a list, stripping the \n newline chars:
with open(infilename, 'r') as ifile:
    ifilelist = ifile.read().splitlines()
ifile.closed

# convert list of strings to array of floats:
ifilearray = np.array([float(item) for item in ifilelist])

# number of time steps in file:
steps = len(ifilearray)/182

# rearrange array to (182 x steps) matrix of t and B:
tbmx = np.reshape(ifilearray, (steps,182))

t = tbmx[:,0]

NpolarB = tbmx[:,1]


tc = 12.0 * ((((t/cycleper)%1)*cycleper)) 
ahat = 0.00185   # "average" cycle, ~ Cycle 17
bhat = 48.7  
chat = 0.71  
sourcescale = 0.003   # merid.flow (1): 0.003; (2): 0.015; (3): 0.0005
ampli = 30*sourcescale * ahat * tc**3 / (np.exp(tc**2/bhat**2) - chat)

#raw_input("Press [enter] to continue.")

# WSO observed polar field observed with finite resolution:
# and dipole moment + mgraph saturation (factor 1.8):
thetadegmax = 35
thetavar = theta[0 : thetadegmax+1]
WSOB = np.zeros(steps)
dipmom = np.zeros(steps)
for j in range(0,(steps)):
    Bvec = tbmx[j,1:N+1]
    integrand = (Bvec*(np.sin(theta))**2)[0 : thetadegmax+1]
    WSOB[j] = 5.0*np.trapz(integrand, thetavar)/1.8
    integrand2 = Bvec*np.sin(2*theta)
    dipmom[j] = 0.75*np.trapz(integrand2, theta)

relWSOB=np.abs(WSOB[0])/np.max(WSOB)
reldipmom=np.abs(dipmom[0])/np.max(dipmom)

plt.ioff()        # set on/off interactive mode, so fig.is redrawn every draw() command or only at the end.
fig = plt.figure(1,figsize=(12,4))
fig.suptitle('Flow ' + str(flowtype) + ', $u_0=$' + str(sys.argv[1]) + ', $\eta=$' + str(sys.argv[2]) + ', $\\tau=$' + str(sys.argv[3]))
fig.subplots_adjust(hspace=-0.1)
plt.subplot(121)
plt.title('Polar field variation')
plt.xlabel('time [years]')
plt.ylabel('$B$')
Brenormf=int(max(dipmom))
plt.plot(t, NpolarB/Brenormf,color='k', label='$B/$'+str(Brenormf));
plt.plot(t, WSOB,color='b', label='WSO polar field');
plt.plot(t, dipmom,color='r', label='dipolar moment');
plt.plot(t,ampli/Brenormf,color='grey', label='SSN (rescaled)')
plt.plot(t,np.zeros(steps))
plt.legend(prop={'size': 8},loc=3,)
fig.canvas.draw()


# generic routines to find zero crossings:
def zcrall(x, y):
    origlen = len(x)
    xtrunc = x[0:origlen-1]
    #for all zero crossings:
    return xtrunc[np.diff(np.sign(y)) != 0]
    #for negative-to-positive zero crossings:
    #return xtrunc[np.diff(np.sign(y)) > 0]
def zcr(x, y):
    origlen = len(x)
    xtrunc = x[0:origlen-1]
    #for all zero crossings:
    #return xtrunc[np.diff(np.sign(y)) != 0]
    #for negative-to-positive zero crossings:
    return xtrunc[np.diff(np.sign(y)) > 0]
# now look for reversal times based on this:
def reversaltime(x, y):
    origlen = len(x)
    xtrunc = x[0:origlen-1]
    #for all zero crossings:
    #return xtrunc[np.diff(np.sign(y)) != 0]
    #for negative-to-positive zero crossings:
    reversals=xtrunc[np.diff(np.sign(y)) > 0]
    reversal=reversals[0]
    reversalno=len(reversals)
    # now take care of triple reversals within 7 years:
    for i in range(reversalno-1):
        thisrev=reversals[i]
        nextrev=reversals[i+1]
        if ( (nextrev-thisrev) < 7 ): reversal=0.5*(nextrev+thisrev)
    return reversal
revBpol = reversaltime(t,NpolarB)-11
revWSOB = reversaltime(t,WSOB)-11
revdipmom = reversaltime(t,dipmom)-11

# look for maxima (minima):
phaseBpol = t[np.argmin(NpolarB)]
phaseWSO = t[np.argmin(WSOB)]
phasedipmom = t[np.argmin(dipmom)]


# plot field profile at max/minimum:
plt.subplot(122)
plt.title('Magnetic field profile at various instants')
plt.xlabel('Latitude')
plt.ylabel('Magnetic field')
plt.subplots_adjust(hspace=0.4)

Bprof = tbmx[0,1:N+1]
# determine latitude at half max.:
halfmaxcycmin = np.abs(zcrall(-latitude,(Bprof-Bprof[0]/2))[0])
plt.plot(latitude, -Bprof,color='gray', label='cycle min.');

Bprof = tbmx[np.argmin(NpolarB),1:N+1]
# determine latitude at half max.:
halfmaxnpolarb = np.abs(zcrall(-latitude,(Bprof-Bprof[0]/2))[0])
plt.plot(latitude, Bprof,color='k',label='max. polar field');

Bprof = tbmx[np.argmin(WSOB),1:N+1]
# determine latitude at half max.:
halfmaxWSOB = np.abs(zcrall(-latitude,(Bprof-Bprof[0]/2))[0])
plt.plot(latitude, Bprof,color='b',label='max. WSO field');

Bprof = tbmx[np.argmin(dipmom),1:N+1]
# determine latitude at half max.:
halfmaxdipmom = np.abs(zcrall(latitude,(Bprof-Bprof[0]/2))[0])
plt.plot(latitude, Bprof,color='r', label='max.dip.mom.');

plt.legend(prop={'size': 8},loc=9,)

plt.savefig(plotfilename)

#print phaseBpol, phaseWSO, phasedipmom, 

#print sys.argv[1], sys.argv[2], sys.argv[3], revBpol, revWSOB, revdipmom, halfmaxnpolarb, halfmaxWSOB
printout =  sys.argv[1] + '\t' + sys.argv[2] + '\t' + sys.argv[3] + '\t' + str(revBpol) + '\t' + str(revWSOB) + '\t' + str(revdipmom) + '\t'  + str(relWSOB) + '\t' + str(reldipmom) + '\t' + str(halfmaxcycmin) + '\t' + str(halfmaxnpolarb) + '\t' + str(halfmaxWSOB) + '\n'
with open(datafilename, 'a') as outtable:
    outtable.write(printout)
outtable.closed


# Mohamed's addition to plot butterfly dgr:
tbmxB=tbmx[:,1:]
MM = zip(*tbmxB)

x = np.array(t)
y = np.array(latitude)
fig=plt.figure()
ax=fig.add_subplot(111)#, projection='3d')
vmaxvalue=min(np.max(tbmxB),20)
plt.imshow(MM,extent=[np.min(x), np.max(x), np.min(y), np.max(y)] , aspect='auto', interpolation='none',cmap='gray', vmin=-vmaxvalue, vmax=vmaxvalue)
plt.xlabel('time [years]')
plt.ylabel('latitude')
fig.suptitle('Flow ' + str(flowtype) + ', $u_0=$' + str(sys.argv[1]) + ', $\eta=$' + str(sys.argv[2]) + ', $\\tau=$' + str(sys.argv[3]))
#plt.title('Polar field variation')
plt.colorbar().set_label('field strength [G]')
ziz=np.transpose(tbmxB)
E=plt.contour(x, y,ziz,levels = [0.0],colors=('k',),linestyles=('--',),linewidths=(0.5))
plt.clabel(E, fmt = '%.2f', colors = 'k', fontsize=14)
#plt.show()
# Mohamed's addition ends here

#raw_input("Press [enter] to continue.")

plt.savefig(bflyfilename, dpi=200)


#raw_input("Press [enter] to terminate.")
