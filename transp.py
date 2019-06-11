#!/usr/bin/env python 
# program transp.py ; invoke with arguments u0 eta tau flowtype in m/s km^2/s year;
# flowtype may take 1 2 3 4 or 5 see below
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as special

# Units:  Megameter ; day ; Gauss

cycleper = 11.0 * 365.25            # cycle period in days
flowtype = int(sys.argv[4])

# grid resolution 1/Nfac degree; time step 1/sl day
# CFL crit -> sl ~ Nfac^2
if flowtype ==1:  
    Nfac=2
    sl=4
if flowtype ==2:  
    Nfac=1
    sl=1
if flowtype ==3:
    Nfac=2
    sl=4
#Nfac=2
N=180*Nfac+1                        # no. of grid points  DO NOT CHANGE!
#sl=4                                # slowness factor; time step = 1 day /sl
# sl=0.5 OK for merid.flows (1) and (2) with N=181 and eta = 300
dt=1.0/sl
nt=int(20*cycleper/dt)              # run for this many cycles
nout = int(0.5*cycleper/dt)         # plot at every nout-th time step
wrinterv = 365.25/12.0              # save output at these intervals
tres=0.0                            # to reset time for last 2 cycles
tres0=tres
#nt=80000                           # time steps to run
#nout=188                           # plot at every nout-th time step 

Rsun = 695.7                        # solar radius 

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

################ Define meridional flow #####################


# (1) Dikpati et al. 2006; used by Cameron et al. 2007:
if flowtype==1:
    latitude0 = 90           # trick to improve flux conservation near poles
    uc = u0*np.sin(np.pi*latitude/latitude0)
    uc[np.where(abs(latitude) > latitude0)]=0

# (2) van Ballegooijen 1998; used by Jiang et al. 2014:
# Jiang 2014:  u0=11  eta=250  tau=infty
if flowtype==2:
    latitude0 = 75.0
    uc = u0*np.sin(np.pi*latitude/latitude0)
    uc[np.where(abs(latitude) > latitude0)]=0

# (3) Lemerle et al. 2017
# u0=18  eta=600  tau=10
if flowtype==3:
    latitude0 = 89  	# trick to improve flux conservation near poles
    V=7.0
    W=1.0
    q=1
    uc=u0*(special.erf(V*np.cos(np.pi/2*latitude/latitude0)))**q * special.erf(W*np.sin(np.pi/2*latitude/latitude0))
    uc[np.where(abs(latitude) > latitude0)]=0

# (4) Whitbread et al. 2017
# 1D: u0=14   eta=350	tau=2.4  p=3.24  or  tau=1.9  (if p=1.87 fixed)
# 2D: u0=11   eta=450	tau=5	    p=2.15  or  p=2.76  (for no tau)  
if flowtype==4:
    p=3.24
    uc =u0*((np.sin(theta))**p)*np.cos(theta)

# (5) Wang 2017
# eta=500    no tau
if flowtype==5:
    u0=13
    uc =u0*np.tanh(np.pi/2*latitude/6)*(np.cos(np.pi/2*latitude))**2

outfilename = 'res_case' + str(flowtype) + '/tr_' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3]) +'.idt'
#outfilename = 'tr_' + str(sys.argv[1]) + '_' + str(sys.argv[2]) + '_'  + str(sys.argv[3]) +'.idt'
#print "Output will be written to  ", outfilename

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
    if flowtype ==1:  
        sourcescale = 0.003   # merid.flow (1): 0.003; (2): 0.015; (3): 0.0005
    if flowtype ==2:  
        sourcescale = 0.015   # merid.flow (1): 0.003; (2): 0.015; (3): 0.0005 
    if flowtype ==3:  
        sourcescale = 0.0005   # merid.flow (1): 0.003; (2): 0.015; (3): 0.0005
    ampli = sourcescale * ahat * tc**3 / (np.exp(tc**2/bhat**2) - chat)
    
    # Latitudinal profile as in Cameron et al. 2007:
    cycleno = int(t // cycleper) + 1
    evenodd = 1 - 2*(cycleno%2)             # 1 for even cycle, -1 for odd
    # Latitudinal profile as in Cameron et al. 2007:
    #lambda0 = 35.0 - 30.0*((t/cycleper)%1)  # (t%cycleper) 
    # Latitudinal profile as in Jiang et al.2011:
    lambda0 = 26.4 - 34.2*((t/cycleper)%1) + 16.1*((t/cycleper)%1)**2  # Jiang et al.2011
    fwhm = 6.0
    #bandn1 = evenodd * ampli * np.exp(-(latitude-lambda0-0.5)**2/2/fwhm**2)
    #bandn2 = -evenodd * ampli * np.exp(-(latitude-lambda0+0.5)**2/2/fwhm**2)
    #bands1 = evenodd * ampli * np.exp(-(latitude+lambda0-0.5)**2/2/fwhm**2)
    #bands2 = -evenodd * ampli * np.exp(-(latitude+lambda0+0.5)**2/2/fwhm**2)
    joynorm = 0.5/np.sin(20.0/180*np.pi)
    bandn1 = evenodd * ampli *np.exp(-(latitude-lambda0-joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    bandn2a = -evenodd * ampli *np.exp(-(latitude-lambda0+joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    bands2a = evenodd * ampli *np.exp(-(latitude+lambda0-joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    bands1 = -evenodd * ampli *np.exp(-(latitude+lambda0+joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    
    ### calculate flux correction due to spherical geom., so net flux is zero:
    ### use a higher resolution here for accuracy
    Nfine=180+1
    thetaf=np.linspace( 0, np.pi, Nfine )    # x = theta: colatitude in radian
    latitudef=90.0-thetaf*180/np.pi	     # latitude in degrees
    bandn1f = evenodd * ampli *np.exp(-(latitudef-lambda0-joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    bandn2af = -evenodd * ampli *np.exp(-(latitudef-lambda0+joynorm*np.sin(lambda0/180*np.pi))**2/2/fwhm**2)
    integrand1=-np.sin(thetaf)*bandn1f
    fluxband1=np.trapz(integrand1, thetaf)
    integrand2=-np.sin(thetaf)*bandn2af
    fluxband2=np.trapz(integrand2, thetaf)
    #print fluxband1, fluxband2
    fluxcorrection=1.0
    if (ampli != 0): fluxcorrection=-fluxband1/fluxband2
    #fluxcorrection=1
    ### end fluxcorrection
    
    # correct amplitude of one ring to ensure zero net flux
    bandn2=fluxcorrection*bandn2a
    bands2=fluxcorrection*bands2a
    
    value=bandn1+bandn2+bands1+bands2
    #print tc, ampli, value[60]
    return value

################ Initialization ##############################

t=0.0
B0=0.001
B=B0*np.sin(np.pi*latitude/180)
B=B0*np.zeros(N)
#B=B0*np.ones(N)
W = Rsun * np.sin(theta) * B   # W: annular flux density

S=source(latitude,t)
#plt.ion()	 # set interactive mode, so fig.is redrawn every draw() command.
#fig = plt.figure()
#plt.title('1D flux transport')
#plt.xlabel('latitude')
#plt.ylabel('$B$')
#plt.plot(latitude, B);
#fig.canvas.draw()
#raw_input("Press [enter] to continue.")

################ Upwind integration (FTBS) ###################

def uproll(blk):
    value = np.copy(blk)
    for j in range(N):
        k=int(hemi[j])
        tmpblk = np.roll(blk,k)
        value[j] = tmpblk[j]
    return value

###############################################################################

#print u0, eta, tau

with open(outfilename, 'w') as ofile:
    for n in range(nt):  #loop for values of n from 0 to nt, so it will run nt times
        t=n*dt
        # Compute divergence dWflux of Wflux, the "flux of flux":
        # Wflux =  Wfladv + Wfldif
        Wfladv = W/Rsun
        Wfladv*=uc
        # use upwind differencing for advective term:
        #Wflupw = uproll(Wfladv)
        #dWflux = hemi*(Wfladv-Wflupw)
        # or use centered differencing for advective term:
        dWflux = (np.roll(Wfladv,-1)-np.roll(Wfladv,1))/2
        #dWflux = np.zeros(N)		       # to switch off advection
        # use centered differencing for diffusive term:
        # Wfldif =  eta/Rsun * np.sin(x) * dB/d(x)
        Wfldifr = np.sin(x+dx/2) * (np.roll(B,-1)-B)/dx
        Wfldifl = np.sin(x-dx/2) * (B-np.roll(B,1))/dx
        dWflux+= eta/Rsun * (Wfldifr-Wfldifl)	 # # to switch off diffusion
        S=source(latitude,t)
        # step W and re-evaluate B:
        mgflux=np.trapz(B[0:(N-1)/2+1]*np.sin(theta[0:(N-1)/2+1]), dx=np.pi/(N-1))/2
        dW = dWflux/dx - W/tau  + S * Rsun * np.sin(theta)
        dW*=dt
        W+=dW
        B[1:N-1] = W[1:N-1]/Rsun/np.sin(x[1:N-1])
        # assuming 3rd derivative =0 :
        B[0]=B[2]+0.5*(B[1]-B[3])
        B[N-1]=B[N-3]+0.5*(B[N-2]-B[N-4])
        # assuming 3rd derivative =const. :
        #B[0]=4*B[1]-6*B[2]+4*B[3]-B[4]   
        #B[N-1]=4*B[N-2]-6*B[N-3]+4*B[N-4]-B[N-5]  
        if (n*dt/cycleper>=nt*dt/cycleper-2):# in last 2 cycles
            if ( int(tres/wrinterv) > (tres0/wrinterv) ):
                ofile.write(str(tres/365.25)+'\n')
                Bsamp =  B[0::Nfac]                       # sample B at 181 pts 
                Bsamplist = np.ndarray.tolist(Bsamp)      # convert to list
                for item in Bsamplist:
                    ofile.write("%s\n" % item )
                    #plt.plot(latitude, B,label='t= %.2f yr' %(t/365.25))
                    #plt.legend()
                    #plt.legend(bbox_to_anchor=(0.98, 1),loc=2,borderaxespad=0.)
                    #fig.canvas.draw()
                    twrite=0.0
            tres0=tres
            tres+=dt
ofile.closed

print wrinterv

#raw_input("Press [enter] to terminate.")
