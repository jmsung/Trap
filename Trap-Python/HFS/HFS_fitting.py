# Kim I. Mortensen
# November 20, 2015

# Import various modules
from pylab import *
import numpy as np
from scipy.special import erfc, i1
from scipy.integrate import romberg
from scipy.optimize import leastsq, fmin

import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 20
mpl.rcParams['ytick.labelsize'] = 20

#--------------------------  USER INPUT  ------------------------------

# Sampling frequency
fsample=20000.0 # [Hz]

# neg/pos force. True means flip Kim's convention
F0_sign = True

# Define molecules
mol1={}
mol1['directory']='C:/Users/chaoliu/Documents/Spudich Lab/Data/20151218_WT808_oscillation/cal1_WT808_m2_HFS'
mol1['fmin']=-6.0 # Lower cut-off on the force axis [pN]
mol1['nfits']=5 # Number of bins on the force axis 
mol1['df']=1.0 # Width of bins on the force axis [pN]
mol1['fdrive']=100.0 # Oscillation frequency [Hz].
#mol1['tlower']=5.0 # Lower cut-off for binding times [ms]. This should be equal to one oscillation period.
mol1['tlower'] = 1000.0/mol1['fdrive']
mol1['tupper']=200.0 # Upper cut-off for binding times [ms]. 
mol1['scut']=5.0 # Cut-off for variances of binding times [ms^2].
mol1['chi2cut']=50.0 # Upper cut-off for average chi-square of binding event (optional).

#mol2={}
#mol2['directory']='Molecule_2'
#mol2['fmin']=
#mol2['nfits']=
#mol2['df']=
#mol2['tlower']=
#mol2['tupper']=
#mol2['fdrive']=
#mol2['scut']=
#mol2['chi2cut']=

# List of molecules to be analyzed
molecules=[mol1] # [mol1, mol2, ...] 

# Modify the bin size and the lower time cut-off
m=1 # The bin size on the time axis is m*tdrive. m is a positive integer.
tlfact=1 # The lower cutoff on the time axis is tdrive*tlfact. tlfact is a positive integer. 

#-------------------------- END USER INPUT --------------------------

# Thermal energy
kBT=4.11 # (pN nm)

# Define functions

def discreteexp(n,a,nmin,nmax,m=1):
    # The probability of unbinding in the m periods following n periods of attachment
    # Note that n takes the values [nmin,nmin+m,nmin+2m,...,nmin+zm]
    z=(nmax-nmin)/m # number of bins -1
    r=exp(-a)
    return r**n*(1-r**m)/(r**nmin*(1-r**((z+1)*m)))


class LogLikelihood:
    """ Class defining the log-likelihood function maximized in MLE."""

    def __init__(self,counts,ns):
        # Counts in the fitting bins
        self.counts=counts
        # Bin number according to theory
        self.ns=ns
        # Number of counts in fitting interval
        self.N=sum(counts)
        print "Total counts",self.N
        
    def Value(self,a):

        ns=self.ns
        nmin=ns[0]
        nmax=ns[-1]

        m=ns[1]-ns[0]

        # Probability in fitting interval
        # by summing over geometric series
        r=exp(-a)

        # Expected counts in each fitting bin
        mu=discreteexp(ns,a,nmin,nmax,m)*self.N

        value=sum(-mu+self.counts*log(mu))
        value*=-1.0
        
        return value

for mol in molecules:

    directory=mol['directory']
    f_min=mol['fmin']
    nfits=mol['nfits']
    df=mol['df']
    tlower=mol['tlower']*tlfact
    tupper=mol['tupper']
    fdrive=mol['fdrive']
    scut=mol['scut']
    chi2cut=mol['chi2cut']

    bindingtimes=array(np.loadtxt(directory+'/'+'bindingtimes_Wavelet.txt'))
    F0=array(np.loadtxt(directory+'/'+'F0_Wavelet.txt'))
    if F0_sign:
        F0 = -F0
    dF=array(np.loadtxt(directory+'/'+'dF_Wavelet.txt'))
    avchi2=array(np.loadtxt(directory+'/'+'avchi2.txt'))

    # convert bindingtimes and period to ms
    bindingtimes/=(fsample/1000.0)
    tdrive=1000.0/fdrive

    # leave outliers out
    sumvar=var(bindingtimes[:,[0,2]],axis=1,ddof=1)+var(bindingtimes[:,[1,3]],axis=1,ddof=1)

    select1=sumvar<scut
    t12=((bindingtimes[:,1]-bindingtimes[:,0])+(bindingtimes[:,3]-bindingtimes[:,2]))/2.0
    select2=t12>(1000.0/fdrive)
    select3=avchi2[:,0]<chi2cut
    select4=avchi2[:,1]<chi2cut
    select=select1*select2*select3*select4
    print "Number of outliers=", shape(F0)[0]-sum(select)

    nevents=arange(shape(dF)[0])
    nevents=nevents[select]
    bindingtimes=bindingtimes[select,:]
    F0=F0[select,:]
    dF=dF[select,:]

    T=fsample/fdrive # units of data points

    # Do some plotting

    figure(1,figsize=(12,8))
    subplot(221)
    plot(F0[:,0],bindingtimes[:,1]-bindingtimes[:,0],'b.')
    plot(F0[:,1],bindingtimes[:,3]-bindingtimes[:,2],'r.')
    ylabel('Binding time (ms)',size=16)
    xlabel('Mean force (pN)',size=16)
    ylim(0,200)

    subplot(222)
    plot(sum(F0,axis=1),((bindingtimes[:,1]-bindingtimes[:,0])+(bindingtimes[:,3]-bindingtimes[:,2]))/2.0,'k.')
    ylabel('Binding time (ms)',size=16)
    xlabel('Mean force (pN)',size=16)
    ylim(0,200)

    subplot(223)
    plot(F0[:,0],abs(dF[:,0]),'b.')
    plot(F0[:,1],abs(dF[:,1]),'r.')
    xlabel('Mean force (pN)',size=16)
    ylabel('Force amplitude (pN)',size=16)

    subplot(224)
    plot(sum(F0,axis=1),abs(dF[:,0])+abs(dF[:,1]),'k.')

    xlabel('Mean force (pN)',size=16)
    ylabel('Force amplitude (pN)',size=16)


    # Summarize forces and binding times
    offset_f=sum(F0,axis=1)
    t12_f=((bindingtimes[:,1]-bindingtimes[:,0])+(bindingtimes[:,3]-bindingtimes[:,2]))/2.0
    deltaf=abs(dF[:,0])+abs(dF[:,1])

    # Eliminate data points outside the chosen time period
    select=(t12_f>tlower)*(t12_f<tupper)
    offset_f=offset_f[select]
    t12_f=t12_f[select]
    deltaf=deltaf[select]
    nevents=nevents[select]

    print "Number of data points=", len(offset_f)

    # Arrays for storage of results
    fs=array([])
    fserror=array([])
    tinvs=array([])
    tinverrors=array([])

    n_included=0

    # Loop over bins on the force axis
    for i in range(nfits):

        # Select data points with mean force in the current bin
        select1=offset_f > f_min+df*i
        select2=offset_f < f_min+df*(i+1)
        select=select1*select2

        ts_temp=t12_f[select]
        deltaf_temp=deltaf[select]
        offset_f_temp=offset_f[select]
        fs_temp=mean(offset_f_temp)

        print "\nForce range =", f_min+df*i,' - ', f_min+df*(i+1)

        ts=ts_temp
        fs=append(fs,mean(offset_f_temp))
        fserror=append(fserror,sqrt(var(offset_f_temp)/len(offset_f_temp)))

        figure()
        title('Range='+str(f_min+df*i)+','+str(f_min+df*(i+1)))
        xlabel('Binding time (ms)',size=16)
        ylabel('Number of events', size=16)
                
        # Define bins on time axis
        lims=arange(tlower,tupper+m*tdrive,m*tdrive)
        nmin=int(tlower/tdrive)
        z=len(lims)-2
        nmax=nmin+z*m
        ns=arange(nmin,nmax+m,m)
        
        # Histogram of binding times
        h=hist(ts,bins=lims,facecolor=(40.0/256,117.0/256,250.0/256))
        d=h[1][1]-h[1][0]
        
        # Number of data points in the current histogram
        counts=h[0]
        n_included+=sum(counts)

        # Perform maximum likelihood estimation of the (inverse) ADP release rate
        mlinit=[0.1]
        ll=LogLikelihood(counts,ns)
        mlfit=fmin(ll.Value,mlinit,maxiter=500,ftol=0.00001)
        a=mlfit 

        # ADP release rate [s^-1]
        tinvs=append(tinvs,1000*a/tdrive)

        # Plot discrete exponential distribution
        r=exp(-a)
        y=sum(counts)*discreteexp(ns,a,nmin,nmax,m)
        trange=(array(ns)+m*0.5)*tdrive
        plot(trange,y,'k-',lw=2.0)
        
        # Calculate errors of the release rate estimates
        n=array(ns)
        mu=discreteexp(n,a,nmin,nmax,m)*sum(counts)
        dmua=sum(counts)*((-n*r**n+(n+m)*r**(n+m))*(r**nmin-r**((z+1)*m+nmin))-(r**n-r**(n+m))*(-nmin*r**nmin+((z+1)*m+nmin)*r**((z+1)*m+nmin)))\
          /(r**nmin-r**((z+1)*m+nmin))**2

        fisherinfo=sum((dmua)**2/mu)
        tinvvar=(1000.0/tdrive)**2/fisherinfo
        tinverrors=append(tinverrors,sqrt(tinvvar))
        

    print "\nNumber of included data points=", n_included

    # Calculate average of deltaf
    deltafav=mean(deltaf)

    # Define load-dependent release rate function
    def kappa(f0,k0,delta):
        val=k0*i0(deltafav*delta/kBT)*exp(-f0*delta/kBT)
        return val

    # Perform maximum likelihood estimation of k0 and delta
    parinit=array([70.0,-1.0])
    parfit=leastsq(lambda par: (tinvs-kappa(fs,par[0],par[1]))/tinverrors,parinit)[0]
    k0=parfit[0]
    delta=parfit[1]    

    # Fisher's information matrix for k0 and delta
    Ikk=sum(1/tinverrors**2*(i0(deltafav*delta/kBT)*exp(-fs*delta/kBT))**2)
    Idd=sum(1/tinverrors**2*k0**2*exp(-2*fs*delta/kBT)*(deltafav/kBT*i1(deltafav*delta/kBT)-\
                                                 fs/kBT*i0(deltafav*delta/kBT))**2)
    Ikd=sum(1/tinverrors**2*(i0(deltafav*delta/kBT)*exp(-2*fs*delta/kBT))*k0\
            *(deltafav/kBT*i1(deltafav*delta/kBT)-fs/kBT*i0(deltafav*delta/kBT)))

    I=zeros((2,2))
    I[0,0]=Ikk
    I[0,1]=I[1,0]=Ikd
    I[1,1]=Idd
    covar=inv(I)

    # Get the variance elements
    k0var=covar[0,0]
    deltavar=covar[1,1]

    # Print the results
    print "\nk0 [1/s] =", k0, '+/-', sqrt(k0var)
    print "delta [nm] =", delta, '+/-', sqrt(deltavar)

    # Plot results
    figure()
    yfit=kappa(fs,parfit[0],parfit[1])
    errorbar(fs,tinvs,yerr=tinverrors,linestyle='none',lw=2.0,marker='o',ms=7.0,color='k')
    fsplot=arange(min(fs),max(fs),0.01)
    yfit=kappa(fsplot,parfit[0],parfit[1])
    plot(fsplot,yfit,'r-',lw=2.0)
    xlabel('Mean force (pN)',size=16)
    ylabel('ADP release rate (s$^{-1}$)',size=16)

    print "\nk(0,delta_F) [1/s] =", kappa(0.0,parfit[0],parfit[1])

show()

