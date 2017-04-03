# Kim I. Mortensen
# September 15, 2015

# Import various modules
import numpy as np
import os
from pylab import *
from scipy import convolve
from scipy.optimize import leastsq
from operator import itemgetter
from itertools import groupby

#--------------------------  USER INPUT  ------------------------------

# Sampling frequency
fsample=20000.0 # [Hz]

# Molecule 1
mol1={}
mol1['directory']='C:/Users/chaoliu/Documents/Spudich Lab/Data/20151218_WT808_oscillation/cal1_WT808_m2_HFS'
mol1['datadir']='cal1_WT808_m2'
mol1['cal_b1']=172.0 # Calibration factor [nm/V]
mol1['cal_b2']=106.0 # Calibration factor [nm/V]
mol1['cal_k1']=0.114 # Trap strength [pN/nm]
mol1['cal_k2']=0.07 # Trap strength [pN/nm]
mol1['fd']=100.0 # Oscillation frequency [Hz]
mol1['phasecut']=0.5 # Threshold for the phase [rad]
mol1['ampthreshold1']=9.0 # Threshold for the amplitude [nm]
mol1['ampthreshold2']=9.0 # Threshold for the amplitude [nm]

## Molecule 2
#mol2={}
#mol2['directory']='Molecule_2'
#mol2['datadir']='Cal1_M2'
#mol2['cal_b1']=
#mol2['cal_b2']=
#mol2['cal_k1']=
#mol2['cal_k2']=
#mol2['fd']=
#mol2['phasecut']=
#mol2['ampthreshold1']=
#mol2['ampthreshold2']=

# List of molecules to be analyzed
molecules=[mol1] # [mol1, mol2, ...]

#-------------------------- END USER INPUT --------------------------


# Auxiliary function
def mergeevents(events):
    merged=True
    while any(merged):
        mergedevents=[]
        merged=repeat(False,len(events)-1)
        n=0
        while n<len(events):
            if n<len(events)-1 and (events[n+1][0]-events[n][1] < 500):
                mergedevents.append((events[n][0],events[n+1][1]))
                merged[n]=True
                n+=2
            else:
                mergedevents.append((events[n][0],events[n][1]))
                n+=1
        events=mergedevents
    return events

# Close all figures at startup
close('all')

nevents_all=0

for mol in molecules:

    # Get all parameters for the relevant molecule
    directory=mol['directory']
    datadir=mol['datadir']
    cal_b1=mol['cal_b1']
    cal_b2=mol['cal_b2']
    cal_k1=mol['cal_k1']
    cal_k2=mol['cal_k2']
    fd=mol['fd']
    phase_threshold=mol['phasecut']
    amp_threshold1=mol['ampthreshold1']
    amp_threshold2=mol['ampthreshold2']

    print directory

    # Make results directory if it does not exist
    resultsdir=directory+'/Results_WaveletAnalysis' #_test8
    if not os.path.exists(resultsdir):
        os.mkdir(resultsdir)

    # Make list of all real data arrays in data directory
    datafiles=os.listdir(directory+'/'+datadir)
    if '.DS_Store' in datafiles:
        datafiles=datafiles[1:]
    if 'processed.txt' in datafiles:
        datafiles.pop()

    # Sort datafiles according to numeric value
    datafiles.sort(key=lambda s: int(s.rsplit('.')[0]))

    # Plot 
    plotcolors=['black','gray']
    highlightcolors=['yellow','cyan']

    # Period of stage osciallation
    T=fsample/fd # [data points]
    
    # List to store binding times, row is event number and columns contain at1,at2,det1,det2
    bindingtimes_all=[]
    F0_all=[]
    dF_all=[]

    bindingtimes_all_Fourier=[]
    F0_all_Fourier=[]
    dF_all_Fourier=[]
    avchi2_all=[]

    nevent_save=0

    # Loop over datafiles
    for datafile in datafiles:
        print 'Analyzing '+datafile
        
        # List for tuples of start and end positions of events
        allevents=[]

        # Load data from data file
        data = np.load(directory+'/'+datadir+'/'+datafile)
        
        # Correct for arbitrary, fixed calibration performed in script
        # reading binary data from trap output files
        data[:,1]*=(cal_b1/(80.0))
        data[:,2]*=(cal_b2/(80.0))

        # Position data from the two beads
        data1=data[:,1]
        data2=data[:,2]
        
        # Correct the sign of the stage position
        data[:,5]*=-1

        # Loop over the two beads
        for beadnumber in [1,2]:
            events=[]
            # Loop over all data in 
            for datai in range(0,int((shape(data)[0]/10000.0))):
                #print 1.0*datai/int((shape(data)[0]/10000.0)), datafile, directory

                # Array of time in units of data points
                t=np.arange(datai*10000,(datai+1)*10000)

                # Fit sine-function to piezo position, to determine phase, assuming known period T
                parsinit=array([20.0,0.0,0.0])
                pars=leastsq(lambda pars: data[t,5]-(pars[0]*sin(t*2*pi/T+pars[1])+pars[2]),parsinit)[0]
                
                # Coerce amplitude into positive number
                if pars[0]<0.0:
                    piezoamp=abs(pars[0])
                    piezophase=pars[1]+pi
                else:
                    piezoamp=pars[0]
                    piezophase=pars[1]

                # Eliminate 2*pi degeneracy 
                piezophase=mod(piezophase,2*pi)

                # Fitted sine function and its corresponding cosine
                piezofit=piezoamp*sin(t*2*pi/T+piezophase)+pars[2]
                piezofit2=piezoamp*sin(t*2*pi/T+piezophase+pi/2.0)+pars[2]
                
                # Calculate wavelet-sine and wavelet-cosine transforms
                # and in turn the amplitude and phase of the signals
                pdata1=data[t,beadnumber]*(piezofit-mean(piezofit))
                pdata2=data[t,beadnumber]*(piezofit2-mean(piezofit2))
                pdata=sqrt(pdata1**2+pdata2**2)

                y=convolve(pdata,ones(T)/T,mode='valid')
                tsub=t[0]+arange(len(y))+(T-1)/2.0

                y1=2*convolve(pdata1,ones(T)/T,mode='valid')
                y2=2*convolve(pdata2,ones(T)/T,mode='valid')
                y=sqrt(y1**2+y2**2)        

                yamp=y/abs(piezoamp)                

                test=zeros(len(y))
                for n in range(len(y)):
                    if y1[n]>0.0 and y2[n]>0.0:
                        test[n]=arctan(y2[n]/y1[n])
                    elif y1[n]<0.0 and y2[n]>0.0:
                        test[n]=arctan(y2[n]/y1[n])+pi
                    elif y1[n]>0.0 and y2[n]<0.0:
                        test[n]=arctan(y2[n]/y1[n])
                    elif y1[n]<0.0 and y2[n]<0.0:
                        test[n]=arctan(y2[n]/y1[n])-pi

                yphase=test

                # Locate events using phase and amplitude
                if beadnumber==1:
                    binding1=yamp > amp_threshold1
                else:
                    binding1=yamp > amp_threshold2
                        
                binding2=yphase < phase_threshold

                # Require that both criteria are satisfied simultaneously
                binding=binding1*binding2
                
                # Get binding times
                tbinding=tsub[binding]
                tbinding=tbinding.astype(int)

                # Find groups of consecutive time points
                groups=[]
                for k, g in groupby(enumerate(tbinding), lambda (i,x):i-x):
                    groups.append(map(itemgetter(1), g))

                for element in groups:
                    tbinding=element
                    events.append((tbinding[0],tbinding[-1]))
            
            # Merge events if they happen to be located over 
            # a 10,000 data point break in the sequence
            n=0
            tempevents=[]
            while n < len(events)-1:
                if mod(events[n][1]-10000+int(T/2+1),10000)==0 and \
                    mod(events[n+1][0]-int(T/2-1),10000)==0:
                    tempevents.append((events[n][0],events[n+1][1]))
                    n+=2
                else:
                    tempevents.append(events[n])
                    n+=1
            events=tempevents

            if events!=[]:
                allevents+=events
            

        events=allevents
        events.sort(key=lambda tup: tup[0])
        
        # Merge events from the two beads if they overlap in time
        if events!=[]:
            merged=True
            while any(merged):
                mergedevents=[]
                merged=repeat(False,len(events)-1)
                n=0
                while n<len(events):
                    if n<len(events)-1 and (events[n+1][0]<events[n][1]):
                        mergedevents.append((events[n][0],max(events[n+1][1],events[n][1])))
                        merged[n]=True
                        n+=2
                    else:
                        mergedevents.append((events[n][0],events[n][1]))
                        n+=1
                events=mergedevents

            # Ignore a possible early event
            if events[0][0]< 5*T:
                events=events[1:]

        nevents_all+=len(events)

        # Loop over possible events
        for nevent in range(len(events)):
            event=events[nevent]

            # Test if threshold criteria is passed for more than one period
            if event[1]-event[0]>T: 

                try:

                    bindingtimes_Fourier=[]
                    F0_Fourier=[]
                    dF_Fourier=[]
                    dFcor_Fourier=[]
                    phi_Fourier=[]

                    avchi2=[]

                    # Prepare to plot results for duration of event +/- 25 periods
                    figure(1,figsize=(16,10))
                    starttime=event[0]
                    endtime=event[1]
                    tevent=np.arange(starttime,endtime)
                    tplot=np.arange(max(0,starttime-25*int(T)),min(endtime+25*int(T),shape(data)[0]))

                    # Plot position of bead 1
                    subplot(511)
                    plot(tplot,data[tplot,1],linestyle='-',color='k',markersize=1,linewidth=0.5)
                    plot(tevent,data[tevent,1],'y-',markersize=1,linewidth=1.0)
                    # Plot position of bead 2
                    subplot(512)
                    plot(tplot,data[tplot,2],linestyle='-',color='gray',markersize=1,linewidth=0.5)
                    plot(tevent,data[tevent,2],'c-',markersize=1,linewidth=1.0)
                    # Plot position of piezo-stage
                    parsinit=array([20.0,0.0,0.0])
                    pars=leastsq(lambda pars: data[tplot,5]-(pars[0]*sin(tplot*2*pi/T+pars[1])+pars[2]),parsinit)[0]

                    if pars[0]<0.0:
                        piezoamp=abs(pars[0])
                        piezophase=pars[1]+pi
                    else:
                        piezoamp=pars[0]
                        piezophase=pars[1]

                    piezophase=mod(piezophase,2*pi)

                    piezofit=piezoamp*sin(tplot*2*pi/T+piezophase)+pars[2]
                    piezofit2=piezoamp*sin(tplot*2*pi/T+piezophase+pi/2.0)+pars[2]

                    subplot(513)
                    plot(tplot,data[tplot,5],'k-',markersize=1,linewidth=0.5)
                    plot(tplot,piezofit,'g-',markersize=1,linewidth=1.0)
                                        
                    # Redo analysis for approved events (not optimal)
                    for beadnumber in [1,2]:

                        if beadnumber==1:
                            cal_k=cal_k1
                            yampthreshold=amp_threshold1
                        elif beadnumber==2:
                            cal_k=cal_k2
                            yampthreshold=amp_threshold2
                                                        
                        t=tplot
                        pdata1=data[t,beadnumber]*(piezofit-mean(piezofit))
                        pdata2=data[t,beadnumber]*(piezofit2-mean(piezofit2))
                        pdata=sqrt(pdata1**2+pdata2**2)

                        y=convolve(pdata,ones(T)/T,mode='valid')
                        tsub=t[0]+arange(len(y))+(T-1)/2.0

                        y1=2*convolve(pdata1,ones(T)/T,mode='valid')
                        y2=2*convolve(pdata2,ones(T)/T,mode='valid')
                        y=sqrt(y1**2+y2**2)        

                        yamp=y/abs(piezoamp)

                        test=zeros(len(y))
                        for n in range(len(y)):
                            if y1[n]>0.0 and y2[n]>0.0:
                                test[n]=arctan(y2[n]/y1[n])
                            elif y1[n]<0.0 and y2[n]>0.0:
                                test[n]=arctan(y2[n]/y1[n])+pi
                            elif y1[n]>0.0 and y2[n]<0.0:
                                test[n]=arctan(y2[n]/y1[n])
                            elif y1[n]<0.0 and y2[n]<0.0:
                                test[n]=arctan(y2[n]/y1[n])-pi

                        yphase=test
                        
                        # Calculate period to use for averaging
                        select=(tsub>starttime+T/2.0)*(tsub<endtime-T/2.0)

                        # Calculate average amplitude and phase of bound state
                        boundamplevel=mean(yamp[select])
                        boundphaselevel=mean(yphase[select])

                        # Plot thresholds
                        subplot(514)
                        plot(tsub,yamp,linestyle='-',color=plotcolors[beadnumber-1])
                        hlines(amp_threshold1,tplot[0],tplot[-1],linestyle='dotted')
                        hlines(amp_threshold2,tplot[0],tplot[-1],linestyle='dotted')
                        plot(tsub[select],yamp[select],color=highlightcolors[beadnumber-1],linestyle='-',markersize=1,linewidth=1.0)
                        hlines(boundamplevel,starttime,endtime,linestyle='dashed',color='black')

                        subplot(515)
                        plot(tsub,yphase,linestyle='-',color=plotcolors[beadnumber-1])
                        plot(tsub[select],yphase[select],color=highlightcolors[beadnumber-1],linestyle='-',markersize=1,linewidth=1.0)
                        hlines(phase_threshold,tplot[0],tplot[-1],linestyle='dotted')
                        hlines(boundphaselevel,starttime,endtime,linestyle='dashed',color='black')
                        ylim(-pi,pi)

                                                
                        # Find surrounding intervals that do not overlap with other events
                        tunbound1=max(starttime-25*int(T),0)
                        if nevent>=1 and tunbound1<events[nevent-1][1]:
                            tunbound1=events[nevent-1][1]

                        tunbound2=min(endtime+25*int(T),shape(data)[0])
                        if nevent<=len(events)-2 and tunbound2>events[nevent+1][0]:
                            tunbound2=events[nevent+1][0]

                        tunbound_before=arange(tunbound1,starttime-int(T))
                        tunbound_after=arange(endtime+int(T),tunbound2)
                        tunbound=arange(max(starttime-4*int(T),0),min(endtime+4*int(T),shape(data)[0]))
                        tunbound_fit=append(tunbound_before,tunbound_after)

                        # Determine parameters using wavelets
                        deltax=mean(yamp[select])
                        phi=mean(yphase[select])
                        
                        dF_Fourier.append(mean(yamp[select])*cal_k)
                        phi_Fourier.append(mean(yphase[select]))

                        # Find average amplitude of unbound states
                        select1=(tsub>tunbound1)*(tsub<starttime-T)
                        select2=(tsub>endtime+T)*(tsub<tunbound2)
                        select=np.logical_or(select1,select2)
                        unboundamplevel=mean(yamp[select])

                        # Locate interval of increasing amplitude
                        select=(tsub>starttime-T)*(tsub<min(starttime+T/2.0,(starttime+endtime)/2.0))
                        rise=yamp[select]
                        trise=tsub[select]

                        t1s=[]
                        zero_crossings = np.where(np.diff(np.sign(rise-(boundamplevel+unboundamplevel)/2.0)))[0]
                        for element in zero_crossings:
                            vlines(element,-10,10)
                            t1s.append(trise[element])
                        
                        # Binding times are calculated as full-width-at-half-max
                        # In case of multiple candidate times for binding, use the one closest to initial value
                        t1_Fourier=t1s[argmin(abs((t1s-starttime)))]

                        # Save binding time
                        bindingtimes_Fourier.append(t1_Fourier)

                        # Locate interval of decreasing amplitude
                        select=(tsub>max(endtime-T/2.0,(starttime+endtime)/2.0))*(tsub<endtime+T)
                        fall=yamp[select]
                        tfall=tsub[select]

                        t2s=[]
                        zero_crossings = np.where(np.diff(np.sign(fall-(boundamplevel+unboundamplevel)/2.0)))[0]
                        for element in zero_crossings:
                            vlines(element,-10,10)
                            t2s.append(tfall[element])

                        # Binding times are calculated as full-width-at-half-max
                        # In case of multiple candidate times for binding, use the one closest to initial value
                        t2_Fourier=t2s[argmin(abs((t2s-endtime)))]
                        
                        for panel in [511,512]:
                            subplot(panel)
                            if beadnumber==1:
                                vlines(t1_Fourier,ylim()[0],ylim()[1],linestyle='dashed',color='k')
                                vlines(t2_Fourier,ylim()[0],ylim()[1],linestyle='dashed',color='k')
                            elif beadnumber==2:
                                vlines(t1_Fourier,ylim()[0],ylim()[1],linestyle='dashed',color='gray')
                                vlines(t2_Fourier,ylim()[0],ylim()[1],linestyle='dashed',color='gray')
                                

                        # Save unbinding time
                        bindingtimes_Fourier.append(t2_Fourier)

                        # Fit harmonic function to binding region
                        parsinit=array([20.0,piezophase,0.0])
                        pars=leastsq(lambda pars: data[tevent,beadnumber]-(pars[0]*sin(tevent*2*pi/T+pars[1])+pars[2]),parsinit)[0]
                        subplot(510+beadnumber)
                        sinebinding=pars[0]*sin(tunbound*2*pi/T+pars[1])+pars[2]
                        plot(tunbound,sinebinding,'g-',lw=1.0)

                        # Coerce amplitude into a positive number and adjust phase accordingly
                        if pars[0]<0.0:
                            amp_bound=abs(pars[0])
                            phase_bound=pars[1]+pi
                        else:
                            amp_bound=pars[0]
                            phase_bound=pars[1]

                        piezophase=mod(piezophase,2*pi)
                        offset_bound=pars[2]


                        # Calculate average chi-squared for the bound state
                        t12=np.arange(int(t1_Fourier),int(t2_Fourier))
                        dev=data[t12,beadnumber]-(amp_bound*sin(t12*2*pi/T+phase_bound)+offset_bound)
                        ssdev=sum(dev**2)
                        avssdev=ssdev/(t2_Fourier-t1_Fourier)
                        avchi2.append(avssdev)
                        
                        # Fit sine-function to early detached state
                        figure(1)
                        parsinit=array([5.0,piezophase,0.0])
                        pars=leastsq(lambda pars: data[tunbound_fit,beadnumber]-(pars[0]*sin(tunbound_fit*2*pi/T+pars[1])+pars[2]),parsinit)[0]
                        subplot(510+beadnumber)
                        sineunbound_fit=pars[0]*sin(tunbound_fit*2*pi/T+pars[1])+pars[2]
                        sineunbound=pars[0]*sin(tunbound*2*pi/T+pars[1])+pars[2]
                        sineunbound_before=pars[0]*sin(tunbound_before*2*pi/T+pars[1])+pars[2]
                        sineunbound_after=pars[0]*sin(tunbound_after*2*pi/T+pars[1])+pars[2]
                        plot(tunbound_before,sineunbound_before,'b-',lw=2.0)
                        plot(tunbound_after,sineunbound_after,'b-',lw=2.0)

                        tmed=arange(starttime-T,endtime+T)
                        sineunbound_med=pars[0]*sin(tmed*2*pi/T+pars[1])+pars[2]
                        plot(tmed,sineunbound_med,'b--',lw=2.0)

                        if pars[0]<0.0:
                            amp_unbound=abs(pars[0])
                            phase_unbound=pars[1]+pi
                        else:
                            amp_unbound=pars[0]
                            phase_unbound=pars[1]

                        offset_unbound=pars[2]


                        # Determine F0 from raw trajectory
                        t=tplot
                        pdata=data[t,beadnumber]
                        y=convolve(pdata,ones(T)/T,mode='valid')
                        tsub=t[0]+arange(len(y))+(T-1)/2.0
                        subplot(510+beadnumber)

                        if endtime-starttime>T:
                            select=(tsub>starttime+T/2.0)*(tsub<endtime-T/2.0)
                        else:
                            select=(tsub>starttime)*(tsub<endtime)
                            
                        boundlevel=mean(y[select])
                        hlines(boundlevel,starttime,endtime,color='k',linestyle='dashed')
                        select1=(tsub>tunbound1)*(tsub<starttime-T)
                        select2=(tsub>endtime+T)*(tsub<tunbound2)
                        select=np.logical_or(select1,select2)
                        unboundlevel=mean(y[select])

                        F0_Fourier.append((boundlevel-unboundlevel)*cal_k)

                    # Polish the plots
                    for panel in [511,512,513,514,515]:
                        subplot(panel)
                        tlim=np.arange(max(0,starttime-5*int(T)),min(endtime+5*int(T),shape(data)[0]))
                        xlim(tlim[0],tlim[-1])
                        if panel==511:
                            text(xlim()[0]+0.95*(xlim()[1]-xlim()[0]),ylim()[1]-20,str(round(avchi2[0],2)))
                        elif panel==512:
                            text(xlim()[0]+0.95*(xlim()[1]-xlim()[0]),ylim()[1]-20,str(round(avchi2[1],2)))
                        
                    subplot(511)
                    ylabel('$x_1$ (nm)')
                    subplot(512)
                    ylabel('$x_2$ (nm)')
                    subplot(513)
                    ylabel('$x_\mathrm{stage}$ (nm)')
                    subplot(514)
                    ylabel('Amplitude (nm)')
                    subplot(515)
                    ylabel('Phase (rad)')
                    xlabel('Time (frames)')

                    # Save the diagnostics figure
                    savefig(resultsdir+'/'+'event'+str(nevent_save)+'.png')
                    close('all')

                    F0_all_Fourier.append(F0_Fourier)
                    dF_all_Fourier.append(dF_Fourier)
                    bindingtimes_all_Fourier.append(bindingtimes_Fourier)

                    avchi2_all.append(avchi2)

                    nevent_save+=1

                except (IndexError,ValueError,TypeError,RuntimeError):
                    pass


    # Save the results to files
    np.savetxt(directory+'/'+'bindingtimes_Wavelet.txt',bindingtimes_all_Fourier)
    np.savetxt(directory+'/'+'F0_Wavelet.txt',F0_all_Fourier)
    np.savetxt(directory+'/'+'dF_Wavelet.txt',dF_all_Fourier)
    np.savetxt(directory+'/'+'avchi2.txt',avchi2_all)

