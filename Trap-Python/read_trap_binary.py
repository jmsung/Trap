#!/usr/bin/env python

#Tural Aksel
#Spudich lab

#05/17/13


import numpy as np
import struct
import os
import sys
import argparse

#Global Variables
PRECISION = 1E-8

def read_binary_data(fname,time_window=30.0): #60.0 sec was original
        '''
                Function for parsing optical trap data written by labview
            
                Inputs
                ------
                        fname           :binary file name
                        time_window     :data will be written out in parts of time_window(in seconds)
                Outputs
                -------
                        Data is written in parts of time_window(s). The output matrix for every time window is stored *.npy format.
                        The matrix has three columns (col1: time(s), col2:channel1 displacement(nm), col3:channel2 displacement(nm))
                        If fname is xxx.bin, the resulting output files will be in the form xxx/%d.sty (%d starting from 0 in folder xxx). 
                
        '''

        print "start parsing"
        
        #Parse the filename
        base,ext = os.path.splitext(fname)
        print base + ext
        
        #Create directory for the extracted data. bin data should have .bin extension else it thinks dir already exists!
        if not os.path.exists(base):
            os.mkdir(base)
            print "made directory" + base
        
        #If the raw file is already processed, do not process
        if os.path.exists(base+'/processed.txt'):
            return
        
        #Read raw data directly coming from trap
        f = open(fname, 'rb') ####### must use 'rb' when opening binary file in windows!
        
        #Extract the useful information from the header
        junk_header=struct.unpack('>i4',f.read(4))[0]
        
        #length of the remainder of the header
        header     = struct.unpack('>i4',f.read(4))[0]
        
        #length of the channel list to follow (M)
        M          = struct.unpack('>i4',f.read(4))[0]
        
        status     = f.seek(12+M, 0)
        #length of the group channel setting to follow (N)
        N = struct.unpack('>i4',f.read(4))[0]
        #Get number of channels
        channels_number = struct.unpack('>i4',f.read(4))[0]-12;
        #channels_number = struct.unpack('>i4',f.read(4))[0];
        print "Number of channels:" + repr(channels_number)
        
        #Read channel detatils
        upper           = np.zeros(channels_number)
        lower           = np.zeros(channels_number)
        lrange          = np.zeros(channels_number)
        coupling        = np.zeros(channels_number)
        term_cfg        = np.zeros(channels_number)
        mult_array_size = np.zeros(channels_number)
        mult0           = np.zeros(channels_number)
        mult1           = np.zeros(channels_number)
        mult2           = np.zeros(channels_number)
        mult3           = np.zeros(channels_number)
        
        for i in range(channels_number):
            channum_string     = struct.unpack('>i4',f.read(4))[0]
            channum            = struct.unpack('>'+channum_string*'c',f.read(channum_string))
            upper[i]           = struct.unpack('>d8',f.read(8))[0]
            lower[i]           = struct.unpack('>d8',f.read(8))[0]
            lrange[i]          = struct.unpack('>d8',f.read(8))[0]
            coupling[i]        = struct.unpack('>i4',f.read(4))[0]
            term_cfg[i]        = struct.unpack('>i4',f.read(4))[0]
            mult_array_size[i] = struct.unpack('>i4',f.read(4))[0]
            mult0[i]           = struct.unpack('>d8',f.read(8))[0]
            mult1[i]           = struct.unpack('>d8',f.read(8))[0]
            mult2[i]           = struct.unpack('>d8',f.read(8))[0]
            mult3[i]           = struct.unpack('>d8',f.read(8))[0]
        
        status = f.seek(12+4+M+N, 0)
        
        #Scan rate
        scan_rate = struct.unpack('>d8',f.read(8))[0]
        print "Scan rate: " + repr(scan_rate) + " Hz"
        
        #Interchannel delay
        delay = struct.unpack('>d8',f.read(8))[0]
        
        #Now read the real data
        data_start = header+12+4                             
        buffer_num = 2000 #2000 for 40kHz sampling, 1000 for 50kHz sampling
        buffer_size = buffer_num/scan_rate  #in sec
        print "Buffer size: " + repr(buffer_size) + " sec"

        #bytes per one block of buffer, data = int16 = 2 bytes, spacer = 8 bytes
        block_size      = channels_number*buffer_num*2 + 8
        status          = f.seek(0,2)
        total_data_size = f.tell() - data_start
        num_block       = total_data_size / block_size
        total_time      = np.floor(buffer_size * num_block)
        print "Total time: " + repr(total_time) + " sec"
    
        #Move back to starting point
        status = f.seek(data_start, 0)                        
        
        #Starting time(offset)
        begin_time = 0
        
        #Shift by the time offset 
        status = f.seek(2*begin_time*scan_rate*channels_number + 8*begin_time/buffer_size, 1)
       
        # Moving 2 bytes for each time point leading up to the begin_time - and multiply by the number of channels  
        # Labview queue injects 8 byte spacer data between each buffer interval. Have to skip those sections
        counter        = buffer_size
        traces = np.resize(np.fromfile(f,'>h2',channels_number*buffer_num),(buffer_num,channels_number))
        
        #Assign also a global counter
        global_counter = buffer_size
        time_left      = total_time - global_counter
        file_counter   = 0
        while global_counter < total_time:
            while counter < time_window - PRECISION and time_left - PRECISION > buffer_size:
                status          = f.seek(8, 1)
                new_traces      = np.resize(np.fromfile(f,'>h2',channels_number*buffer_num),(buffer_num,channels_number))
                
                #If it is the first time we read a buffer, initialize traces to new_traces
                if counter == 0:
                    traces = new_traces
                else:
                    traces = np.vstack((traces,new_traces))
                
                #Update the counter variables
                counter        += buffer_size
                global_counter += buffer_size
                time_left      -= buffer_size
            
            print "Time left: " + repr(time_left) + " sec"
            
            #If nothing is read(the final round), exit the loop
            if counter == 0.0:
                np.savetxt(base+'/processed.txt',np.zeros(1))
                break
            time = begin_time+1.0/scan_rate*np.arange(np.round(counter*scan_rate))
            
            #Now, divide the data among the channels
            true_traces = np.zeros(traces.shape)
            for i in range(channels_number):
                true_traces[:,i] = mult0[i]+traces[:,i]*mult1[i] + traces[:,i]**2 * mult2[i]+traces[:,i]**3 * mult3[i]
            
            x1_raw  = true_traces[:,0]
            x2_raw = true_traces[:,2]
            
            x1_raw   = x1_raw  - np.mean(x1_raw)
            x2_raw  = x2_raw - np.mean(x2_raw)
        
            #Piezo stage position in nm -10000 (for volt to nm conversion - 1 volt = 10 um)
            pzt_x    = true_traces[:,6]*10000
            
            #Multiplication factor b1,b2 :Volt to nm conversion
            b1=80
            b2=80
            
            x1_nm = x1_raw*b1
            x2_nm = x2_raw*b2
            
            #Multiplation for k1,k2 :nm to force(pN) conversion
            k1 = 0.1
            k2 = 0.1
            
            x1_pN = x1_nm*k1
            x2_pN = x2_nm*k2 
            
            #Combine the data file in matrix form
            #time, bead1-position_x, bead2-position_x, bead1-force, bead2-force, piezo-stage1-position_x, piezo-stage2-position_x
            print "Saving " + repr(file_counter) + ".npy..."
            final_data = np.hstack((np.vstack(time),np.vstack(x1_nm),np.vstack(x2_nm),np.vstack(x1_pN),np.vstack(x2_pN),np.vstack(pzt_x)))
            np.save(base+'/'+str(file_counter)+'.npy',final_data)
            
            #Update the counter variables
            begin_time   += counter
            counter       = 0.0
            file_counter +=1
        
        f.close()
        
if __name__=="__main__":

        # Example values
        directory='C:/Users/chaoliu/Documents/Spudich Lab/Data/20151218_WT808_oscillation'
        fname='cal1_WT808_m2.bin'
        read_binary_data(directory+'/'+fname)


        
