"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Trap calibration (Jongmin Sung)

class Data()
- path, name


"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

from __future__ import division, print_function, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import os
from PIL import Image

pixel_nm = 135 
Bead_nm = 1000
Bead_px = Bead_nm / pixel_nm
PM_mV = 20 


class Data(object):
    def __init__(self):
        data_path = os.getcwd()
        path_split = data_path.split('\\')      
        self.data_name = path_split[len(path_split)-1]      
        print('Data name = %r' %(self.data_name))
        file_list = os.listdir(data_path) 

        for i in range(len(file_list)):
            if file_list[i][-3:] == 'tif':
                movie_name = file_list[i]

        movie_path = data_path + '\\' + movie_name   
        movie = Image.open(movie_path)
        self.n_frame = movie.n_frames
        print('Number of frames = %d' %(self.n_frame))
        self.n_row = movie.size[1]
        self.n_col = movie.size[0]
    
    
        self.I = np.zeros((self.n_frame, self.n_row, self.n_col), dtype=int)
        for i in range(self.n_frame): 
            movie.seek(i) # Move to i-th frame
            self.I[i,] = np.array(movie, dtype=int)
        


    def analysis(self):
        pass
        
    def plot(self):
        pass



# Start  
plt.close('all')
data = Data()
data.analysis()
data.plot()

