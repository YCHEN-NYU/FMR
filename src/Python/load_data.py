
# load Packages & Settings
import os
import numpy as np


def get_files(rootdir, suffix):
    file_list = []
    for subdir, dirs, files in os.walk(rootdir):
        for file in files:
            if file.endswith(suffix):
                file_list.append(os.path.join(subdir, file))
    file_list.sort(key=lambda x: x)
    return file_list

# load data from file and return it as numpy-array
def load_data(filename):
    arr = []
    temp = []
    for line in [line for line in open(filename, 'r')]:
        # parse each line and convert to floating
        temp = [float(element) for element in line.split(',') if (element != '' and element != ' ')] 
        temp =  np.array(temp[:]) 
        # shape into columns
        temp.reshape(-1,1)
        arr.append(temp)
    # release the temp variable
    del temp
    return np.array(arr)