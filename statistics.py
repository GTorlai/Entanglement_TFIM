import gzip
import cPickle
import numpy as np
import argparse
import glob
import math as m

#-------------------------------------------------

def get_average(data):

    """ Compute averages of data """
    
    n_meas = len(data)
    
    if len(data.shape) != 1:
        n_obs  = len(data[0])
        avg = np.zeros((n_obs))
        for i in range(n_meas):
            for j in range(n_obs):
                avg[j] += data[i,j]
 
    else:
        avg = 0.0
        for i in range(n_meas):
            avg += data[i]
        
    avg /= 1.0*n_meas
    return avg

#------------------------------------------------

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-L', help='Linear size of the model',type=int)
    args = parser.parse_args()

    outputName = 'data/randomIsing/collapse_L'
    outputName += str(args.L)
    outputName += '.dat'
    
    outputFile = open(outputName,'w')
    
    h = ["0.9500","0.9550","0.9600","0.9650","0.9700",
         "0.9750","0.9800","0.9825","0.9850","0.9875",
         "0.9900","0.9910","0.9920","0.9930","0.9940",
         "0.9950","0.9960","0.9970","0.9980","0.9990",
         "1.0000",
         "1.0010","1.0020","1.0030","1.0040","1.0050",
         "1.0060","1.0070","1.0080","1.0090","1.0100",
         "1.0125","1.0150","1.0175","1.0200","1.0250",
         "1.0300","1.0350","1.0400","1.0450","1.0500"]


    for i in range(41):

        #h = 0.0 + 0.05*i

        dataFileName = 'data/randomIsing/L'
        dataFileName += str(args.L)
        dataFileName += '/'
        dataFileName += 'randomIsing_L'
        dataFileName += str(args.L)
        dataFileName += '_h'
        dataFileName += h[i]
        dataFileName += "collapse_data.txt"
        #dataFileName += str(h)
        #if ((i%2) == 0):
        #    dataFileName += '0'
        #dataFileName += '00collapsedata.txt'
        
        print dataFileName
     
        dataFile = open(dataFileName,'r')
        data = np.loadtxt(dataFile)
        avg = get_average(data)
        #err = np.std(data,0)/m.sqrt(data[:].shape[0]-1)
        err = np.std(data,axis=0)/m.sqrt(len(data))
 
        outputFile.write('%f   ' % float(h[i]))

        for k in range(2):
            outputFile.write('%E   ' % avg[k])
            outputFile.write('%E   ' % err[k])
            
        outputFile.write('\n') 
 



