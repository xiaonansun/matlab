import numpy as np
import os
import scipy
import scipy.io
import sys

def combined(savepath):
    stat = np.load(os.path.join(savepath,'stat.npy'), allow_pickle=True)
    ops = np.load(os.path.join(savepath,'ops.npy'), allow_pickle=True)

    statpath = os.path.join(savepath,'stat.mat')
    opspath = os.path.join(savepath,'ops.mat')

    ops=ops.tolist()
    ops['date_proc']=[]
    ops=np.array(ops)

    scipy.io.savemat(statpath, {'stat': stat})
    scipy.io.savemat(opspath, {'ops': ops})

savepath = sys.argv[1]
print('Parsing {0}'.format(savepath))
combined(savepath)