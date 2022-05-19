import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pysam
from scipy.sparse import coo_matrix
import time
import numpy as np
import tensorflow as tf
import pysam
import matplotlib.pyplot as plt
import multiprocessing
import os
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Value, Array
import time

import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pysam
from scipy.sparse import coo_matrix
import time
import numpy as np
import tensorflow as tf
import pysam
import matplotlib.pyplot as plt

def fx(alist, blist, clist, rowcount):
    for b in blist:
        alist.append(b)
        clist.append(rowcount)
def chioce_top18(tensor):
    batch_size, window_size, rowcount = tensor.shape[0], tensor.shape[1], tensor.shape[2]
    tensor = tf.concat([tensor, tf.zeros([batch_size, window_size, 18])], axis = 2)
    return tf.reshape(tf.gather(tensor, tf.argsort(tf.reduce_sum(tensor, 1, keepdims = True), axis = 2), axis=2, batch_dims=1)[:,:,:,-18:], [tensor.shape[0], tensor.shape[1], 18, 1])
    

def labeldata(vcfpath, contig, start, end):
  goldl = []
  window_size = 200
  index = start + np.column_stack((np.arange(0, window_size * (int((end - start - 1) / window_size) + 1), window_size).reshape((int(( end - start - 1) / window_size) + 1), 1), np.arange(0, window_size * (int((end - start - 1) / window_size) + 1), window_size).reshape((int((end - start - 1) / window_size) + 1), 1) + window_size))
  if('chr' in contig):
    contig = contig[3:]
  for rec in pysam.VariantFile(vcfpath).fetch():

    if(rec.contig != contig):
      continue            
    if((rec.info['SVTYPE'] == 'DEL')):
      goldl.append([rec.start, rec.stop, rec.stop - rec.start, 1])
        
  
    
  goldl = (pd.DataFrame(goldl).sort_values([0, 1]).values).astype('float64')


  y = []
  for rec in index:
        
    if(((goldl[:,1:2] > rec[0]) & (goldl[:,:1] < rec[1])).sum() != 0):
      y.append((((goldl[:,1:2] > rec[0]) & (goldl[:,:1] < rec[1])) * goldl[:,3:]).sum())


    else:
      y.append(0)
  return (np.array(y)>0).astype('float32')




def cigarinfo(cigararray, refstartarray, start, end, reccount, cigarweight):

    a = (~((cigararray[:,0] == 1) | (cigararray[:,0] == 4))).reshape(cigararray.shape[0], 1) * cigararray
    a1 = a[:,1].reshape((reccount, cigarweight))
    a1 = (a1 @ np.triu(np.ones((a1.shape[1] , a1.shape[1])), k=0)) + refstartarray
    a2 = (a1 != 0) * (np.arange(a1.shape[0]).reshape(a1.shape[0], 1))
    a = np.concatenate([cigararray, a2.reshape((cigararray.shape[0], 1)), a1.reshape([cigararray.shape[0], 1])], axis = 1)
    

    return a[(start <= a[:,-1]) & (a[:,-1] < end)]



def fx(alist, blist, clist, rowcount):
    for b in blist:
        alist.append(b)
        clist.append(rowcount)



def chioce_top18(tensor, geno = False):
    batch_size, window_size, rowcount = tensor.shape[0], tensor.shape[1], tensor.shape[2]
    if(rowcount < 32):
        tensor = np.concatenate([tensor, np.zeros([batch_size, window_size, 32-rowcount], dtype = np.float32)], axis = 2)
    if(geno == True):
        return tensor[:,:,:32].reshape(tensor.shape[0], tensor.shape[1], 32, 1)
    index = np.argsort(tensor.sum(axis = 1, keepdims = True), axis = 2)[:,:,-32:].reshape(-1, 32)
    l = []
    for loc in range(batch_size):
        l.append(tensor[loc][:,index[loc]])
    return np.concatenate(l, axis = 0).reshape(tensor.shape[0], tensor.shape[1], 32, 1)

def delinfo(bamfile, contig, start, end, window_size = 200):
    collist = []
    rowlist = []
    rowend = []
    count = 0
    maxend = 0
    minstart = 999999999
    rowcount = 0
    cigararray = []
    readstartandend = []
    maxlen = 0
    nooverlap = False
    for AlignedSegment in bamfile.fetch(contig, start, end):


        seqqosition = np.array(AlignedSegment.get_reference_positions())
        if((seqqosition[-1] - seqqosition[0])<1000):
            continue
        nooverlap = True
        count += 1
        
        readstartandend.append([seqqosition[0] - start])
        seqqosition = seqqosition[(seqqosition>=start) & (seqqosition<end)]
        cstart = seqqosition[0]
        cend = seqqosition[-1]
        seqqosition = set(range(cstart, cend+1)) - set(seqqosition)
        newrow = True
        loc = -1
        
        cigararray.append((np.array(AlignedSegment.cigartuples).flatten()))
        if(maxlen<cigararray[-1].size):
            maxlen = cigararray[-1].size
            
        for oneend in rowend:
            loc += 1
            if(oneend<cstart):
                fx(collist, seqqosition, rowlist, loc)
                rowend[loc] = cend
                newrow = False
                break
        if(newrow == True):
            rowcount += 1
            rowend.append(cend)
            fx(collist, seqqosition, rowlist, len(rowend)-1)
        if(maxend<cend):
            maxend = cend
        if(minstart>cstart):
            minstart=cstart
    if(nooverlap == False):
        return False, 0, 0
    for loc in range(len(cigararray)):
        cigararray[loc].resize(maxlen, refcheck=False)
    cigararray = np.array(cigararray).astype('float32')
    reccount, cigarweight = cigararray.shape[0], int(cigararray.shape[1] / 2)
    cigararray = cigararray.reshape(int(cigararray.size / 2), 2)
    #print('cigarinfo')
    cigararray = cigarinfo(cigararray, np.array(readstartandend), 0, end - start, reccount, cigarweight)

    a = cigararray[(cigararray[:,0] == 1)]
    
    row  = np.array(rowlist)
    col  = np.array(collist)-minstart
    data = np.ones(col.size, dtype = np.float32)
    return True, (coo_matrix((data, (row, col)), shape=(len(rowend), end-start)).toarray()).T, coo_matrix((np.log(a[:, 1]), (a[:,2], a[:,3])), shape=(int((a[:,2].max()))+1, (end-start))).toarray().T
def insloc(rs, start, end, cigartuples, ins_collist, ins_rowlist, ins_data, row):
    loc = rs
    rec = False
    for item in cigartuples:
        if(item[0] not in [1, 4, 5, 6]):
            loc += item[1]
            if(rec == False and loc>=start):
                rec = True
            if(loc>=end):
                break
        elif(rec == True and item[0] == 1):
            ins_collist.append(loc)
            ins_rowlist.append(row)
            ins_data.append(item[1])
def fx(alist, blist, clist, rowcount):
    for b in blist:
        alist.append(b)
        clist.append(rowcount)
def delinfo(bamfile, contig, start, end, window_size = 200):
    collist = []
    rowlist = []
    rowend = []
    ins_collist = []
    ins_rowlist = []
    ins_data = []
    count = 0
    maxend = 0
    minstart = 999999999
    rowcount = 0

    maxlen = 0
    nooverlap = False
    for AlignedSegment in bamfile.fetch(contig, start, end):


        seqqosition = np.array(AlignedSegment.get_reference_positions())
        if((seqqosition[-1] - seqqosition[0])<1000):
            continue
        nooverlap = True
        count += 1
        
        rs = seqqosition[0]
        seqqosition = seqqosition[(seqqosition>=start) & (seqqosition<end)]
        cstart = seqqosition[0]
        cend = seqqosition[-1]
        if(maxend<cend):
            maxend = cend
        if(minstart>cstart):
            minstart=cstart
        seqqosition = set(range(cstart, cend+1)) - set(seqqosition)
        newrow = True
        loc = -1

        for oneend in rowend:
            loc += 1
            if(oneend<cstart):
                fx(collist, seqqosition, rowlist, loc)
                insloc(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, loc)
                rowend[loc] = cend
                newrow = False
                break
        if(newrow == True):
            rowcount += 1
            rowend.append(cend)
            fx(collist, seqqosition, rowlist, len(rowend)-1)
            insloc(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, len(rowend)-1)
        
    if(nooverlap == False):
        return False, 0, 0

    
    row  = np.array(rowlist)
    col  = np.array(collist)-minstart
    data = np.ones(col.size, dtype = np.float32)
    return True, (coo_matrix((data, (row, col)), shape=(len(rowend), end-start)).toarray()).T, coo_matrix((np.log(ins_data), (ins_rowlist, ins_collist)), shape=((len(rowend), end-start))).toarray().T
def insloc(rs, start, end, cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist,  row):
    loc = rs
    rec = False
    for item in cigartuples:
        if(item[0] not in [1, 4, 5, 6]):
            if(rec == False and loc>=start):
                rec = True
            if(rec == True and item[0] == 2):
                for bias in range(min(item[1], end - loc)):
                    del_collist.append(loc+bias)
                    del_rowlist.append(row)

                    
            loc += item[1]
            
            
            if(loc>=end):
                break
        elif(rec == True and item[0] == 1):
            ins_collist.append(loc)
            ins_rowlist.append(row)
            ins_data.append(item[1])
def fx(alist, blist, clist, rowcount):
    for b in blist:
        alist.append(b)
        clist.append(rowcount)
def delinfo(bamfile, contig, start, end, window_size = 200):
    collist = []
    rowlist = []
    del_collist = []
    del_rowlist = []
    del_data = []
    rowend = []
    ins_collist = []
    ins_rowlist = []
    ins_data = []
    count = 0
    maxend = 0
    minstart = 999999999
    rowcount = 0

    maxlen = 0
    nooverlap = False
    for AlignedSegment in bamfile.fetch(contig, start, end):


        nooverlap = True
        count += 1
        
        rs = AlignedSegment.reference_start

        cstart = max(rs, start)
        cend = min(AlignedSegment.reference_end, end)
        if(maxend<cend):
            maxend = cend
        if(minstart>cstart):
            minstart=cstart

        newrow = True
        loc = -1

        for oneend in rowend:
            loc += 1
            if(oneend<cstart):
                insloc(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist, loc)
                rowend[loc] = cend
                newrow = False
                break
        if(newrow == True):
            rowcount += 1
            rowend.append(cend)
            insloc(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist, len(rowend)-1)
        
    if(nooverlap == False):
        return False, 0, 0

    return True, (coo_matrix((np.ones(len(del_rowlist)), (del_rowlist, del_collist)), shape=(len(rowend), end-start)).toarray()).T.astype('float32'), coo_matrix((np.log(ins_data), (ins_rowlist, ins_collist)), shape=((len(rowend), end-start))).toarray().T.astype('float32')


def indelinfo(bamfile, contig, start, end, window_size=200, geno = False):

    st = time.time()
    SIGNAL, t1, t2 = delinfo(bamfile, contig, start, end)
    if(SIGNAL == False):
        return False, 0, 0

    return True, np.concatenate((chioce_top18(t1.reshape(-1, window_size, t1.shape[1]), geno), chioce_top18(t2.reshape(-1, window_size, t2.shape[1]), geno)), axis = -1), time.time() - st


def chioce_top18(tensor, geno = False):
    k = 32
    batch_size, window_size, rowcount = tensor.shape[0], tensor.shape[1], tensor.shape[2]
    if(rowcount < k):
        tensor = np.concatenate([tensor, np.zeros([batch_size, window_size, k-rowcount], dtype = np.float32)], axis = 2)
    if(geno == True):
        return tensor[:,:,:k].reshape(tensor.shape[0], tensor.shape[1], k, 1)
    index = np.argsort(tensor.sum(axis = 1, keepdims = True), axis = 2)[:,:,-k:].reshape(-1, k)
    l = []
    for loc in range(batch_size):
        l.append(tensor[loc][:,index[loc]])
    return np.concatenate(l, axis = 0).reshape(tensor.shape[0], tensor.shape[1], k, 1)
###########################
def insloc(rs, start, end, cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist,  row):
    loc = rs
    rec = False
    for item in cigartuples:
        if(item[0] not in [1, 4, 5, 6]):
            if(rec == False and loc>=start):
                rec = True
            if(rec == True and item[0] == 2):
                for bias in range(min(item[1], end - loc)):
                    del_collist.append(loc+bias)
                    del_rowlist.append(row)

                    
            loc += item[1]
            
            
            if(loc>=end):
                break
        elif(rec == True and item[0] == 1):
            ins_collist.append(loc)
            ins_rowlist.append(row)
            ins_data.append(item[1])
############################

def myshow(contig, bamfilepath, start, end, outputpath = '', vcfpath = '', window_size = 200):
       
    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 64)

    
    samplelocation = start + np.column_stack((np.arange(0, window_size * (int((end - start - 1) / window_size) + 1), window_size).reshape((int(( end - start - 1) / window_size) + 1), 1), np.arange(0, window_size * (int((end - start - 1) / window_size) + 1), window_size).reshape((int((end - start - 1) / window_size) + 1), 1) + window_size))
    end = samplelocation[-1, 1]  
    ######################################
    maxread = 32
    collist = []
    rowlist = []
    del_collist = []
    del_rowlist = []
    del_data = []
    rowend = []
    ins_collist = []
    ins_rowlist = []
    ins_data = []
    count = 0
    maxend = 0
    minstart = 999999999
    rowcount = 0

    maxlen = 0
    nooverlap = False
    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)
    for AlignedSegment in bamfile.fetch(contig, start, end):



        nooverlap = True
        count += 1
        
        rs = AlignedSegment.reference_start

        cstart = max(rs, start)
        cend = min(AlignedSegment.reference_end, end)
        
        if(maxend<cend):
            maxend = cend
        if(minstart>cstart):
            minstart=cstart

        newrow = True
        loc = -1

        for oneend in rowend:
            loc += 1
            if(oneend<cstart):
                insloc(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist, loc)
                rowend[loc] = cend
                newrow = False
                break
        if(newrow == True and (rowcount <= 500)):
            rowcount += 1
            rowend.append(cend)
            insloc(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist, len(rowend)-1)
        
    geno = False
    t1, t2 = (coo_matrix((np.ones(len(del_rowlist)), (del_rowlist, del_collist)), shape=(len(rowend), end-start)).toarray()).T, coo_matrix((np.log(ins_data), (ins_rowlist, ins_collist)), shape=((len(rowend), end-start))).toarray().T
    data = np.concatenate((chioce_top18(t1.reshape(-1, window_size, t1.shape[1]), geno), chioce_top18(t2.reshape(-1, window_size, t2.shape[1]), geno)), axis = -1)
    
    
    
    ##################################################

    mask = data.reshape(((end-start)//window_size, window_size * maxread*2)).sum(axis =1) != 0
    data = data[mask]
    index = np.arange(start, end, window_size)[mask]
    np.savez_compressed(outputpath+contig+','+str(start)+','+str(end)+',data', data = data, label = np.array([0]), index = index)
    
    



def labeldata(vcfpath, contig, start, end, window_size, index):
  goldl = []
  if('chr' in contig):
    contig = contig[3:]
  for rec in pysam.VariantFile(vcfpath).fetch():

    if(rec.contig != contig):
      continue            
    if((rec.info['SVTYPE'] == 'DEL')):
      goldl.append([rec.start, rec.stop, rec.stop - rec.start, 1])
        
  
    
  goldl = (pd.DataFrame(goldl).sort_values([0, 1]).values).astype('float64')


  y = []
  for rec in index:
        
    if(((goldl[:,1:2] > rec) & (goldl[:,:1] < (rec+window_size))).sum() != 0):
      y.append((((goldl[:,1:2] > rec) & (goldl[:,:1] < (rec+window_size))) * goldl[:,3:]).sum())


    else:
      y.append(0)
  return (np.array(y)>0).astype('float32')


def feature_matrix(bamfilepath, outputpath = '', vcfpath = '', window_size = 200, max_worker = 16):
    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)
    contig2length = {}
    block = 5000000
    includecontig = [str(i) for i in range(12, 23)]
    for count in range(len(bamfile.get_index_statistics())):
        contig2length[bamfile.get_index_statistics()[count].contig] = bamfile.lengths[count]
    for contig in includecontig:
        for AlignedSegment in bamfile.fetch(contig):
            start = AlignedSegment.reference_start
            break
        while((start + 2*block) <  contig2length[contig]):
            while(True):
                if(len(multiprocessing.active_children()) < max_worker):
                    print('working on contig = ', contig, start, start + block)
                    contig = str(contig)
                    multiprocessing.Process(target=myshow, args=(contig, bamfilepath, start, start + block, outputpath, vcfpath, window_size)).start()
                    break
                else:
                    time.sleep(2)
            start +=  block
        while(True):
            if(len(multiprocessing.active_children()) < max_worker):
                print('working on contig = ', contig, start, start + block)
                contig = str(contig)
                multiprocessing.Process(target=myshow, args=(contig, bamfilepath, start, start + block, outputpath, vcfpath, window_size)).start()
                break
            else:
                time.sleep(2)


import time
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Input, Dense, Conv2D, DepthwiseConv2D, BatchNormalization, Dropout, GlobalAveragePooling2D, Reshape, multiply, add, Activation
import multiprocessing
import os
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Value, Array
import time
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
import tensorflow as tf
import threading
#tf.config.set_visible_devices([], 'GPU')

tf.random.set_seed(123)


class convse(keras.layers.Layer):
    
    def __init__(self,filters, kernel_size, strides=(1, 1), padding='valid', activation='relu', se_ratio = 0.25):

        super(convse, self).__init__()
        self.conv = layers.Conv2D(filters, kernel_size, strides=strides, padding=padding, activation=activation)
        
        self.filters = filters
        self.kernel_size = kernel_size
        self.strides = strides
        
        self.filters_se = max(1, int(filters*se_ratio))
        self.conv_1 = Conv2D(filters=self.filters_se,
                    kernel_size=1,
                    padding='same',
                    use_bias=True)
        self.conv_2 = Conv2D(filters=self.filters,
                    kernel_size=1,
                    padding='same',
                    activation='sigmoid',
                    use_bias=True)

    def call(self, inputs):
        
        x = self.conv(inputs)



        se = GlobalAveragePooling2D()(x)
        se = Reshape((1, 1, self.filters))(se)


        se = self.conv_1(se)

        se = Activation(tf.nn.swish)(se)

        se = self.conv_2(se)


        return multiply([x, se])
    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[1], 1, self.filters)



def cnn_layer():
    model = keras.models.Sequential()

    model.add(tf.keras.layers.AveragePooling2D([2, 4]))
    model.add(convse(108, (1, 8), strides=(1, 8), padding='valid', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(64, (2, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(64, (3, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(80, (4, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([3,1]))

    model.add(convse(80, (3, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(64, (2, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))
    model.add(layers.Flatten())


    return model

def init_model():
    inputs = tf.keras.Input((None, 200, 32, 2))

    cnn_layer_object = cnn_layer()
    encoded_frames = tf.keras.layers.TimeDistributed(cnn_layer_object)(inputs)
    encoded_sequence = layers.Bidirectional(tf.keras.layers.LSTM(64, return_sequences = True))(encoded_frames)
    encoded_sequence = layers.Bidirectional(tf.keras.layers.LSTM(64, return_sequences = True))(encoded_sequence)

    hidden_layer = tf.keras.layers.Dense(units=64, activation="relu")(encoded_sequence)
    hidden_layer = layers.Dropout(0.4)(hidden_layer)
    hidden_layer = tf.keras.layers.Dense(units=64, activation="relu")(hidden_layer)
    hidden_layer = layers.Dropout(0.4)(hidden_layer)

    outputs = tf.keras.layers.Dense(units=2, activation="sigmoid")(hidden_layer)
    model = tf.keras.Model(inputs, outputs)
    model.compile(loss=tf.keras.losses.binary_crossentropy, optimizer=tf.optimizers.Adam(),metrics=[keras.metrics.TruePositives(name='tp'),
              keras.metrics.FalsePositives(name='fp'),
              keras.metrics.TrueNegatives(name='tn'),
              keras.metrics.FalseNegatives(name='fn'), 
              keras.metrics.BinaryAccuracy(name='accuracy'),
              keras.metrics.Precision(name='precision'),
              keras.metrics.Recall(name='recall')])
    model.summary()
    return model





def chioce_top18(tensor, geno = False):
    batch_size, window_size, rowcount = tensor.shape[0], tensor.shape[1], tensor.shape[2]
    if(rowcount < 32):
        tensor = np.concatenate([tensor, np.zeros([batch_size, window_size, 32-rowcount], dtype = np.float32)], axis = 2)
    if(geno == True):
        return tensor[:,:,:32].reshape(tensor.shape[0], tensor.shape[1], 32, 1)
    index = np.argsort(tensor.sum(axis = 1, keepdims = True), axis = 2)[:,:,-32:].reshape(-1, 32)
    l = []
    for loc in range(batch_size):
        l.append(tensor[loc][:,index[loc]])
    return np.concatenate(l, axis = 0).reshape(tensor.shape[0], tensor.shape[1], 32, 1)
###########################
def insloc(rs, start, end, cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist,  row):
    loc = rs
    rec = False
    for item in cigartuples:
        if(item[0] not in [1, 4, 5, 6]):
            if(rec == False and loc>=start):
                rec = True
            if(rec == True and item[0] == 2):
                for bias in range(min(item[1], end - loc)):
                    del_collist.append(loc+bias)
                    del_rowlist.append(row)

                    
            loc += item[1]
            
            
            if(loc>=end):
                break
        elif(rec == True and item[0] == 1):
            ins_collist.append(loc)
            ins_rowlist.append(row)
            ins_data.append(item[1])
############################

def baseinfo_AlignedSegment(bamfilepath, contig, start, end, window_size, maxread, outputpath, samplelocation = np.array(0)):


    if(samplelocation.size == 1):
        samplelocation = start + np.column_stack((np.arange(0, window_size * (int((end - start - 1) / window_size) + 1), window_size).reshape((int(( end - start - 1) / window_size) + 1), 1), np.arange(0, window_size * (int((end - start - 1) / window_size) + 1), window_size).reshape((int((end - start - 1) / window_size) + 1), 1) + window_size))
        end = samplelocation[-1, 1]  
    ######################################

    collist = []
    rowlist = []
    del_collist = []
    del_rowlist = []
    del_data = []
    rowend = []
    ins_collist = []
    ins_rowlist = []
    ins_data = []
    count = 0
    maxend = 0
    minstart = 999999999
    rowcount = 0

    maxlen = 0
    nooverlap = False
    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)
    for AlignedSegment in bamfile.fetch(contig, start, end):



        nooverlap = True
        count += 1
        
        rs = AlignedSegment.reference_start

        cstart = max(rs, start)
        cend = min(AlignedSegment.reference_end, end)
        
        if(maxend<cend):
            maxend = cend
        if(minstart>cstart):
            minstart=cstart

        newrow = True
        loc = -1

        for oneend in rowend:
            loc += 1
            if(oneend<cstart):
                insloc(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist, loc)
                rowend[loc] = cend
                newrow = False
                break
        if(newrow == True and (rowcount <= 500)):
            rowcount += 1
            rowend.append(cend)
            insloc(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist, len(rowend)-1)
        
    geno = False
    t1, t2 = (coo_matrix((np.ones(len(del_rowlist)), (del_rowlist, del_collist)), shape=(len(rowend), end-start)).toarray()).T, coo_matrix((np.log(ins_data), (ins_rowlist, ins_collist)), shape=((len(rowend), end-start))).toarray().T
    data = np.concatenate((chioce_top18(t1.reshape(-1, window_size, t1.shape[1]), geno), chioce_top18(t2.reshape(-1, window_size, t2.shape[1]), geno)), axis = -1)
    
    
    
    ##################################################

    mask = data.reshape(((end-start)//window_size, window_size * maxread*2)).sum(axis =1) != 0
    data = data[mask]
    index = np.arange(start, end, window_size)[mask]


    print(data.shape, index.shape)
    if(outputpath[-1] != '/'):
        outputpath += '/'
    np.savez_compressed(outputpath+contig+','+str(start)+','+str(end)+',data', data = data, label = np.array([0]), index = index)
    print(contig, str(start), str(end), 'completed')
def call_deletion(svlist, contig, index, data, predict, genotype = False, genotypelist = [], window_size = 200):
    if(genotype == False):
        loc = -1
        insv = False
        for p in predict:
            loc += 1
            if(p > 0):
                if(insv == False):
                    svstartloc = index[loc]
                    qc = p
                    insv = True

                qc = max(qc, p)
            else:
                if(insv == True):

                    svlist.append(['DEL', contig, svstartloc, index[loc], qc, ''])
                    insv = False


                    continue


    else:
        locingenotype = 0
        genotypecache = []
        loc = -1
        insv = False
        for p in predict:
            loc += 1
            if(p > 0):
                if(insv == False):
                    svstartloc = index[loc]
                    qc = p
                    genotypecache = []
                    insv = True

                qc = max(qc, p)
                genotypecache.append(np.argsort(genotypelist[locingenotype])[::-1][0])
                locingenotype += 1
            else:
                if(insv == True):

                    svlist.append(['DEL', contig, svstartloc, index[loc], qc, genotypecache])
                    genotypecache = []
                    insv = False

  
                    continue



        
    return svlist

def call_Insertion(svlist, contig, index, data, predict, genotype = False, genotypelist = [],  window_size = 200):

    if(genotype == False):
        loc = 0
        insv = False
        for p in predict:
            if(p > 0):
                if(insv == False):
                    #svstartloc = index[loc]
                    qc = p
                    insv = True
                    svstartloc = index[loc]
                qc = max(qc, p)
                if(qc > p):
                    svstartloc = index[loc]
            else:
                if(insv == True):
                    svlist.append(['INS', contig, svstartloc, index[loc], qc, ''])
                    insv = False
            loc += 1  
    else:
        
        locingenotype = 0
        genotypecache = []
        
        loc = 0
        insv = False
        for p in predict:
            if(p > 0):
                if(insv == False):
                    #svstartloc = index[loc]
                    qc = p
                    insv = True
                    svstartloc = index[loc]
                    genotypecache = []
                qc = max(qc, p)
                if(qc > p):
                    svstartloc = index[loc]
                genotypecache.append(np.argsort(genotypelist[locingenotype])[::-1][0])
                locingenotype += 1
            else:
                if(insv == True):
                    svlist.append(['INS', contig, svstartloc, index[loc], qc, genotypecache])
                    genotypecache = []
                    insv = False
            loc += 1
        
        
    return svlist

def call_genotype(model, prediction, data):
    #model.load_weights(genotypeweight)
    mask = np.where(prediction>0)[0]

    tmpdata = tf.data.Dataset.from_tensor_slices(data).batch(8048)
    return model.predict(tmpdata)



def batchdata(data, timestep, window_size, maxread):
    numoftimestep = ((data.shape[0])//timestep) - 1
    tailtimestep = data.shape[0] - numoftimestep*timestep

    return data[:numoftimestep*timestep].reshape((numoftimestep, timestep, window_size, maxread, data.shape[-1])), data[numoftimestep*timestep:].reshape((1, tailtimestep, window_size, maxread, data.shape[-1]))

  
def feature_matrix(bamfilepath, outputpath = './data/', max_worker = 16, vcfpath = '', window_size = 200, maxread = 32, block = 5e6, includecontig = [], genotype = False):

    if __name__ == '__main__':
        max_worker = min(max_worker, (len(os.sched_getaffinity(0))-1))
        print('Process number limit to ', max_worker-1)
        print('maxread ', maxread)
        

        bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)

        
        block = int(block)
        contig2length = {}
        for count in range(len(bamfile.get_index_statistics())):
            contig2length[bamfile.get_index_statistics()[count].contig] = bamfile.lengths[count]
        if(includecontig == []):
            for contig in contig2length:
                includecontig.append(contig)
        for contig in includecontig:
            for AlignedSegment in bamfile.fetch(contig):
                start = AlignedSegment.reference_start
                break
            while((start + 2*block) <  contig2length[contig]):
                while(True):
                    if(len(multiprocessing.active_children()) < max_worker):
                        print('working on contig = ', contig, start, start + block)
                        contig = str(contig)                                          
                        multiprocessing.Process(target=baseinfo_AlignedSegment, args=(bamfilepath, contig,start, start + block, window_size, maxread, outputpath)).start()
                        break
                    else:
                        time.sleep(2)
                start +=  block
            while(True):
                if(len(multiprocessing.active_children()) < max_worker):
                    print('working on contig = ', contig, start, start + block)
                    contig = str(contig)
                    multiprocessing.Process(target=baseinfo_AlignedSegment, args=( bamfilepath, contig,start, contig2length[contig], window_size, maxread, outputpath)).start()
                    break
                else:
                    time.sleep(2)


import time
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Input, Dense, Conv2D, DepthwiseConv2D, BatchNormalization, Dropout, GlobalAveragePooling2D, Reshape, multiply, add, Activation
import multiprocessing
import os
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Value, Array
import time
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
import tensorflow as tf
import threading
from os import listdir
from os.path import isfile, join
#tf.config.set_visible_devices([], 'GPU')

tf.random.set_seed(123)


class convse(keras.layers.Layer):
    
    def __init__(self,filters, kernel_size, strides=(1, 1), padding='valid', activation='relu', se_ratio = 0.25):

        super(convse, self).__init__()
        self.conv = layers.Conv2D(filters, kernel_size, strides=strides, padding=padding, activation=activation)
        
        self.filters = filters
        self.kernel_size = kernel_size
        self.strides = strides
        
        self.filters_se = max(1, int(filters*se_ratio))
        self.conv_1 = Conv2D(filters=self.filters_se,
                    kernel_size=1,
                    padding='same',
                    use_bias=True)
        self.conv_2 = Conv2D(filters=self.filters,
                    kernel_size=1,
                    padding='same',
                    activation='sigmoid',
                    use_bias=True)

    def call(self, inputs):
        
        x = self.conv(inputs)



        se = GlobalAveragePooling2D()(x)
        se = Reshape((1, 1, self.filters))(se)


        se = self.conv_1(se)

        se = Activation(tf.nn.swish)(se)

        se = self.conv_2(se)


        return multiply([x, se])
    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[1], 1, self.filters)



def cnn_layer():
    model = keras.models.Sequential()

    model.add(tf.keras.layers.AveragePooling2D([2, 4]))
    model.add(convse(108, (1, 8), strides=(1, 8), padding='valid', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(64, (2, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(64, (3, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(80, (4, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([3,1]))

    model.add(convse(80, (3, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(64, (2, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))
    model.add(layers.Flatten())


    return model

def init_model():
    inputs = tf.keras.Input((None, 200, 32, 2))

    cnn_layer_object = cnn_layer()
    encoded_frames = tf.keras.layers.TimeDistributed(cnn_layer_object)(inputs)
    encoded_sequence = layers.Bidirectional(tf.keras.layers.LSTM(64, return_sequences = True))(encoded_frames)
    encoded_sequence = layers.Bidirectional(tf.keras.layers.LSTM(64, return_sequences = True))(encoded_sequence)

    hidden_layer = tf.keras.layers.Dense(units=64, activation="relu")(encoded_sequence)
    hidden_layer = layers.Dropout(0.4)(hidden_layer)
    hidden_layer = tf.keras.layers.Dense(units=64, activation="relu")(hidden_layer)
    hidden_layer = layers.Dropout(0.4)(hidden_layer)

    outputs = tf.keras.layers.Dense(units=2, activation="sigmoid")(hidden_layer)
    model = tf.keras.Model(inputs, outputs)
    model.compile(loss=tf.keras.losses.binary_crossentropy, optimizer=tf.optimizers.Adam(),metrics=[keras.metrics.TruePositives(name='tp'),
              keras.metrics.FalsePositives(name='fp'),
              keras.metrics.TrueNegatives(name='tn'),
              keras.metrics.FalseNegatives(name='fn'), 
              keras.metrics.BinaryAccuracy(name='accuracy'),
              keras.metrics.Precision(name='precision'),
              keras.metrics.Recall(name='recall')])
    model.summary()
    return model





def chioce_top18(tensor, geno = False):
    batch_size, window_size, rowcount = tensor.shape[0], tensor.shape[1], tensor.shape[2]
    if(rowcount < 32):
        tensor = np.concatenate([tensor, np.zeros([batch_size, window_size, 32-rowcount], dtype = np.float32)], axis = 2)
    if(geno == True):
        return tensor[:,:,:32].reshape(tensor.shape[0], tensor.shape[1], 32, 1)
    index = np.argsort(tensor.sum(axis = 1, keepdims = True), axis = 2)[:,:,-32:].reshape(-1, 32)
    l = []
    for loc in range(batch_size):
        l.append(tensor[loc][:,index[loc]])
    return np.concatenate(l, axis = 0).reshape(tensor.shape[0], tensor.shape[1], 32, 1)
###########################
def insloc_1(rs, start, end, cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist,  row):
    loc = rs
    rec = False
    for item in cigartuples:
        if(item[0] not in [1, 4, 5, 6]):
            if(rec == False and loc>=start):
                rec = True
            if(rec == True and item[0] == 2):
                for bias in range(min(item[1], end - loc)):
                    del_collist.append(loc+bias)
                    del_rowlist.append(row)

                    
            loc += item[1]
            
            
            if(loc>=end):
                break
        elif(rec == True and item[0] == 1):
            ins_collist.append(loc)
            ins_rowlist.append(row)
            ins_data.append(item[1])
############################

def baseinfo_AlignedSegment(bamfilepath, contig, start, end, window_size, maxread, outputpath, samplelocation = np.array(0)):


    if(samplelocation.size == 1):
        samplelocation = start + np.column_stack((np.arange(0, window_size * (int((end - start - 1) / window_size) + 1), window_size).reshape((int(( end - start - 1) / window_size) + 1), 1), np.arange(0, window_size * (int((end - start - 1) / window_size) + 1), window_size).reshape((int((end - start - 1) / window_size) + 1), 1) + window_size))
        end = samplelocation[-1, 1]  
    ######################################

    collist = []
    rowlist = []
    del_collist = []
    del_rowlist = []
    del_data = []
    rowend = []
    ins_collist = []
    ins_rowlist = []
    ins_data = []
    count = 0
    maxend = 0
    minstart = 999999999
    rowcount = 0

    maxlen = 0
    nooverlap = False
    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)
    for AlignedSegment in bamfile.fetch(contig, start, end):



        nooverlap = True
        count += 1
        
        rs = AlignedSegment.reference_start

        cstart = max(rs, start)
        cend = min(AlignedSegment.reference_end, end)
        
        if(maxend<cend):
            maxend = cend
        if(minstart>cstart):
            minstart=cstart

        newrow = True
        loc = -1

        for oneend in rowend:
            loc += 1
            if(oneend<cstart):
                insloc_1(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist, loc)
                rowend[loc] = cend
                newrow = False
                break
        if(newrow == True and (rowcount <= 500)):
            rowcount += 1
            rowend.append(cend)
            insloc_1(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist, len(rowend)-1)
        
    geno = False
    t1, t2 = (coo_matrix((np.ones(len(del_rowlist)), (del_rowlist, del_collist)), shape=(len(rowend), end-start)).toarray()).T, coo_matrix((np.log(ins_data), (ins_rowlist, ins_collist)), shape=((len(rowend), end-start))).toarray().T
    data = np.concatenate((chioce_top18(t1.reshape(-1, window_size, t1.shape[1]), geno), chioce_top18(t2.reshape(-1, window_size, t2.shape[1]), geno)), axis = -1)
    
    
    
    ##################################################

    mask = data.reshape(((end-start)//window_size, window_size * maxread*2)).sum(axis =1) != 0
    data = data[mask]
    index = np.arange(start, end, window_size)[mask]


    print(data.shape, index.shape)
    np.savez_compressed(outputpath+contig+','+str(start)+','+str(end)+',data', data = data, label = np.array([0]), index = index)
    print(contig, str(start), str(end), 'completed')
def call_deletion(svlist, contig, index, data, predict, genotype = False, genotypelist = [], window_size = 200):
    if(genotype == False):
        loc = -1
        insv = False
        for p in predict:
            loc += 1
            if(p > 0):
                if(insv == False):
                    svstartloc = index[loc]
                    qc = p
                    insv = True

                qc = max(qc, p)
            else:
                if(insv == True):

                    svlist.append(['DEL', contig, svstartloc, index[loc], qc, ''])
                    insv = False


                    continue


    else:
        locingenotype = 0
        genotypecache = []
        loc = -1
        insv = False
        for p in predict:
            loc += 1
            if(p > 0):
                if(insv == False):
                    svstartloc = index[loc]
                    qc = p
                    genotypecache = []
                    insv = True

                qc = max(qc, p)
                genotypecache.append(np.argsort(genotypelist[locingenotype])[::-1][0])
                locingenotype += 1
            else:
                if(insv == True):

                    svlist.append(['DEL', contig, svstartloc, index[loc], qc, genotypecache])
                    genotypecache = []
                    insv = False

  
                    continue



        
    return svlist

def call_Insertion(svlist, contig, index, data, predict, genotype = False, genotypelist = [],  window_size = 200):

    if(genotype == False):
        loc = 0
        insv = False
        for p in predict:
            if(p > 0):
                if(insv == False):
                    #svstartloc = index[loc]
                    qc = p
                    insv = True
                    svstartloc = index[loc]
                qc = max(qc, p)
                if(qc > p):
                    svstartloc = index[loc]
            else:
                if(insv == True):
                    svlist.append(['INS', contig, svstartloc, index[loc], qc, ''])
                    insv = False
            loc += 1  
    else:
        
        locingenotype = 0
        genotypecache = []
        
        loc = 0
        insv = False
        for p in predict:
            if(p > 0):
                if(insv == False):
                    #svstartloc = index[loc]
                    qc = p
                    insv = True
                    svstartloc = index[loc]
                    genotypecache = []
                qc = max(qc, p)
                if(qc > p):
                    svstartloc = index[loc]
                genotypecache.append(np.argsort(genotypelist[locingenotype])[::-1][0])
                locingenotype += 1
            else:
                if(insv == True):
                    svlist.append(['INS', contig, svstartloc, index[loc], qc, genotypecache])
                    genotypecache = []
                    insv = False
            loc += 1
        
        
    return svlist

def call_genotype(model, prediction, data):
    #model.load_weights(genotypeweight)
    mask = np.where(prediction>0)[0]

    tmpdata = tf.data.Dataset.from_tensor_slices(data).batch(8048)
    return model.predict(tmpdata)



def batchdata(data, timestep, window_size, maxread):
    numoftimestep = ((data.shape[0])//timestep) - 1
    tailtimestep = data.shape[0] - numoftimestep*timestep

    return data[:numoftimestep*timestep].reshape((numoftimestep, timestep, window_size, maxread, data.shape[-1])), data[numoftimestep*timestep:].reshape((1, tailtimestep, window_size, maxread, data.shape[-1]))

  
def feature_matrix(bamfilepath, outputpath = './data/', max_worker = 16, vcfpath = '', window_size = 200, maxread = 32, block = 5e6, includecontig = [], genotype = False):

    if __name__ == '__main__':
        max_worker = min(max_worker, (len(os.sched_getaffinity(0))-1))
        print('Process number limit to ', max_worker-1)
        print('maxread ', maxread)
        
        if(outputpath[-1] != '/'):
            outputpath += '/'
        if(os.path.isdir(outputpath) == True):
            for path in [outputpath+f for f in listdir(outputpath) if(isfile(join(outputpath, f)) and 'npz' in f)]:
                os.remove(path)
        else:
            os.mkdir(outputpath[:-1])
            
        bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)

        
        block = int(block)
        contig2length = {}
        for count in range(len(bamfile.get_index_statistics())):
            contig2length[bamfile.get_index_statistics()[count].contig] = bamfile.lengths[count]
        if(includecontig == []):
            for contig in contig2length:
                includecontig.append(contig)
        for contig in includecontig:
            for AlignedSegment in bamfile.fetch(contig):
                start = AlignedSegment.reference_start
                break
            while((start + 2*block) <  contig2length[contig]):
                while(True):
                    if(len(multiprocessing.active_children()) < max_worker):
                        print('working on contig = ', contig, start, start + block)
                        contig = str(contig)                                          
                        multiprocessing.Process(target=baseinfo_AlignedSegment, args=(bamfilepath, contig,start, start + block, window_size, maxread, outputpath)).start()
                        break
                    else:
                        time.sleep(2)
                start +=  block
            while(True):
                if(len(multiprocessing.active_children()) < max_worker):
                    print('working on contig = ', contig, start, start + block)
                    contig = str(contig)
                    multiprocessing.Process(target=baseinfo_AlignedSegment, args=( bamfilepath, contig,start, contig2length[contig], window_size, maxread, outputpath)).start()
                    break
                else:
                    time.sleep(2)

                


import time
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Input, Dense, Conv2D, DepthwiseConv2D, BatchNormalization, Dropout, GlobalAveragePooling2D, Reshape, multiply, add, Activation
import multiprocessing
import os
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import Value, Array
import time
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
import matplotlib.pyplot as plt
import tensorflow as tf
import threading
from os import listdir
from os.path import isfile, join
#tf.config.set_visible_devices([], 'GPU')

tf.random.set_seed(123)


class convse(keras.layers.Layer):
    
    def __init__(self,filters, kernel_size, strides=(1, 1), padding='valid', activation='relu', se_ratio = 0.25):

        super(convse, self).__init__()
        self.conv = layers.Conv2D(filters, kernel_size, strides=strides, padding=padding, activation=activation)
        
        self.filters = filters
        self.kernel_size = kernel_size
        self.strides = strides
        
        self.filters_se = max(1, int(filters*se_ratio))
        self.conv_1 = Conv2D(filters=self.filters_se,
                    kernel_size=1,
                    padding='same',
                    use_bias=True)
        self.conv_2 = Conv2D(filters=self.filters,
                    kernel_size=1,
                    padding='same',
                    activation='sigmoid',
                    use_bias=True)

    def call(self, inputs):
        
        x = self.conv(inputs)



        se = GlobalAveragePooling2D()(x)
        se = Reshape((1, 1, self.filters))(se)


        se = self.conv_1(se)

        se = Activation(tf.nn.swish)(se)

        se = self.conv_2(se)


        return multiply([x, se])
    def compute_output_shape(self, input_shape):
        return (input_shape[0], input_shape[1], 1, self.filters)



def cnn_layer():
    model = keras.models.Sequential()

    model.add(tf.keras.layers.AveragePooling2D([2, 4]))
    model.add(convse(108, (1, 8), strides=(1, 8), padding='valid', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(64, (2, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(64, (3, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(80, (4, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([3,1]))

    model.add(convse(80, (3, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(64, (2, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))
    model.add(layers.Flatten())


    return model

def init_model():
    inputs = tf.keras.Input((None, 200, 32, 2))

    cnn_layer_object = cnn_layer()
    encoded_frames = tf.keras.layers.TimeDistributed(cnn_layer_object)(inputs)
    encoded_sequence = layers.Bidirectional(tf.keras.layers.LSTM(64, return_sequences = True))(encoded_frames)
    encoded_sequence = layers.Bidirectional(tf.keras.layers.LSTM(64, return_sequences = True))(encoded_sequence)

    hidden_layer = tf.keras.layers.Dense(units=64, activation="relu")(encoded_sequence)
    hidden_layer = layers.Dropout(0.4)(hidden_layer)
    hidden_layer = tf.keras.layers.Dense(units=64, activation="relu")(hidden_layer)
    hidden_layer = layers.Dropout(0.4)(hidden_layer)

    outputs = tf.keras.layers.Dense(units=2, activation="sigmoid")(hidden_layer)
    model = tf.keras.Model(inputs, outputs)
    model.compile(loss=tf.keras.losses.binary_crossentropy, optimizer=tf.optimizers.Adam(),metrics=[keras.metrics.TruePositives(name='tp'),
              keras.metrics.FalsePositives(name='fp'),
              keras.metrics.TrueNegatives(name='tn'),
              keras.metrics.FalseNegatives(name='fn'), 
              keras.metrics.BinaryAccuracy(name='accuracy'),
              keras.metrics.Precision(name='precision'),
              keras.metrics.Recall(name='recall')])
    model.summary()
    return model





def chioce_top18(tensor, geno = False):
    batch_size, window_size, rowcount = tensor.shape[0], tensor.shape[1], tensor.shape[2]
    k = 18
    if(rowcount < k):
        tensor = np.concatenate([tensor, np.zeros([batch_size, window_size, k-rowcount], dtype = np.float32)], axis = 2)
    if(geno == True):
        return tensor[:,:,:k].reshape(tensor.shape[0], tensor.shape[1], k, 1)
    index = np.argsort(tensor.sum(axis = 1, keepdims = True), axis = 2)[:,:,-k:].reshape(-1, k)
    l = []
    for loc in range(batch_size):
        l.append(tensor[loc][:,index[loc]])
    return np.concatenate(l, axis = 0).reshape(tensor.shape[0], tensor.shape[1], k, 1)
###########################
def insloc_1(rs, start, end, cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist,  row):
    loc = rs
    rec = False
    for item in cigartuples:
        if(item[0] not in [1, 4, 5, 6]):
            if(rec == False and loc>=start):
                rec = True
            if(rec == True and item[0] == 2):
                for bias in range(min(item[1], end - loc)):
                    del_collist.append(loc+bias)
                    del_rowlist.append(row)

                    
            loc += item[1]
            
            
            if(loc>=end):
                break

############################

def baseinfo_AlignedSegment(bamfilepath, contig, start, end, window_size, maxread, outputpath, samplelocation = np.array(0)):


    if(samplelocation.size == 1):
        samplelocation = start + np.column_stack((np.arange(0, window_size * (int((end - start - 1) / window_size) + 1), window_size).reshape((int(( end - start - 1) / window_size) + 1), 1), np.arange(0, window_size * (int((end - start - 1) / window_size) + 1), window_size).reshape((int((end - start - 1) / window_size) + 1), 1) + window_size))
        end = samplelocation[-1, 1]  
    ######################################

    collist = []
    rowlist = []
    del_collist = []
    del_rowlist = []
    del_data = []
    rowend = []
    ins_collist = []
    ins_rowlist = []
    ins_data = []
    count = 0
    maxend = 0
    minstart = 999999999
    rowcount = 0

    maxlen = 0
    nooverlap = False
    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)
    blocksize = int(1e6)
    data = []

    for AlignedSegment in bamfile.fetch(contig, start, end):



        nooverlap = True
        count += 1

        rs = AlignedSegment.reference_start

        cstart = max(rs, start)
        cend = min(AlignedSegment.reference_end, end)

        if(maxend<cend):
            maxend = cend
        if(minstart>cstart):
            minstart=cstart

        newrow = True
        loc = -1

        for oneend in rowend:
            loc += 1
            if(oneend<cstart):
                insloc_1(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist, loc)
                rowend[loc] = cend
                newrow = False
                break
        if(newrow == True and (rowcount <= 70)):
            rowcount += 1
            rowend.append(cend)
            insloc_1(rs-minstart, start-minstart, end-minstart, AlignedSegment.cigartuples, ins_collist, ins_rowlist, ins_data, del_collist, del_rowlist, len(rowend)-1)

    geno = False
    t1 = (coo_matrix((np.ones(len(del_rowlist)), (del_rowlist, del_collist)), shape=(len(rowend), end-start)).toarray()).T
    data = (chioce_top18(t1.reshape(-1, window_size, t1.shape[1]), geno))
    
    
    ##################################################

    mask = data.reshape(((end-start)//window_size, window_size * maxread)).sum(axis =1) != 0
    data = data[mask]
    index = np.arange(start, end, window_size)[mask]


    np.savez_compressed(outputpath+contig+','+str(start)+','+str(end)+',data', data = data, label = np.array([0]), index = index)
    print(contig, str(start), str(end), 'completed')
def call_deletion(svlist, contig, index, data, predict, genotype = False, genotypelist = [], window_size = 200):
    if(genotype == False):
        loc = -1
        insv = False
        for p in predict:
            loc += 1
            if(p > 0):
                if(insv == False):
                    svstartloc = index[loc]
                    qc = p
                    insv = True

                qc = max(qc, p)
            else:
                if(insv == True):

                    svlist.append(['DEL', contig, svstartloc, index[loc], qc, ''])
                    insv = False


                    continue


    else:
        locingenotype = 0
        genotypecache = []
        loc = -1
        insv = False
        for p in predict:
            loc += 1
            if(p > 0):
                if(insv == False):
                    svstartloc = index[loc]
                    qc = p
                    genotypecache = []
                    insv = True

                qc = max(qc, p)
                genotypecache.append(np.argsort(genotypelist[locingenotype])[::-1][0])
                locingenotype += 1
            else:
                if(insv == True):

                    svlist.append(['DEL', contig, svstartloc, index[loc], qc, genotypecache])
                    genotypecache = []
                    insv = False

  
                    continue



        
    return svlist

def call_Insertion(svlist, contig, index, data, predict, genotype = False, genotypelist = [],  window_size = 200):

    if(genotype == False):
        loc = 0
        insv = False
        for p in predict:
            if(p > 0):
                if(insv == False):
                    #svstartloc = index[loc]
                    qc = p
                    insv = True
                    svstartloc = index[loc]
                qc = max(qc, p)
                if(qc > p):
                    svstartloc = index[loc]
            else:
                if(insv == True):
                    svlist.append(['INS', contig, svstartloc, index[loc], qc, ''])
                    insv = False
            loc += 1  
    else:
        
        locingenotype = 0
        genotypecache = []
        
        loc = 0
        insv = False
        for p in predict:
            if(p > 0):
                if(insv == False):
                    #svstartloc = index[loc]
                    qc = p
                    insv = True
                    svstartloc = index[loc]
                    genotypecache = []
                qc = max(qc, p)
                if(qc > p):
                    svstartloc = index[loc]
                genotypecache.append(np.argsort(genotypelist[locingenotype])[::-1][0])
                locingenotype += 1
            else:
                if(insv == True):
                    svlist.append(['INS', contig, svstartloc, index[loc], qc, genotypecache])
                    genotypecache = []
                    insv = False
            loc += 1
        
        
    return svlist

def call_genotype(model, prediction, data):
    #model.load_weights(genotypeweight)
    mask = np.where(prediction>0)[0]

    tmpdata = tf.data.Dataset.from_tensor_slices(data).batch(8048)
    return model.predict(tmpdata)



def batchdata(data, timestep, window_size, maxread):
    numoftimestep = ((data.shape[0])//timestep) - 1
    tailtimestep = data.shape[0] - numoftimestep*timestep

    return data[:numoftimestep*timestep].reshape((numoftimestep, timestep, window_size, maxread, data.shape[-1])), data[numoftimestep*timestep:].reshape((1, tailtimestep, window_size, maxread, data.shape[-1]))

  
def feature_matrix(bamfilepath, outputpath = './data/', max_worker = 8, vcfpath = '', window_size = 200, maxread = 18, block = 1e7, includecontig = [], genotype = False):

    max_worker = min(max_worker, (len(os.sched_getaffinity(0))-1))
    print('Process number limit to ', max_worker-1)

    if(outputpath[-1] != '/'):
        outputpath += '/'
    if(os.path.isdir(outputpath) == True):
        for path in [outputpath+f for f in listdir(outputpath) if(isfile(join(outputpath, f)) and 'npz' in f)]:
            os.remove(path)
    else:
        os.mkdir(outputpath[:-1])

    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)


    block = int(block)
    contig2length = {}
    for count in range(len(bamfile.get_index_statistics())):
        contig2length[bamfile.get_index_statistics()[count].contig] = bamfile.lengths[count]
    if(includecontig == []):
        for contig in contig2length:
            includecontig.append(contig)
    for contig in includecontig:
        start = 0
        b = False
        no = True
        while(True):
            for AlignedSegment in bamfile.fetch(contig, start, start + block):
                no = False
                if(AlignedSegment.reference_start<start):
                    b = True
                    continue
                if(b != True):
                    start = AlignedSegment.reference_start
                break
            if((no == True) and (start < contig2length[contig])):
                start += block
            else:
                break

        for AlignedSegment in bamfile.fetch(contig, start, start + block):
            pass
        end = min(start + block, AlignedSegment.reference_end)
        while((start) <  contig2length[contig]):


            while(True):
                if(len(multiprocessing.active_children()) < max_worker):
                    print('working on contig = ', contig, start, start + block)
                    contig = str(contig)                                          
                    multiprocessing.Process(target=baseinfo_AlignedSegment, args=(bamfilepath, contig,start, start + block, window_size, maxread, outputpath)).start()
                    break
                else:
                    time.sleep(2)

            start +=  block
            b = False
            no = True
            while(True):
                for AlignedSegment in bamfile.fetch(contig, start, start + block):
                    no = False
                    if(AlignedSegment.reference_start<start):
                        b = True
                        continue
                    if(b != True):
                        start = AlignedSegment.reference_start
                    break
                if((no == True) and (start < contig2length[contig])):
                    start += block
                else:
                    break
            for AlignedSegment in bamfile.fetch(contig, start, start + block):
                pass
            end = min(start + block, AlignedSegment.reference_end)
                
        
                            
                    
                
                

       
                            
                    
                
                



                
        
                            
                    
                
                


            
    