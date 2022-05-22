import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Input, Dense, Conv2D, DepthwiseConv2D, BatchNormalization, Dropout, GlobalAveragePooling2D, Reshape, multiply, add, Activation
#tf.config.set_visible_devices([], 'GPU')
tf.random.set_seed(123)
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Input, Dense, Conv2D, DepthwiseConv2D, BatchNormalization, Dropout, GlobalAveragePooling2D, Reshape, multiply, add, Activation
from feature import *
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



import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras


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
from os import listdir
from os.path import isfile, join
import time
from collections import Counter
import numpy as np
import pysam
import time
import time
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
def cigarinfo(cigararray, refstartarray, start, end, reccount, cigarweight):

    a = (~((cigararray[:,0] == 1) | (cigararray[:,0] == 4))).reshape(cigararray.shape[0], 1) * cigararray
    a1 = a[:,1].reshape((reccount, cigarweight))
    a1 = (a1 @ np.triu(np.ones((a1.shape[1] , a1.shape[1])), k=0)) + refstartarray
    a2 = (a1 != 0) * (np.arange(a1.shape[0]).reshape(a1.shape[0], 1))
    a = np.concatenate([cigararray, a2.reshape((cigararray.shape[0], 1)), a1.reshape([cigararray.shape[0], 1])], axis = 1)
    

    return a[(start <= a[:,-1]) & (a[:,-1] < end)]



def baseinfo(bamfile, contig, start, end, SVtype):
    

    cigararray = []
    readstartandend = []
    refpositonlist = []
    refpositonweight = []
    substitionarray, deletionarray, substitionweight, deletionweight = [], [], [], []
    nooverlap = True
    qualityarray = []

    qualityarray = []

    for AlignedSegment in bamfile.fetch(contig, start, end):



        nooverlap = False
        cigararray.append(tf.keras.backend.flatten(tf.constant(AlignedSegment.cigartuples)))
        readstartandend.append([AlignedSegment.reference_start-start, AlignedSegment.reference_end-start, AlignedSegment.mapping_quality, (1 - (AlignedSegment.query_alignment_length / AlignedSegment.infer_read_length()))**2])


    if(nooverlap == True):
        return []
    readstartandend = np.array(readstartandend)
    cigararray = tf.keras.preprocessing.sequence.pad_sequences(cigararray)
    reccount, cigarweight = cigararray.shape[0], int(cigararray.shape[1] / 2)
    cigararray = cigararray.reshape(int(cigararray.size / 2), 2)

    cigararray = cigarinfo(cigararray, readstartandend[:,:1], 0, end - start, reccount, cigarweight)
    if(SVtype == 'DEL'):
        a = cigararray[(cigararray[:,0] == 2) & (cigararray[:,1] > 20)]
    else:
        a = cigararray[(cigararray[:,0] == 1) & (cigararray[:,1] > 20)]
    if(a.size == 0):
        return []
    a[:,-1] = a[:,-1] - a[:,-2] 
    delsig = np.column_stack((a[:,-1:], a[:,-3:-2]))

    loc = np.array(list(delsig))[:,0]
    binnum = 20
    binresult = (loc//binnum)
    mc = Counter(binresult).most_common(1)[0]
    sp = mc[1]
    minv, maxv = mc[0]-1, mc[0]+1
    tmp = np.median(np.array(list(delsig))[(minv<=binresult) *  (maxv>= binresult)], axis = 0).astype('int64')
    tmp[0] = tmp[0]+start
    return [contig] + tmp.tolist()+[sp, SVtype]
def insloc(rs, start, end, cigartuples, delinfo, insinfo):
    loc = rs
    rec = False
    insstart = -1
    inssize = 0
    insloc = 0
    delstart = -1
    delsize = 0
    for item in cigartuples:
        if(item[0] not in [1, 4, 5, 6]):
            if(rec == False and loc>=start):
                rec = True
            if(rec == True and item[0] == 2):
                if(item[1]>49):
                    if(delstart == -1):
                        delstart = loc
                    delsize = loc +item[1] - delstart 

                    
            loc += item[1]
            
            
            if(loc>=end):
                break
        elif(rec == True and item[0] == 1):
            if(item[1]>49):
                if(insstart == -1):
                    insstart = loc
                    insloc = loc
                inssize = item[1] + loc - insloc
                insloc = loc
    if(delstart != -1):
        delinfo.append([delstart, delsize])
    if(insstart != -1):
        insinfo.append([insstart, inssize])
def baseinfo(bamfile, contig, start, end, SVtype):
    delinfo = []
    insinfo = []
    for AlignedSegment in bamfile.fetch(contig, start, end):
        insloc(AlignedSegment.reference_start, start, end, AlignedSegment.cigartuples, delinfo, insinfo)
    if(SVtype == 'INS'):
        delsig = np.array(insinfo)
    else:
        delsig = np.array(delinfo)
    if(delsig.size == 0):
        return []

    loc = (delsig)[:,0]    
    binnum = 20
    binresult = (loc//binnum)
    mc = Counter(binresult).most_common(1)[0]
    #print(mc)
    #print(delsig)
    sp = mc[1]
    minv, maxv = mc[0]-1, mc[0]+1
    tmp = np.median(np.array(list(delsig))[(minv<=binresult) *  (maxv>= binresult)], axis = 0).astype('int64')
    return [contig] + tmp.tolist()+[sp, SVtype]                

    
def baseinfo_main_binsaver(bamfilepath, delloc):




    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)


    delsig = []
    for rec in delloc:
        contig, start, end = str(rec[1]), int(rec[2]), int(rec[3])
        bk = baseinfo(bamfile, contig, max(0, start), end, rec[0])
        if(len(bk) == 0):
            continue
        delsig.append(bk)

    return delsig

def tovcf(rawsvlist, contig2length, outputpath):
    top = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">\n"""
    body = ''
    for contig in contig2length:
        body += "##contig=<ID="+contig+",length="+str(int(contig2length[contig]))+">\n"
    tail = """##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, INS=Insertion">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=RE,Number=1,Type=Integer,Description="Number of read support this record">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"""

    myvcf = top+body+tail
    #          contig    start        svlen         re
    svlist = [[rec[0], int(rec[1]), int(rec[2]), int(rec[3]), '.', rec[4]] for rec in rawsvlist]
        
    for rec in pd.DataFrame(svlist).sort_values([0, 1]).values:


        contig = str(rec[0])


        geno = rec[4]
        if(rec[5] == 'DEL'):
            recinfo = 'SVLEN=' + str(int(rec[2]))+';SVTYPE=' + 'DEL'+';END='+str(int(rec[1])+abs(int(rec[2])))+';RE='+str(int(rec[3])+1)+'\tGT\t'+geno+'\n'
            myvcf += (contig +'\t'+ str(int(rec[1]))+'\t'+ '.'+'\t'+ '.'+'\t'+ '.'+'\t'+ str(int(rec[3])+1)+'\t'+ 'PASS'+'\t'+recinfo)
        elif(rec[5] == 'INS'):
            recinfo = 'SVLEN=' + str(int(rec[2]))+';SVTYPE=' + 'INS'+';END='+str(int(rec[1])+1)+';RE='+str(int(rec[3])+1)+'\tGT\t'+geno+'\n'
            myvcf += (contig +'\t'+ str(int(rec[1]))+'\t'+ '.'+'\t'+ '.'+'\t'+ '.'+'\t'+ str(int(rec[3])+1)+'\t'+ 'PASS'+'\t'+recinfo)


    with open(outputpath, "w") as f:
        f.write(myvcf)

def submergeSVcall(predict, contig, index, bamfilepath, SVtype, region, window_size = 200):
    #contig: contig name in bam file, predict: 1d array, startloc: first prediction location in reference, window_size: the sub-region of a feature matrix included
    svlist = []
    loc = -1
    insv = False
    for p in predict:
        loc += 1
        if(p > 0):
            if(insv == False):
                svstartloc = index[loc]
                insv = True

        else:
            if(insv == True):

                svlist.append([SVtype, contig, svstartloc, index[loc], ])
                insv = False


                continue
    if(insv == True):
        svlist.append([SVtype, contig, svstartloc, index[count]])
    print('svlist size = ', len(svlist))
    print('predict sum = ', predict.sum())
    if(region != True):

        return baseinfo_main_binsaver(bamfilepath, svlist)
    else:
        return svlist


from collections import Counter
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

    model.add(tf.keras.layers.AveragePooling2D([1, 2]))
    model.add(convse(128, (1, 9), strides=(1, 9), padding='valid', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(80, (2, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(80, (3, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(80, (4, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([3,1]))

    model.add(convse(80, (3, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))

    model.add(convse(80, (2, 1), padding='same', activation='relu'))
    model.add(layers.MaxPooling2D([2,1]))
    model.add(layers.Flatten())


    return model
def init_model():
    inputs = tf.keras.Input((None, 200, 18, 1))

    cnn_layer_object = cnn_layer()
    encoded_frames = tf.keras.layers.TimeDistributed(cnn_layer_object)(inputs)
    encoded_sequence = layers.Bidirectional(tf.keras.layers.LSTM(64, return_sequences = True))(encoded_frames)
    encoded_sequence = layers.Bidirectional(tf.keras.layers.LSTM(64, return_sequences = True))(encoded_sequence)

    hidden_layer = tf.keras.layers.Dense(units=64, activation="relu")(encoded_sequence)
    hidden_layer = layers.Dropout(0.4)(hidden_layer)
    hidden_layer = tf.keras.layers.Dense(units=64, activation="relu")(hidden_layer)
    hidden_layer = layers.Dropout(0.4)(hidden_layer)

    outputs = tf.keras.layers.Dense(units=1, activation="sigmoid")(hidden_layer)
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
def insloc_2(rs, start, end, cigartuples, delinfo, insinfo):
    loc = rs
    rec = False
    insstart = -1
    inssize = 0
    insloc = 0
    delstart = -1
    delsize = 0
    for item in cigartuples:
        if(item[0] not in [1, 4, 5, 6]):
            if(rec == False and loc>=start):
                rec = True
            if(rec == True and item[0] == 2):
                if(item[1]>49):
                    if(delstart == -1):
                        delstart = loc
                    delsize = loc +item[1] - delstart 

                    
            loc += item[1]
            
            
            if(loc>=end):
                break
        elif(rec == True and item[0] == 1):
            if(item[1]>49):
                if(insstart == -1):
                    insstart = loc
                    insloc = loc
                inssize = item[1] + loc - insloc
                insloc = loc
    if(delstart != -1):
        delinfo.append([delstart, delsize])
    if(insstart != -1):
        insinfo.append([insstart, inssize])
def baseinfo(bamfile, contig, start, end, SVtype):
    delinfo = []
    insinfo = []
    for AlignedSegment in bamfile.fetch(contig, start, end):
        insloc_2(AlignedSegment.reference_start, start, end, AlignedSegment.cigartuples, delinfo, insinfo)
    if(SVtype == 'INS'):
        delsig = np.array(insinfo)
    else:
        delsig = np.array(delinfo)
    if(delsig.size == 0):
        return []

    loc = (delsig)[:,0]    
    binnum = 20
    binresult = (loc//binnum)
    mc = Counter(binresult).most_common(1)[0]
    #print(mc)
    #print(delsig)
    sp = mc[1]
    minv, maxv = mc[0]-1, mc[0]+1
    tmp = np.median(np.array(list(delsig))[(minv<=binresult) *  (maxv>= binresult)], axis = 0).astype('int64')
    return [contig] + tmp.tolist()+[sp, SVtype]                

def c_pos(cigar, refstart):
    number = ''
    numlist = [str(i) for i in range(10)]
    readstart = False
    readend = False
    refend = False
    readloc = 0
    refloc = refstart
    for c in cigar:
        if(c in numlist):
            number += c
        else:
            number = int(number)
            if(readstart == False and c in ['M', 'I', '=', 'X']):
                readstart = readloc
            if(readstart != False and c in ['H', 'S']):
                readend = readloc
                refend = refloc
                break

            if(c in ['M', 'I', 'S', '=', 'X']):
                readloc += number

            if(c in ['M', 'D', 'N', '=', 'X']):
                refloc += number
            number = ''
    if(readend == False):
        readend = readloc
        refend = refloc

    return refstart, refend, readstart, readend 
def decode_flag(Flag):
    signal = {1 << 2: 0, 1 >> 1: 1, 1 << 4: 2, 1 << 11: 3, 1 << 4 | 1 << 11: 4}
    return signal[Flag] if(Flag in signal) else 0
def baseinfo(bamfile, contig, start, end, SVtype):
    delinfo = []
    insinfo = []
    for AlignedSegment in bamfile.fetch(contig, start, end):
        insloc_2(AlignedSegment.reference_start, start, end, AlignedSegment.cigartuples, delinfo, insinfo)
        if(AlignedSegment.has_tag('SA') == True):
            code = decode_flag(AlignedSegment.flag)
            sapresent = True
            rawsalist = AlignedSegment.get_tag('SA').split(';')
            for sa in rawsalist[:-1]:
                sainfo = sa.split(',')
                tmpcontig, tmprefstart, strand, cigar = sainfo[0], int(sainfo[1]), sainfo[2], sainfo[3]
                if(tmpcontig != contig):
                    continue
                if((strand == '-' and (code %2) ==0) or (strand == '+' and (code %2) ==1)):
                    refstart_1, refend_1, readstart_1, readend_1 =  AlignedSegment.reference_start, AlignedSegment.reference_end, AlignedSegment.query_alignment_start, AlignedSegment.query_alignment_end
                    refstart_2, refend_2, readstart_2, readend_2 = c_pos(cigar, tmprefstart)
                    a = readend_1 - readstart_2
                    b = refend_1 - refstart_2
                    if(abs(a)<2000):
                        if(abs(b-a)<50):
                            continue
                        if((b-a)<0):
                            if(max(refend_1, refstart_2)<end and max(refend_1, refstart_2)<start):
                                delinfo.append([max(refend_1, refstart_2), abs(b-a)])
                        else:
                            if(refend_1<end and refend_1<start):
                                insinfo.append([refend_1, (b-a)])
                            if(refstart_2<end and refstart_2<start):
                                insinfo.append([refstart_2, (b-a)])
    if(SVtype == 'INS'):
        delsig = np.array(insinfo)
    else:
        delsig = np.array(delinfo)
    if(delsig.size == 0):
        return []

    loc = (delsig)[:,0]    
    binnum = 20
    binresult = (loc//binnum)
    mc = Counter(binresult).most_common(1)[0]
    #print(mc)
    #print(delsig)
    sp = mc[1]
    minv, maxv = mc[0]-1, mc[0]+1
    tmp = np.median(np.array(list(delsig))[(minv<=binresult) *  (maxv>= binresult)], axis = 0).astype('int64')
    return [contig] + tmp.tolist()+[sp, SVtype]    
def baseinfo_main_binsaver(bamfilepath, delloc, count, tmppath):




    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)


    delsig = []
    for rec in delloc:
        contig, start, end = str(rec[1]), int(rec[2]), int(rec[3])

        bk = baseinfo(bamfile, contig, max(0, start), end, rec[0])
        if(len(bk) == 0):
            continue
        delsig.append(bk)
    np.save(tmppath+str(count), np.array(delsig))
    return delsig
def tovcf(rawsvlist, contig2length, outputpath):
    top = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">\n"""
    body = ''
    for contig in contig2length:
        body += "##contig=<ID="+contig+",length="+str(int(contig2length[contig]))+">\n"
    tail = """##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV:DEL=Deletion, INS=Insertion">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=RE,Number=1,Type=Integer,Description="Number of read support this record">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"""

    myvcf = top+body+tail
    #          contig    start        svlen         re
    svlist = [[rec[0], int(rec[1]), int(rec[2]), int(rec[3]), '.', rec[4]] for rec in rawsvlist]
        
    for rec in pd.DataFrame(svlist).sort_values([0, 1]).values:


        contig = str(rec[0])


        geno = rec[4]
        if(rec[5] == 'DEL'):
            recinfo = 'SVLEN=' + str(int(rec[2]))+';SVTYPE=' + 'DEL'+';END='+str(int(rec[1])+abs(int(rec[2])))+';RE='+str(int(rec[3])+1)+'\tGT\t'+geno+'\n'
            myvcf += (contig +'\t'+ str(int(rec[1]))+'\t'+ '.'+'\t'+ '.'+'\t'+ '.'+'\t'+ str(int(rec[3])+1)+'\t'+ 'PASS'+'\t'+recinfo)
        elif(rec[5] == 'INS'):
            recinfo = 'SVLEN=' + str(int(rec[2]))+';SVTYPE=' + 'INS'+';END='+str(int(rec[1])+1)+';RE='+str(int(rec[3])+1)+'\tGT\t'+geno+'\n'
            myvcf += (contig +'\t'+ str(int(rec[1]))+'\t'+ '.'+'\t'+ '.'+'\t'+ '.'+'\t'+ str(int(rec[3])+1)+'\t'+ 'PASS'+'\t'+recinfo)


    with open(outputpath, "w") as f:
        f.write(myvcf)
def predict_fn(datapath = './data/', weightpath = 'detect', bamfilepath = '', tmppath = './prediction/', window_size = 200, max_worker = 16, vcfname = '', includecontig = []):
    model = init_model()
    model.load_weights(weightpath)
    window_size = 200
    maxread = 18
    genotype = False
    svlist = []
    predictiondict = {'DEL':0, 'INS':0, 'DEL_GENOTYPE':0, 'INS_GENOTYPE':0, 'DELINS': 0}
    count = 0

    if(datapath[-1] != '/'):
        datapath += '/'
    filelist = [datapath+f for f in listdir(datapath) if(isfile(join(datapath, f)) and 'npz' in f)]

    if(tmppath[-1] != '/'):
        tmppath += '/'
    if(os.path.isdir(tmppath) == True):
        for path in [tmppath+f for f in listdir(tmppath) if(isfile(join(tmppath, f)) and 'npy' in f)]:
            os.remove(path)
    else:
        os.mkdir(tmppath[:-1])
    filecount = 0
    oldpercentage = 0
    st = time.time()
    for filepath in filelist:
        contig = filepath.split(',')[0].split('/')[-1]
        if((includecontig != []) and (contig not in includecontig)):
            continue
        filecount += 1
        percentage = int((filecount/len(filelist))*100)
        if(percentage > oldpercentage):
            oldpercentage = percentage
            usedtime = (time.time() - st)
            st = time.time()
            print(str(oldpercentage)+'% Finished. ' +'ETA: ' + str(usedtime*(100-oldpercentage)//60)+' minutes.')
        Xy = np.load(filepath)

        contig = filepath.split(',')[0].split('/')[-1]
        index = Xy['index']
        data = Xy['data']
        if(data.size % (100* window_size*maxread) != 0):
            tmpdata, taildata = batchdata(data, 100, window_size, maxread)
            tmpdata = tf.data.Dataset.from_tensor_slices(tmpdata).batch(100)
            taildata = tf.data.Dataset.from_tensor_slices(taildata).batch(100)
            prediction = tf.concat([(model.predict(tmpdata) > 0.5).reshape(-1, 1), (model.predict(taildata) > 0.5).reshape(-1, 1)], axis = 0).numpy().astype('float32')
        else:
            tmpdata = tf.data.Dataset.from_tensor_slices(data.reshape(-1, 100, window_size, maxread, 1)).batch(100)
            prediction = (model.predict(tmpdata) > 0.5).reshape(-1, 1).astype('float32')
            
        tmpsvlist = []
        tmpsvlist = call_deletion(tmpsvlist, contig, index, data, prediction[:,0], genotype = genotype, genotypelist = predictiondict['DEL_GENOTYPE'], window_size = window_size)


        while(True):
            if(len(multiprocessing.active_children()) < max_worker):
                #print('Predicting on ', datapath)
                multiprocessing.Process(target=baseinfo_main_binsaver, args=(bamfilepath, tmpsvlist, count, tmppath)).start()
                count += 1
                break
            else:
                time.sleep(2)
    while(True):
        if(len(multiprocessing.active_children()) > 0):
            time.sleep(30)
        else:
            break
    count = 0
    l = []
    for path in os.listdir(tmppath):
        if('npy' in path):
            tmparray = np.load(tmppath + path)
            if(tmparray.size == 0):
                continue
            l.append(np.load(tmppath + path))
            count += 1
    bamfile = pysam.AlignmentFile(bamfilepath, 'rb')
    contig2length = {}
    for count in range(len(bamfile.get_index_statistics())):
        contig2length[bamfile.get_index_statistics()[count].contig] = bamfile.lengths[count]
    if(vcfname == ''):
        vcfname = 'Breaknet_'+bamfilepath.split('/')[-1][:-3] + 'vcf'
    tovcf(np.concatenate(l), contig2length, vcfname)
