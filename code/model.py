import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.layers import Input, Dense, Conv2D, DepthwiseConv2D, BatchNormalization, Dropout, GlobalAveragePooling2D, Reshape, multiply, add, Activation
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

    a = tf.reshape(tf.cast(~((cigararray[:,0] == 1) | (cigararray[:,0] == 4)), 'float32'), [cigararray.shape[0], 1]) * cigararray
    a1 = tf.reshape(a[:,1], [reccount, cigarweight])

    a = tf.concat([cigararray, tf.reshape(tf.matmul(a1, tf.linalg.band_part(tf.ones([a1.shape[1] , a1.shape[1]], tf.float32), 0, -1)) + refstartarray, [cigararray.shape[0], 1])], axis = 1)


    return tf.boolean_mask(a, (start <= a[:,-1]) & (a[:,-1] < end))
    return tf.boolean_mask(a, (start <= a[:,-2]) & (a[:,-2] < end) & (a[:,0] == 1))[:,1:]



def baseinfo(bamfile, contig, start, end):
    

    cigararray = []
    readstartandend = []
    refpositonlist = []
    refpositonweight = []
    substitionarray, deletionarray, substitionweight, deletionweight = [], [], [], []
    nooverlap = True
    qualityarray = []

    for AlignedSegment in bamfile.fetch(contig, start, end):
    



        cigararray.append(tf.keras.backend.flatten(tf.constant(AlignedSegment.cigartuples)))
        readstartandend.append([AlignedSegment.reference_start-start, AlignedSegment.reference_end-start, AlignedSegment.mapping_quality, (1 - (AlignedSegment.query_alignment_length / AlignedSegment.infer_read_length()))**2])

        nooverlap = False

    
    if(nooverlap):
        print(dsahjdaj)
    readstartandend = tf.constant(readstartandend, tf.float32)
    cigararray = tf.keras.preprocessing.sequence.pad_sequences(cigararray)
    reccount, cigarweight = cigararray.shape[0], int(cigararray.shape[1] / 2)
    cigararray = cigararray.reshape(int(cigararray.size / 2), 2)

    cigararray = cigarinfo(cigararray, readstartandend[:,:1], 0, end - start, reccount, cigarweight).numpy().astype('int64')
    a = cigararray[(cigararray[:,0] == 2) & (cigararray[:,1] > 20)]
    if(a.size == 0):
        return []
    a[:,-1] = a[:,-1] - a[:,-2] 
    delsig = np.column_stack((a[:,-1:], a[:,-2:-1]))



    loc = np.array(list(delsig))[:,0]
    binnum = 20
    binresult = (loc//binnum)
    mc = Counter(binresult).most_common(1)[0]
    sp = mc[1]
    minv, maxv = mc[0]-1, mc[0]+1
    tmp = np.median(np.array(list(delsig))[(minv<=binresult) *  (maxv>= binresult)], axis = 0).astype('int64')
    tmp[0] = tmp[0]+start
    return [contig] + tmp.tolist()+[sp, 'DEL']
    
def baseinfo_main_binsaver(bamfilepath, delloc):




    bamfile = pysam.AlignmentFile(bamfilepath, 'rb', threads = 20)


    delsig = []
    for rec in delloc:
        contig, start, end = str(rec[1]), int(rec[2]), int(rec[3])
        bk = baseinfo(bamfile, contig, start, end)
        if(len(bk) == 0):
            continue
        delsig.append(bk)
    return delsig



def mergedeletioncall(predict, contig, index, bamfilepath, window_size = 200):
    #contig: contig name in bam file, predict: 1d array, startloc: first prediction location in reference, window_size: the sub-region of a feature matrix included
    predict = (predict > 0.5).flatten()
    svlist = []
    count = -1
    insv = False
    for p in predict:
        count += 1
        if(p > 0):
            if(insv == False):
                svstartloc = index[count]
                insv = True

        else:
            if(insv == True):
                svlist.append(['DEL', contig, svstartloc, index[count]])
                insv = False

                continue
    if(bamfilepath != ''):
        return baseinfo_main_binsaver(bamfilepath, svlist)
    else:
        return [[rec[1], rec[2], rec[2] - rec[3], 'NA', 'DEL']]
def train_fn(traindatapath = './',  testdatapath = './', trainedweightspath = 'savedweight', epochs = 1000):
    model = init_model()
    trainfilelist = [traindatapath+f for f in listdir(traindatapath) if(isfile(join(traindatapath, f)) and 'npz' in f)]
    count = 0
    while(count<len(trainfilelist)):
        Xy = np.load(trainfilelist[count])
        count += 1
        if(Xy['data'].shape[1] != 100):
            continue
    train = tf.data.Dataset.from_tensor_slices((Xy['data'].astype('float32'), Xy['label'].astype('float32').reshape(-1, 100))).batch(80)
    for loc in range(count, len(trainfilelist)):
        
        Xy = np.load(trainfilelist[loc])
        if(Xy['data'].shape[1] != 100):
            continue
        tmpds = tf.data.Dataset.from_tensor_slices((Xy['data'].astype('float32'), Xy['label'].astype('float32').reshape(-1, 100))).batch(80)
        train = train.concatenate(tmpds)

    testfilelist = [testdatapath+f for f in listdir(testdatapath) if(isfile(join(testdatapath, f)) and 'npz' in f)]
    count = 0
    while(count<len(testfilelist)):
        Xy = np.load(testfilelist[count])
        count += 1
        if(Xy['data'].shape[1] != 100):
            continue
    test = tf.data.Dataset.from_tensor_slices((Xy['data'].astype('float32'), Xy['label'].astype('float32').reshape(-1, 100))).batch(80)
    for loc in range(count, len(testfilelist)):
        Xy = np.load(testfilelist[loc])
        if(Xy['data'].shape[1] != 100):
            continue
        tmpds = tf.data.Dataset.from_tensor_slices((Xy['data'].astype('float32'), Xy['label'].astype('float32').reshape(-1, 100))).batch(80)
        test = train.concatenate(tmpds)    
    bestf1 = 0
    for epoch in range(epochs):
        model.fit(train)
        h = model.evaluate(test)
        r, p = h[-1], h[-2]
        f1 = 2*(r*p)/(r+p+1e-10)
        if(bestf1<f1):
            model.save_weights(trainedweightspath)
            bestf1 = f1
            print()
            print('New best F1: ', bestf1)
            print()
    return model
def predict_fn(datapath = './', weightpath = 'trainedweights', bamfilepath = '', window_size = 200):
    model = init_model()
    model.load_weights(weightpath) 
    filelist = [datapath+f for f in listdir(datapath) if(isfile(join(datapath, f)) and 'npz' in f)]
    resultlist = [['CONTIG', 'START', 'SVLEN', 'READ_SUPPORT', 'SVTYPE']]
    for datapath in filelist:
        Xy = np.load(datapath)
        index = Xy['index']
        data = tf.data.Dataset.from_tensor_slices((Xy['data'])).batch(80)
        predict = model.predict(data)>0.5
        datapath = datapath.split('/')[-1].split(',')
        contig, start = datapath[0], datapath[1]
        resultlist += mergedeletioncall(predict, contig, index, bamfilepath, window_size)
    df = pd.DataFrame(resultlist)
    new_header = df.iloc[0] 
    df = df[1:] 
    df.columns = new_header
    (df.sort_values(['CONTIG', 'START']).to_csv('result.csv'))
