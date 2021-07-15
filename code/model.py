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
    model.compile(loss=tf.keras.losses.binary_crossentropy, optimizer=tf.optimizers.Adam())

    return tf.keras.Model(inputs, outputs)

def mergedeletioncall(predict, contig, startloc, window_size = 200):
    #contig: contig name in bam file, predict: 1d array, startloc: first prediction location in reference, window_size: the sub-region of a feature matrix included
    predict = (predict > 0.5).flatten()
    svlist = [['contig', 'start', 'end', 'SVTYPE']]
    loc = startloc-window_size
    insv = False
    for p in predict:
        loc += window_size
        if(p > 0):
            if(insv == False):
                svstartloc = loc 
                insv = True

        else:
            if(insv == True):
                svlist.append([contig, svstartloc, loc, 'DEL'])
                insv = False

                continue
    return pd.DataFrame(svlist)

