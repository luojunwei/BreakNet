import sys
from model import train_fn, predict_fn
from feature import feature_matrix
import tensorflow as tf
tf.config.set_visible_devices([], 'GPU')
mode = sys.argv[1]
if(mode == 'data_mode'):
    print('Produce data')
    try:
        bamfilepath, contig, start, end, outputpath, vcfpath = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]
    except:
        bamfilepath, contig, start, end, outputpath, vcfpath = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], ''
    print('bamfile path ', bamfilepath)
    print('output path ', outputpath)
    print('vcf path ', vcfpath)
    feature_matrix(bamfilepath, str(contig), int(start), int(end), outputpath, vcfpath)
    print('Completed')
elif(mode == 'train_mode'):
    print('training')
    print('Best weights as saved as "trainedweights"')
    traindatapath,  testdatapath, epochs = sys.argv[2], sys.argv[3], sys.argv[4]
    print('training data path ', traindatapath)
    print('evaluation data path ', testdatapath)
    train_fn(traindatapath,  testdatapath, int(epochs))
    print('Completed')
elif(mode == 'call_mode'):
    print('testing')
    datapath, weightpath = sys.argv[2], sys.argv[3]
    print('testing bamfile path ', datapath)
    print('weight path ', weightpath)
    predict_fn(datapath, weightpath)
    print('Completed, Result saved at current folder as result.csv')
else:
    print('Useage')
    print('To produce data for training')
    print('python breaknet.py data_mode bamfilepath contig start end outputpath vcfpath')
    print('To produce data for call sv')
    print('python breaknet.py data_mode bamfilepath contig start end outputpath')
    print('Train a new model')
    print('python breaknet.py train_mode traindatapath  evaluationdatapath epochs')
    print('To call sv')
    print('python breaknet.py call_mode datapath trainedweightspath')
    print('BAM file should be sorted and indexed\ncontig is name of contig in bam file\nstart/end are positions in reference\noutputpath store information of reads\nvcf file use to label data for training\ntraindata,  evaluationdata used to train and evaluate model')
    
