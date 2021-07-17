import sys
from model import train_fn, predict_fn
from feature import feature_matrix
import tensorflow as tf
tf.config.set_visible_devices([], 'GPU')
mode = sys.argv[1]
debug = 0
if(mode == 'data_mode'):
    if(len(sys.argv) != 8 and len(sys.argv) != 7):
        debug = 1
    else:
        print('Produce data')
        try:
            bamfilepath, contig, start, end, outputpath, vcfpath = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7]
        except:
            bamfilepath, contig, start, end, outputpath, vcfpath = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], ''
        print('bamfile path ', bamfilepath)
        print('output path ', outputpath)
        print('vcf path ', vcfpath)
        feature_matrix(bamfilepath, str(contig), int(start), int(end), outputpath, vcfpath)
        print('\n\n')
        print('Completed')
        print('\n\n')
elif(mode == 'train_mode'):
    if(len(sys.argv) != 6):
        debug = 1
    else:
        print('training')
        
        traindatapath,  testdatapath, trainedweightspath, epochs = sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]
        print('Best weights as saved as '+trainedweightspath)
        print('training data path ', traindatapath)
        print('evaluation data path ', testdatapath)
        train_fn(traindatapath,  testdatapath, trainedweightspath, int(epochs))
        print('\n\n')
        print('Completed')
        print('\n\n')
elif(mode == 'call_mode'):
    if(len(sys.argv) != 4):
        debug = 1
    else:
        print('testing')
        datapath, weightpath = sys.argv[2], sys.argv[3]
        print('testing bamfile path ', datapath)
        print('weight path ', weightpath)
        predict_fn(datapath, weightpath)
        print('\n\n')
        print('Completed, Result saved at current folder as result.csv')
        print('\n\n')
else:
    debug = 1
if(debug ==1):
    print('\n\n')
    print('Useage')
    print('To produce data for training')
    print('python breaknet.py data_mode bamfile_path chromosome_name start_position end_position output_data_folder vcf_path')
    print('To produce data for call sv')
    print('python breaknet.py data_mode bamfile_path chromosome_name start_position end_position output_data_folder')
    print('Train a new model')
    print('python breaknet.py train_mode training_data_folder evaluation_data_folder trained_weight_path epochs')
    print('To call sv')
    print('python breaknet.py data_mode bamfile_path chromosome_name start_position end_position call_folder')
                 
