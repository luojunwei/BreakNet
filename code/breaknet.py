import sys
from model import train_fn, predict_fn, init_model
from feature import feature_matrix
import tensorflow as tf
tf.config.set_visible_devices([], 'GPU')

mode = sys.argv[1]
debug = 0
model = init_model()
if(mode == 'data_mode'):
    if(len(sys.argv) not in [4, 5, 6]):
        debug = 1
    else:
        print('Produce data')
        if(len(sys.argv) == 6):
            bamfilepath, outputpath, max_worker, includecontig = sys.argv[2], sys.argv[3], int(sys.argv[4]), [str(contig) for contig in eval(sys.argv[5])]
        if(len(sys.argv) == 5):
            bamfilepath, outputpath, max_worker, includecontig = sys.argv[2], sys.argv[3], int(sys.argv[4]), []
        if(len(sys.argv) == 4):
            bamfilepath, outputpath, max_worker, includecontig = sys.argv[2], sys.argv[3], 8, []
        print('bamfile path ', bamfilepath)
        print('output path ', outputpath)
        print('max_worker set to ', str(max_worker))
        if(includecontig == []):
            print('All chromosomes within bamfile will be used')
        else:
            print('Following chromosomes will be used')
            print(includecontig)
        feature_matrix(bamfilepath=bamfilepath, outputpath=outputpath, max_worker=max_worker, includecontig=includecontig)
        print('\n\n')
        print('Completed')
        print('\n\n')

elif(mode == 'call_mode'):
    if(len(sys.argv) not in  [5, 6]):
        debug = 1
    else:
        print('testing')
        if(len(sys.argv) == 6):
            datapath, weightpath, bamfilepath, includecontig = sys.argv[2], sys.argv[3], sys.argv[4], [str(contig) for contig in eval(sys.argv[5])]
        else:
            datapath, weightpath, bamfilepath, includecontig = sys.argv[2], sys.argv[3], sys.argv[4], []
        
        print('bamfile path ', bamfilepath)
        print('weight path ', weightpath)
        print('data file path ', datapath)
        if(includecontig == []):
            print('All chromosomes within bamfile will be used')
        else:
            print('Following chromosomes will be used')
            print(includecontig)
            


        predict_fn(datapath = datapath, weightpath = weightpath, bamfilepath = bamfilepath, includecontig=includecontig )
        print('\n\n')
        print('Completed, Result saved in current folder')
        print('\n\n')
else:
    debug = 1
if(debug ==1):
    print('\n\n')
    print('Useage')
    print('Produce data for call sv')
    print('python breaknet.py data_mode bamfile_path output_data_folder max_worker(default:8), includecontig(default:[](all chromosomes))')
    print('Call sv')
    print('python breaknet.py call_mode datapath, weightpath, bamfilepath, includecontig(default:[](all chromosomes)')
                 
