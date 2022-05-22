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
                
        
                            
                    
                
                

       
                            
                    
                
                



                
        
                            
                    
                
                


            
    
