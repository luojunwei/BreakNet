
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pysam

def allocation(cigartuples, posinref, readarray, qualityarray, start, end): 
  posinread = 0

  padinsertinfo = []

  softclipinfo = []

  cigararray = np.array(-1)

  itemcount = 0


  slicestart = 0
  reachstart = False
  reachend = False
  softclippresent = False
  softclipcount = 0
  totalcount = readarray.shape[0]

  if(posinref <= start): 


    itemcount = -1

    for item in cigartuples:

      itemcount = itemcount + 1

      if(item[0] in [0, 2, 7, 8]):

        if(posinref + item[1] > start): # Meet start. Slice required.

          #print('front early stopping')
          #print('posinref is', posinref)
          #cigararray = np.append(cigararray, np.ones(posinref + item[1] - start) * item[0])
          cigartuples[itemcount] = [item[0], int(posinref + item[1] - start)] 
          #print(item[1])
          #print('push in cigararray ', cigartuples[itemcount])

          if(item[0] in [0, 7, 8]): # For reference match
            
            slicestart = posinread + start - posinref
            posinread = posinread + start - posinref
            posinref = start

          else: # For deletion

            slicestart = posinread
            posinref = start

          reachstart = True

          break


        if(item[0] in [0, 7, 8]): 

          posinread = posinread + item[1]
          posinref = posinref + item[1]

        else: # For deletion

          posinref = posinref + item[1]

        continue

      
      if(item[0] in [1, 4]): 

        posinread = posinread + item[1]
        if(item[0] == 4):
          softclippresent = True
          softclipcount = softclipcount + item[1]
        continue
      
      if(item[0] == 5):
        
        softclippresent = True
        softclipcount = softclipcount + item[1]
        totalcount = totalcount + item[1]
        continue
        
      print(str(item[0]),' CIGAR TAG DONT SUPPORT ')
      return None
        

  else: 

    if(posinref >= end): 

      return None, np.ones(int(end - start)) * 6, np.column_stack((np.arange(start, end).reshape(int(end - start), 1), np.zeros((int(end - start), 1)))), np.ones(int(end - start)).reshape(int(end - start), 1) * 0
    
    else: 

      readfrontpadding = np.column_stack((np.arange(start, posinref).reshape(int(posinref - start), 1), np.zeros((int(posinref - start), 1))))
      readarray = np.row_stack((readfrontpadding, readarray))
      qualityarray = np.row_stack((np.zeros((int(posinref - start), 1)), qualityarray))
      cigararray = np.ones(int(posinref - start + 1)) * 6 # add 1 to deal with last slice
      posinread = posinread + posinref - start
      reachstart = True


  if(reachstart == False): 
    
    return None, np.ones(int(end - start)) * 6, np.column_stack((np.arange(start, end).reshape(int(end - start), 1), np.zeros((int(end - start), 1)))), np.ones(int(end - start)).reshape(int(end - start), 1) * 0
  
  #print('posinref is', posinref)
  #print('read position is ', readarray[slicestart: slicestart +10])

  
  if(posinref < end):
    
    usefront = False
    usetail = True
    
    locincigartuples = 0
    
    cliploc = []
    
    for item in cigartuples[itemcount:]:
      
      locincigartuples = locincigartuples + 1
      #print('current cigar is ', item, 'posinref is ', posinref, 'posinread is ', posinread, 'start is ', start, 'end is ', end)
      #For reference match: posinref is match part first base
      if(item[0] in [0, 7, 8]):
        
        usefront = True
        
        if(end <= item[1] + posinref): 

          #print('posinref, end ', posinref, end)
          cigararray = np.append(cigararray, np.ones(end - posinref) * item[0])
          readarray = readarray[: posinread + end - posinref]
          qualityarray = qualityarray[: posinread + end - posinref]
          posinref = end
          break
          
        else: 

          cigararray = np.append(cigararray, np.ones(item[1]) * item[0])
            
          posinref = posinref + item[1]

          posinread = posinread + item[1]
          usefront = True

          continue
        
      #For insertion
      if(item[0] == 1):

        cigararray = np.append(cigararray, np.ones(item[1]) * item[0])

        posinread = posinread + item[1]

        padinsertinfo.append([int(posinref - 1), int(item[1])])


        continue

   
      if(item[0] == 2): 
        
        usefront = True
        
        if(end <= item[1] + posinref): 
            
          #print('del posinref, end ', posinref, end)
          cigararray = np.append(cigararray, np.ones(end - posinref) * item[0])
          readarray = np.row_stack((readarray[: posinread], np.column_stack((np.arange(posinref, end).reshape(end - posinref, 1), np.zeros(end - posinref).reshape(end - posinref, 1)))))
          qualityarray = np.row_stack((qualityarray[: posinread], np.zeros(end - posinref).reshape(end - posinref, 1)))

          posinref = end 
          break
          
        else: #Not yet. Following parameter need update: posinref

          cigararray = np.append(cigararray, np.ones(item[1]) * item[0])
          paddelarray = np.column_stack((np.arange(posinref, posinref + item[1]).reshape(item[1], 1), np.zeros(item[1]).reshape(item[1], 1)))
          #print(paddelarray.shape)
          readarray = np.row_stack((np.row_stack((readarray[:posinread], paddelarray)), readarray[posinread:]))
          qualityarray = np.row_stack((qualityarray[: posinread], np.row_stack((np.zeros(item[1]).reshape(item[1], 1), qualityarray[posinread: ]))))

          posinref = posinref + item[1]
          posinread = posinread + item[1]

          usefront = True


          continue

      #For soft clip: posinref is first base follow with this softclip part
      if(item[0] == 4): # We have not meet the end. Following parameter need update: posinread.
        softcigartmp = np.ones(item[1]) * item[0]
        if(locincigartuples == len(cigartuples[itemcount:])):
         
          usetail = False
        
        else:
          
          usetail = False
          for cigarinfo in cigartuples[itemcount + locincigartuples:]:
            if(cigarinfo[0] in [0, 2, 7, 8]):
              usetail = True
              

        if(usefront == True and usetail == True):

          softclipinfo.append([[int(posinref - 1), int(posinref)], readarray[posinread: posinread + item[1], 1]])
          
          softcigartmp[0], softcigartmp[-1] = softcigartmp[0] + 0.1, softcigartmp[-1] + 0.1
          
          cliploc.append(int(posinref - 1))
          cliploc.append(int(posinref))
        
        elif(usefront == True):

          softclipinfo.append([[int(posinref - 1)], readarray[posinread: posinread + item[1], 1]])
          
          softcigartmp[0] = softcigartmp[0] + 0.1
          
          cliploc.append(int(posinref - 1))
        
        else:
          
          softclipinfo.append([[-int(posinref)], readarray[posinread: posinread + item[1], 1]])
          
          softcigartmp[-1] = softcigartmp[-1] + 0.1
          
          cliploc.append(int(posinref))
          
          


        
        cigararray = np.append(cigararray, softcigartmp)
          
        posinread = posinread + item[1]

        padinsertinfo.append([int(posinref - 1), int(item[1])])

        softclippresent = True

        softclipcount = softclipcount + item[1]


        continue
        
      if(item[0] == 5): # Pass the hardclip cigar string.
        
        totalcount = totalcount + item[1]
        
        softclippresent = True

        softclipcount = softclipcount + item[1]
        
        if(locincigartuples == len(cigartuples[itemcount:])):
         
          usetail = False
        
        else:
          
          usetail = False
          for cigarinfo in cigartuples[itemcount + locincigartuples:]:
            if(cigarinfo[0] in [0, 2, 7, 8]):
              usetail = True
              
        hardlen = 1
        tmppadnone = np.array([[None]])
        
        
        if(usefront and usetail):
          hardlen = 2
          tmppadnone = np.array([[None], [None]])
          
          cliploc.append(int(posinref - 1))
          cliploc.append(int(posinref))
        elif(usefront == True):
          
          cliploc.append(int(posinref - 1))
        
        else:
          
          cliploc.append(int(posinref))
          
          
        cigararray = np.append(cigararray, np.ones(hardlen) * item[0])

        padinsertinfo.append([int(posinref - 1), hardlen])
        
        
        paddelarray = np.column_stack((tmppadnone, np.zeros(hardlen).reshape(hardlen, 1)))
          #print(paddelarray.shape)
        if(posinread == 0):
          readarray = np.row_stack((paddelarray, readarray[posinread:]))
          qualityarray = np.row_stack((np.zeros(hardlen).reshape(hardlen, 1), qualityarray[posinread: ]))
        elif(posinread == readarray.shape[0]):
          readarray = np.row_stack((readarray, paddelarray))
          qualityarray = np.row_stack((qualityarray, np.zeros(hardlen).reshape(hardlen, 1)))
        else:
          readarray = np.row_stack((np.row_stack((readarray[:posinread], paddelarray)), readarray[posinread:]))
          qualityarray = np.row_stack((qualityarray[: posinread], np.row_stack((np.zeros(hardlen).reshape(hardlen, 1), qualityarray[posinread: ]))))

        
        posinread = posinread + hardlen

        continue
        
      print(str(item[0]),' CIGAR TAG DONT SUPPORT ')
      return None


  #readarray = readarray[slicestart: posinread + end - posinref]
  
  #print(cigarinfodict)
  if(posinref < end): #Need padding.
    readrearpadding = np.column_stack((np.arange(posinref, end).reshape(end - posinref, 1), np.zeros(end - posinref).reshape(end - posinref, 1)))
    readarray = np.row_stack((readarray, readrearpadding))
    qualityarray = np.row_stack((qualityarray, np.zeros(end - posinref).reshape(end - posinref, 1)))
    cigararray = np.append(cigararray, np.ones(end - posinref) * 6)
  
  return np.array(padinsertinfo), cigararray[1:], readarray[slicestart:], qualityarray[slicestart:], [softclippresent, softclipcount / totalcount , softclipinfo, cliploc]

def paddel(AlignedSegment, start, end):

  grp = np.array(AlignedSegment.get_reference_positions(True))#[AlignedSegment.query_alignment_start: AlignedSegment.query_alignment_end]
  grp = grp.reshape(grp.size, 1)
  read = pd.DataFrame(list(AlignedSegment.query_sequence)).replace({'A': 10, 'T': 15, 'G': 20, 'C': 25}).values.reshape(grp.size,1)
  readarray = np.column_stack((grp, read))

  pos = 0
  i  = 0
  delpos = 0

  '''print(readarray.T[0].tolist())
  print(readarray.T[1].tolist())'''
  #print(grp[-10:].T)

  cigartuples = AlignedSegment.cigartuples
  #print(grp[-cigartuples[-1][1] - 10:].T)
  if(cigartuples == None):
    print('The alignment does not give CIGAR string')
    return None

  while(cigartuples[i][0] not in [0, 7, 8]): # Get TRUE reference start location
    if(cigartuples[i][0] != 5):
      pos = cigartuples[i][1] + pos
      if(cigartuples[i][0] == 2):
        delpos = cigartuples[i][1] + delpos

    i = i + 1

  '''if(grp[0] == None):#test for insertion or softclip on begin
    tmpbase = np.array([readarray[0][0] - 1, 0]).reshape(1,2)
    readarray = np.row_stack((tmpbase, readarray))'''
  #print('origin readarray is ', readarray.T.tolist())
  qqarray = np.array(AlignedSegment.query_qualities)
  if(qqarray.size == 1):
    qqarray = np.zeros((readarray.shape[0], 1))
  else:
    qqarray = qqarray.reshape((readarray.shape[0], 1))

  iarray, cigararray, readarray, qualityarray, softclipinfo = allocation(cigartuples, int(grp[pos] - delpos), readarray, qqarray, start, end)
  if(softclipinfo[0] == True):
    softclipinfo.append(read)
  else:
    softclipinfo.append([])

  return readarray, iarray, cigararray, qualityarray, softclipinfo


def dropselectcigartag(pcigararray, selectedtag, refseq, filter = False, preadarray = None, cutvalue = 1, keepfrontsoftclip = True):
  for tag in selectedtag:
    try:
      outputset = set(np.nonzero((pcigararray == tag).sum(axis = 0) == 0)[0].tolist()) & outputset
    except:
      outputset = set(np.nonzero((pcigararray == tag).sum(axis = 0) == 0)[0].tolist())

  if(filter == True):
    outputset = set(np.nonzero((preadarray > 0).sum(axis = 0, keepdims = True) > (cutvalue))[1].tolist()) & outputset

  if(keepfrontsoftclip):

    outputset = set(np.nonzero((pcigararray == 4.1).sum(axis = 0) != 0)[0].tolist()) & outputset

  outputset = set(np.nonzero(refseq > -1)[0].tolist()) | outputset
  return np.sort(np.array(list(outputset)))

def myshow(preadarray, pcigararray, qualityarray, refseq, seqbase, excludelist = [1, 4, 4.1], includecagtag = [2], filter = False, cutvalue = 1, maxdot = 200, softread = [], showpic = False, minrow = 18, maxrow = 18):
  #print(preadarray.shape, pcigararray.shape, qualityarray.shape)
  pltsoft = False
  if(len(softread) > 0):
    psoftarray = (np.array(softread[:,1]).reshape(preadarray.shape[0], 1).astype('float32') * (preadarray > 0).astype('float32'))
    pltsoft = True
  if(len(excludelist) > 0):

    columnkeeped = dropselectcigartag(pcigararray, excludelist, refseq, filter, preadarray, cutvalue)

    preadarray = preadarray[:, columnkeeped]
    qualityarray = qualityarray[:, columnkeeped]
    if(pltsoft == True):
      psoftarray = psoftarray[:, columnkeeped]

    slicedpcigararray = pcigararray[:, columnkeeped]
    for cigartag in includecagtag:
      try:
        spotarray = spotarray + (slicedpcigararray == cigartag) * cigartag
      except:
        spotarray = (slicedpcigararray == cigartag) * cigartag
  else:

    slicedpcigararray = pcigararray
    for cigartag in includecagtag:
      try:
        spotarray = spotarray + (slicedpcigararray == cigartag) * cigartag
      except:
        spotarray = (slicedpcigararray == cigartag) * cigartag

  cspotarray = spotarray[:,:200][np.argsort(spotarray[:200].sum(axis = 1))[::-1]]
  for loc in range(200, spotarray.shape[1], 200):
    tmp = spotarray[:,loc:loc+200][np.argsort(spotarray[:,loc:loc+200].sum(axis = 1))[::-1]]
    cspotarray = np.column_stack((cspotarray, tmp))
  spotarray=cspotarray 
  if(preadarray.shape[0] < minrow):
    pt = [preadarray, spotarray, psoftarray, qualityarray]
    padrownumber = minrow - preadarray.shape[0]
    column_number = preadarray.shape[1]
    for locinpt in range(4):
      
      pt[locinpt] = np.row_stack((pt[locinpt], np.zeros((padrownumber, column_number))))

    preadarray, spotarray, psoftarray, qualityarray = pt[0], pt[1], pt[2], pt[3]
  
  if(preadarray.shape[0] > maxrow and showpic == False):

    preadarray, spotarray, psoftarray, qualityarray = preadarray[:maxrow], spotarray[:maxrow], psoftarray[:maxrow], qualityarray[:maxrow]

    
  cliparray = np.ones((preadarray.shape[0], 1)) * seqbase

  if(showpic):
    loc = 0
    print('Read sequence')
    plt.matshow(preadarray)
    plt.show()
  else:
    fm = (spotarray>0).astype('float32')
    return fm.T.reshape(1, fm.size//(200 * 18), 200, 18, 1)





def pileupf(bamfile, contig, start, end, droplq = False, dropvalue = 0.8):
  end = start + ((end-start)//200+1)*200
  totalstarttime = time.time()
  locationlist = []
  readlist = []
  insertlist = []
  insertinfo = dict()
  keylist = []
  cigarlist = []
  qualitylist = []
  softcliplist = []

  debug = []
  
  depth = 0
  paddeltime = 0
  fetchtime = time.time()
  seqbase = np.zeros((1, end - start)).astype('float32')
  overlap = False
  
  for AlignedSegment in bamfile.fetch(contig, start, end):
    
    #debug.append(AlignedSegment)
    if(AlignedSegment.reference_start <= start):
        frontloc = 0
        whichstart = start
    else:
      frontloc = AlignedSegment.reference_start - start
      whichstart = AlignedSegment.reference_start
      
    if(AlignedSegment.reference_end >= (end - 1)):
      tailloc = end - start - 1
    else:
      tailloc = AlignedSegment.reference_end - start
    paddelstarttime = time.time()
    
    read, iarray, cigararray, qualityarray, softclipinfo = paddel(AlignedSegment, start, end)
    ratio = softclipinfo[1]
    tmp = seqbase[:,frontloc: tailloc + 1] + ratio
    for cclip in set(softclipinfo[3]):
      tmp[:,cclip - whichstart] = tmp[:,cclip - whichstart] + ratio
    fillblank = tmp.mean()
    seqbase[:,frontloc: tailloc + 1] = (np.column_stack((np.array([[fillblank]]), tmp))[:,:-1] + tmp + np.column_stack((tmp, np.array([[fillblank]])))[:,1:])/3
    
    
    if(droplq and softclipinfo[1] > dropvalue):
      continue
    overlap = True
    paddeltime = - paddelstarttime + time.time() + paddeltime
    locationlist.append(read[:,0])
    readlist.append(np.array(read[:,1]).astype('float32'))
    insertlist.append(iarray)
    cigarlist.append(cigararray.astype('float32'))
    qualitylist.append(qualityarray.flatten().astype('float32'))
    softcliplist.append(softclipinfo)
    depth = depth + 1
    '''print(read.shape)
    print(cigararray.shape)'''

    if(type(iarray) != np.ndarray):
      continue

    for item in iarray:
      cloc = item[0]
      if((item[0]) in insertinfo):
        if(cloc == lastloc):
          insertinfo[item[0]] = insertinfo[item[0]] + item[1]
        else:
          insertinfo[item[0]] = max(insertinfo[item[0]], int(item[1]))
      else:
        insertinfo[item[0]] = item[1]
        keylist.append(item[0])
      lastloc = cloc
  #print('fetch time = ', time.time() - fetchtime, time.time() - totalstarttime)
  
  keylist = np.sort(np.array(keylist))

  #print(keylist)
  '''print(readlist)
  print()
  print(insertinfo)'''

  #print(insertinfo)
  #readposlist = [0 for i in range(len(readlist))]

  refseq = np.arange(start, end)
  bias = 0
  for key in keylist:
    if(key == (end - 1)):
      print('end in insertioninfo')
      return 0
      
    insert = - np.ones(insertinfo[key])
    refseq = np.append(np.append(refseq[:key+1 + bias-start], insert), refseq[key+1+bias-start: ])
    bias = bias + insertinfo[key]
  
  if(overlap == False):
    return np.zeros((1, end - start)), np.ones((1, end - start)) * 6, np.zeros((1, end - start)), refseq, np.array([[False, 0]]), seqbase
  readcount = 0
  state = False
  readlisttime = time.time()
  
  cc = 0
  
  if(True):
    pallarray = 'None'
    for read in readlist:
    
      refloc = start
      locinread = 0
      insertcount = 0
      parray = 'None'


      for key in keylist:

        #print(insertlist[readcount][insertcount])
        while(True):
          slicelength = key + 1 - refloc
          refloc = refloc + slicelength
          insertsizeofrfortkey = 0
          insertpresentonkey = False
          if(insertlist[readcount].shape[0] > insertcount and key == insertlist[readcount][insertcount][0]):
            
            onreadinsert = insertlist[readcount][insertcount][1]

            insertpresentonkey = True
            while(insertlist[readcount].shape[0] > (insertcount + 1) and key == insertlist[readcount][insertcount + 1][0]):
              insertcount = insertcount + 1
              onreadinsert = onreadinsert + insertlist[readcount][insertcount][1]
            slicelength = slicelength + onreadinsert


          tmpparray = np.array([read[locinread: locinread + slicelength], cigarlist[readcount][locinread: locinread + slicelength], qualitylist[readcount][locinread: locinread + slicelength]])
          locinread = locinread + slicelength

          if(insertpresentonkey):

            remainlengh = insertinfo[key] - onreadinsert
            insertcount = insertcount + 1
            if(remainlengh > 0):
              
              tmpparray = np.column_stack((tmpparray, np.ones((3, remainlengh)) * np.array([[0.], [6.], [0.]])))
          
          else:

            tmpparray = np.column_stack((tmpparray, np.ones((3, insertinfo[key])) * np.array([[0.], [6.], [0.]])))


          try:

            parray = np.column_stack((parray, tmpparray))
          
          except:

            parray = tmpparray
          
  
          break

      
      if(locinread !=  len(read) and len(keylist) > 0):

        parray = np.column_stack((parray, np.array([read[locinread: ], cigarlist[readcount][locinread: ], qualitylist[readcount][locinread: ]])))
        
        #print()
      else:

        parray = np.array([read[locinread: ], cigarlist[readcount][locinread: ], qualitylist[readcount][locinread: ]])
      
      readcount = readcount + 1
        
      if(state):

        pallarray = np.column_stack((pallarray, parray))
        #print(parray.shape)


          
      else:

        pallarray = parray
        readlength = parray.shape[1]
        state = True
        
  
  #print('readlisttime',time.time() - readlisttime)
  #print(time.time() - totalstarttime, paddeltime)
  
  return pallarray[0].reshape(int(pallarray.shape[1] / readlength), readlength), pallarray[1].reshape(int(pallarray.shape[1] / readlength), readlength), pallarray[2].reshape(int(pallarray.shape[1] / readlength), readlength), refseq, np.array(softcliplist), seqbase

def feature_matrix(bamfile, contig, start, end, showpic = False):

  preadarray, pcigararray, qualityarray, refseq, softcliplist, seqbase = pileupf(bamfile, contig, start, end)
  return myshow(preadarray, pcigararray, qualityarray, refseq, seqbase, softread =  softcliplist, showpic = showpic)

