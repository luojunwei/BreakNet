# BreakNet

A deep learning method to detect deletion from long reads alignment. It is
built with **Tensorflow** and **Python 3**.


## Installation

### Requirements
  * python 3.7, numpy, scipy, pandas, Matplotlib, TensorFlow 2.4, pysam

#### 1. Create a virtual environment

```bash
# create
conda create -n BreakNet python=3.6
# activate
conda activate BreakNet
# deactivate
conda deactivate
```

#### 2. clone BreakNet
- After creating and activating the BreakNet virtual environment, download breakNet from github:
```bash
git clone https://github.com/luojunwei/BreakNet.git
cd BreakNet
``` 
#### 3. Install

```bash
conda activate Breaknet
conda install numpy, scipy, pandas, Matplotlib, TensorFlow, pysam

``` 

## Tested data
The example data can be downloaded from 
### HG002
https://ftp.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/HG002_PB_70x_RG_HP10XtrioRTG.bam
### HG00514
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160905_smithm_pacbio_aligns/HG00514_bwamem_GRCh38DH_CHS_20160905_pacbio.bam
### HG00733
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160905_smithm_pacbio_aligns/HG00733_bwamem_GRCh38DH_PUR_20160905_pacbio.bam
### NA19240
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160905_smithm_pacbio_aligns/NA19240_bwamem_GRCh38DH_YRI_20160905_pacbio.bam


## Usage



#### 1. extract features
```bash
from feature import feature_matrix
bamfile = pysam.AlignmentFile("*/HG002_PB_70x_RG_HP10XtrioRTG.bam", "rb")
contig, start, end = '1', 2431832-1300, 2431832+1000 
fm = feature_matrix(bamfile, contig, start, end)
```

#### 2. load model and make prediction
```bash
from model import init_model
model = init_model()
model.load_weights('timedist_w2down10_17andfullna19240hg002')
prediction = model.predict(fm)>0.5
```
