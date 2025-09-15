import os
import json

from MOFA_functions import *

## LEARNING DATASET FACTORIZATION  WITH MOFA

## INIT PARAMETERS
Prjct = 'Lrn_4test'
TopD = '5000D'
Prior_K = 50
omics = ['mRNA','miRNA', 'DNAme', 'SNV']
likelihoods = ['gaussian', 'gaussian', 'gaussian', 'bernoulli']
M = len(omics)
G = 1 ## groups - always set to 1
data_options = {
    "scale_views": False,
    "scale_groups": False,
    "center_groups": True
}
training_options = {
    "TrainingIter": 1000, ## max training iterations
    "freqELBOChk": 5, ## how often to check the elbo, this is the default for R so using this
    "drop_factor_threshold": 0.01,
    "seed": 1234567
}

## DIRECTORIES
BasePath = os.getcwd()
InputDir = os.path.join(BasePath, Prjct + '_' + TopD)
OutputDir = os.path.join(InputDir,'Fctrzn_' + str(Prior_K) + 'K_' + str(training_options["drop_factor_threshold"])[2:]+'TH')
if not os.path.exists(OutputDir):
    os.mkdir(OutputDir)
    
## IMPORT METADATA
Lrn_meta_in = os.path.join(InputDir,"Lrn_meta.json")
with open(Lrn_meta_in) as infile:
    Lrn_meta = json.load(infile)

## IMPORT OMICS DATA
data = import_data(G = G, M = M, omics = omics, InputDir = InputDir, metadata = Lrn_meta)

## RUN MOFA
script_start_time = time()
ent = runMOFA(data, data_options, likelihoods, omics, Prior_K, training_options)
script_end_time = time()

## SAVE MODEL AND RESULTS
saveResultsModel(ent, OutputDir, training_options, data_options, G, likelihoods, script_start_time, script_end_time)

