
# find ranges to run simulated annealing
# then run SetupBatchRuns.py

import json
import numpy as np
import math
import yaml
import os

# get from config: output_dir, registration_start_frame
config_name = 'CurrentConfig'
config_path = '../';
print('Config ',config_name)
with open(os.path.join(config_path,config_name,'config.yaml'), 'r') as file:
    config_opts = yaml.safe_load(file)

output_path = config_opts["output_dir"]
# read Features.json from output_dir
feature_path = os.path.join(output_path,'Features.json')
fid = open(feature_path, "r")
data = json.load(fid)
# find number of nuclei per frame
start_frame = 0
register_start_frame = 0
nframes = len(data['centroids'])
min_frame = 1000
max_frame = 0


# if bad segmentation or polar bodies - graph has exclusions
number_of_excludes_per_frame = 0
exclude_labels = {}  # key frame, value is label to exclude (just one here)
if ('exclude_start_frame' in config_opts.keys()):
    exclude_start_frame = config_opts['exclude_start_frame']
else:
    exclude_start_frame = -1
if (exclude_start_frame > 0):
    number_of_excludes_per_frame = 1
else:
    number_of_excludes_per_frame = 0
print('number of excluded labels per frame ',number_of_excludes_per_frame)



for iframe in range(nframes):
    s = np.asarray(data['mean_intensity'][iframe])
    ind = np.nonzero(s)
    nNuc = int((ind[0].shape[0])) - 1 - number_of_excludes_per_frame # don't include label 0
    if nNuc > 0:
        min_frame = iframe
        break

for iframe in range(nframes):
    s = np.asarray(data['mean_intensity'][iframe])
    ind = np.nonzero(s)
    nNuc = int((ind[0].shape[0])) - 1 - number_of_excludes_per_frame  # don't include label 0
    if nNuc > 0:
        max_frame = iframe

print('min and max frames of sequence ',min_frame,max_frame)



nframes = max_frame - min_frame + 1
# check that number of nuclei always increasing
prevnNuc = 0
badFrames = []
for iframe in range(nframes):
    s = np.asarray(data['sphericity'][iframe])
    ind = np.nonzero(s)
    nNuc = int((ind[0].shape[0])) - 1 - number_of_excludes_per_frame # don't include label 0
    #print(nNuc)
    if (nNuc < prevnNuc): # check that number of nuclei always increasing
        print('WARNING: number of nuclei decreased from ',prevnNuc, 'to ',nNuc,' at frame ',iframe)
        print('will not run annealing past this frame - there may be issues in making tasks')
        badFrames.append(iframe-1)
        prevnNuc = nNuc
        max_frame = iframe - 1
        badnNuc = nNuc
        break
    if (nNuc > 0) and (nNuc > prevnNuc):
        max_frame = iframe
        prevnNuc = nNuc

nframes = max_frame - min_frame + 1
print('min and max frames of sequence ',min_frame,max_frame)

# get last frame where there are 2^n nuclei
# get first frame where there are 2^n nuclei

last_frame ={} # key is 2^n
first_frame = {} # key is 2^n
for n in range(2,7):
    test = pow(2,n)
    for iframe in range(nframes):
        s = np.asarray(data['sphericity'][iframe])
        ind = np.nonzero(s)
        nNuc = int((ind[0].shape[0])) - 1 - number_of_excludes_per_frame# don't include label 0
        # last frame where there are 2^n
        if (nNuc == test):
            last_frame[test] = iframe
    for iframe in range(nframes):
        s = np.asarray(data['sphericity'][iframe])
        ind = np.nonzero(s)
        nNuc = int((ind[0].shape[0])) - 1 - number_of_excludes_per_frame # don't include label 0
        if (nNuc == test):
            first_frame[test] = iframe
            break
            


print('first frames ',first_frame)
print('last frames ',last_frame)

# make dictionary with start,nframes of each 1-to-1 segment
NoSplitSeg = {} # key is number of nuclei
for n in range(2,7):
    test = pow(2,n)
    if test in first_frame.keys():
        start_frame = first_frame[test]
        end_frame = last_frame[test]
        nframes = end_frame - start_frame + 1
        if nframes > 1:
            NoSplitSeg[test] = [start_frame,end_frame, nframes]
print('NoSplitSeg ',NoSplitSeg)

# make dictionary with start/end for each doubling segment
SplitSeg = {} # key is number of nuclei at start
for n in range(2,7):
    test = pow(2,n)
    if test in last_frame.keys():
        start_frame = last_frame[test]
        end_test = test*2
        if end_test in first_frame.keys():
            end_frame = first_frame[test*2]
            nframes = end_frame - start_frame + 1
            SplitSeg[test] = [start_frame,end_frame, nframes]
print('SplitSeg ',SplitSeg)


# make dictionary with start/end for each partial-doubling segment
PartialSplitSeg = {} # key is number of nuclei at start of segment
# could be at either beginning and/or end of sequence
# first test if at beginning
bBeginClean = False
s = np.asarray(data['volumes'][min_frame])
ind = np.nonzero(s)
nNuc = int(ind[0].shape[0]) - 1 - number_of_excludes_per_frame
if (math.log2(nNuc) == int(math.log2(nNuc))):
        bBeginClean = True
if (not bBeginClean):
    start_frame = min_frame
    # doubling ends when nNuc at start_frame reaches a multiple of two
    s = np.asarray(data['volumes'][start_frame])
    ind = np.nonzero(s)
    nNuc = int(ind[0].shape[0]) - 1 - number_of_excludes_per_frame
    nPow = math.floor(math.log2(nNuc)) 
    nMult = 2**nPow
    #print('partial doubling from ',start_frame,last_frame[nMult])
    if (nMult in last_frame):
        end_frame = last_frame[nMult]
        nframes = end_frame - start_frame + 1
        if (nframes > 1):
            PartialSplitSeg[nNuc] = [start_frame, end_frame, nframes]
            #print('Begin PartialSplitSeg ',PartialSplitSeg)
            
# now looking for partial doubling at end
bEndClean = False
s = np.asarray(data['volumes'][max_frame])
ind = np.nonzero(s)
nNuc = int(ind[0].shape[0]) - 1 - number_of_excludes_per_frame
if (math.log2(nNuc) == int(math.log2(nNuc))):
    bEndClean = True

# is last seg a split or nosplit
largest_frame_split = 0
for split in SplitSeg.values():
    if (split[1] > largest_frame_split):
        largest_frame_split = split[1]
largest_frame_nosplit = 0
for nosplit in NoSplitSeg.values():
    if (nosplit[1] > largest_frame_split):
        largest_frame_nosplit = nosplit[1]
if (largest_frame_split > largest_frame_nosplit) and (not bEndClean):
    partial_end_frame = max_frame
    # doubling started when nNuc at partial_end_frame was a multiple of two
    s = np.asarray(data['volumes'][partial_end_frame])
    ind = np.nonzero(s)
    nNuc = int(ind[0].shape[0]) - 1 - number_of_excludes_per_frame 
    nPow = math.floor(math.log2(nNuc))
    nMult = 2**nPow
    start_frame = last_frame[nMult]
    end_frame = partial_end_frame
    nframes = end_frame - last_frame[nMult]  + 1
    PartialSplitSeg[nMult] = [start_frame, end_frame, nframes]

print('PartialSplitSeg ',PartialSplitSeg)

# make taskfile.txt for running in disBatch
fid = open('./taskfile.txt','w')
ijob = 1
for k in NoSplitSeg.values():
    job_name = 'job' + str(ijob)
    cmd_str = './RunUnbalanced.sh ' + job_name + ' -c ' + config_name + ' -s ' + str(k[0]) + ' -e ' + str(k[1]) + ' > Results/' + config_name + '_' + str(k[0]) + '_' + str(k[1]) +  '.log' + '\n'
    fid.write(cmd_str)
    ijob = ijob + 1

for k in SplitSeg.values():
    job_name = 'job' + str(ijob)
    cmd_str = './RunBalanced1.sh ' + job_name + ' -c ' + config_name + ' -s ' + str(k[0]) + ' -e ' + str(k[1]) + ' > Results/' + config_name + '_' + str(k[0]) + '_' + str(k[1]) +  '.log' + '\n'
    fid.write(cmd_str)
    ijob = ijob + 1

for k in PartialSplitSeg.values():
    job_name = 'job' + str(ijob)
    cmd_str = './RunUnbalanced.sh ' + job_name + ' -c ' + config_name + ' -s ' + str(k[0]) + ' -e ' + str(k[1]) + ' > Results/' + config_name + '_' + str(k[0]) + '_' + str(k[1]) +  '.log' + '\n'
    fid.write(cmd_str)
    ijob = ijob + 1

fid.close()
