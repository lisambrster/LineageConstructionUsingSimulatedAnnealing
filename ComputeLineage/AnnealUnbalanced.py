# main simulated annealing for frame segment with balanced splits
# input start/end and config file

import math
import random
import numpy as np

import sim_anneal
from sim_anneal_utils import *
from NewCostFunction import  DisplayAllSplits, MyCostParams, GeneralCostFunc, DisplayTree
from EvaluateCost import *
import json
import os
import yaml
import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_name', type=str,  default='.',
                        help='configuration name')
    parser.add_argument('-s', '--start_frame', type=int,  default='0',
                        help='first frame for annealing')
    parser.add_argument('-e', '--end_frame', type=int,  default='10',
                        help='last frame for annealing')
    args = parser.parse_args()

start_frame = args.start_frame
end_frame = args.end_frame
nframes = end_frame - start_frame + 1

config_path = '../';
config_name = args.config_name
print('Config ',config_name)
Output_info = {}
Output_info['config_name'] = config_name
with open(os.path.join(config_path,config_name,'config.yaml'), 'r') as file:
    config_opts = yaml.safe_load(file)
#print(config_opts)

register_start_frame = 0 #config_opts['register_begin_frame'] # used in InitGT (so need to set start_frame accordingly)
data_path = config_opts["output_dir"]
output_path = data_path
output_path = os.path.join(data_path,'SimGraphs')
if (not os.path.exists(output_path)):
    os.mkdir(output_path)
feature_path = os.path.join(data_path,'Features.json')
fid = open(feature_path, "r")
data = json.load(fid)

# make NucCounts (sequence of number of nuclei in this frame segment)
nucCounts = []
list_nuclei_labels = {} # key is frame number
for iframe in range(start_frame, start_frame + nframes):
    s = np.asarray(data['volumes'][iframe])
    ind = np.nonzero(s)
    nNuc = int(ind[0].shape[0])
    nucCounts.append(nNuc)
    labels =  list(ind[0])
    #print(len(labels),nNuc)
    list_nuclei_labels[iframe] = labels

print('lineage start_frame, nframes', start_frame, nframes)
print('Nuclear Counts ',nucCounts)

print('constructing tracks ...',start_frame ,start_frame + nframes)
tracks = sim_anneal.Tracks (nucCounts, supportSplit = False)

print ('Load features into tracks')
tracksInit (tracks, data, start_frame, list_nuclei_labels,register_start_frame)
sim_anneal.tracksShuffle (tracks)
#tracks.validate() # only if balanced

nepochs = 2000
Output_info = {}
Output_info['nepochs'] = nepochs
print('number of epochs ',nepochs)

splitWt = 0
solWt = 170 # this is symmetry
angWt = 575 #
symWt = 100 #  this is aspect ratio
meanIWt = 10
stdIWt = 0

centWt = 1; centNoSplitWt = 8; centSplitMDWt = 0; centSplitDDWt =  0;
volWt = 1; volNoSplitWt = 0.0005; volSplitWt = 0.002; # was .002
volNoSplitMult = 0.9; volSplitMult = 0.8
MyCP = MyCostParams(solWt,angWt,symWt,meanIWt,stdIWt,   centWt,centNoSplitWt,centSplitMDWt,centSplitDDWt,17,33, volWt,volNoSplitWt,volSplitWt,volNoSplitMult,volSplitMult,splitWt)
print('Unbalanced Weights ',solWt,angWt,symWt,meanIWt,stdIWt,   centWt,centNoSplitWt,centSplitMDWt,centSplitDDWt,17,33, volWt,volNoSplitWt,volSplitWt,volNoSplitMult,volSplitMult,splitWt)


cost = sim_anneal.tracksCost (tracks, daughterCostFunc = GeneralCostFunc(MyCP), splitWt = splitWt)
print('cost before annealing ',cost)

history = sim_anneal.anneal (tracks, epochs = nepochs, daughterCostFunc = GeneralCostFunc(MyCP), splitWt = None ) #, costParams = CP, history = True)   # anneal and collect history
cost = sim_anneal.tracksCost (tracks, daughterCostFunc = GeneralCostFunc(MyCP), splitWt = splitWt)
print('cost after annealing ',cost)

# display all splits
print('Stats of final lineage')
DisplayAllSplits(tracks,nframes, start_frame)

Output_info['Config_Params'] = MyCP.__dict__
Output_info['Sim_full_cost'] = cost

costs = EvaluateSimCost(tracks, MyCP, splitWt)
Output_info['costs'] = costs

# output each lineage for input into matlab based on tracks
mat_graph = OutputGraph(tracks, start_frame)
start_str = sprintf('%03d',start_frame)
end_str = sprintf('%03d',end_frame)
fid = open(os.path.join(output_path,'sim_graph_' + start_str + '_' + end_str + '.json'),'w')
json.dump(mat_graph,fid, indent = 4)
fid.close()
#print(mat_graph['Edges'])
