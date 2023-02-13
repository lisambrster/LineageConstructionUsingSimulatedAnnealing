# Lineage Construction Using Simulated Annealing
matlab and python code to register, extract features and compute lineage using simulated annealing



Lineage Construction

(1)

-- make a configuration directory in the LineageConstruction directory
mkdir YourConfigName
-- make a configuration file config.yaml (see example: 220827_stack1)
-- and copy it to the CurrentConfig directory
cp config.yaml ../CurrentConfig/.


(2)

cd RegisterForLineage
sbatch -p ccb RunRegistration.sh
-- registration usually takes about 10-15 minutes
-- this will save the registration transforms in your output_dir using your prefix (from your config.yaml)

-- optionally, you can make a plot of the sigma error in registration and a video of the registration
sbatch -p ccb RunViz.sh
-- NOTE: Currently you need to run on this on a local machine (or maybe with VNC display?)


(3)

cd ../FeaturesForLineage
sbatch -p ccb RunFeatures.sh
-- this can take a few hours since it needs to read each of your label and raw images
-- this will make the Features.json file in your Output_Path

(4)

cd ../ComputeLineage
python MakeTasks.py
-- this makes taskfile.txt (in this directory)
-- Note if any frame has less nuclei than the frame before it, then processing stops at that frame
module add disBatch/beta
sbatch -p ccb -n 8 --ntasks-per-node 1 disBatch --logfile ./disBatchLogs/log.txt taskfile.txt
-- this runs simulated annealing on different sets of frames in your embryo
-- this runs quickly for no-split segments (<5min), and more slowly for segments with splits depending on
-- the number of frames (15-30min)
-- the outputs are in your Output_Path/SimGraphs

-- a little clean up now is good
rm slurm*
rm *disB*

(5)

cd ../VisualizeLineage
sbatch -p ccb RunGraphLineage.sh
-- this makes the full lineage graph here:
    Output_Path/LineageVisualizations/Full_Sim_Graph.m -- matlab formatted version of lineage graph
    Output_Path/LineageVisualizations/Full_Sim_Graph.json -- json formatted version of lineage graph
    Output_Path/LineageVisualizations/Full_Sim_Graph.fig -- matlab fig of lineage for visualization


sbatch -p ccb RunCheckSplits.sh
-- output point clouds of nuclei for each of the two frames involved in splits with labels/links of matches
Output_Path/LineageVisualizations/Point_Clouds_FrameXX.fig (XX is the frame number)
-- output centroids of nuclei for each of the two frames involved in splits with labels/links of matches
Output_Path/LineageVisualizations/Just_Points_FrameXX.fig (XX is the frame number)
-- where Output_Path is the output_dir in your config file
