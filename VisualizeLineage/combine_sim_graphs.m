
function [] = view_frame_match(config_name)
% This is Masha's version of tracking, to run on early mouse embryos.
% This code attempts to construct the full lineage tree, 
% in particular, identify mitotic events.
% Pre-registration is advised but not necessary.
close all;
% SETUP starts here.
config_path = '../';
config_path = strcat(config_path,config_name);

addpath(genpath('../CPD2/core'));
addpath(genpath('../CPD2/data'));
addpath(genpath('../YAMLMatlab_0.4.3'));
addpath(genpath('../klb_io'));
addpath(genpath('../common'));

config_opts = ReadYaml(fullfile(config_path,'config.yaml'));
reg_start_frame = config_opts.register_begin_frame

% get all files in SimGraph subdir
out_dir = config_opts.output_dir;
sim_dir = strcat(out_dir,'SimGraphs/');
graph_list = dir(sim_dir)
ngraphs = size(graph_list,1)
full_sim_graph = graph;
for igraph = 3:ngraphs
    sim_name = graph_list(igraph).name;
    disp(sim_name);
end


for igraph = 3:ngraphs
    sim_name = graph_list(igraph).name;
    disp(sim_name);
    sim_graph_file = fullfile(sim_dir,sim_name);
    % check that this a file name with sim_graph_XX_XX.json
    [filepath,name,ext] = fileparts(sim_graph_file)
    if strcmp(ext,'.json')
        disp('json file')
        if strcmp(sim_name(1:8),'sim_grap')
            disp('processing ')
            fid = fopen(sim_graph_file);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            sim_data = jsondecode(str);
            val = sim_data;

            %% turn into matlab graph
            gg1 = graph;
            nNodes = size(val.Nodes,1);
            nEdges = size(val.Edges,1);
            start_graph_frame = 1000;
            end_graph_frame = 0;
            for iNode=1:nNodes
                frame = val.Nodes(iNode,1);
                if ((frame + reg_start_frame) < start_graph_frame)
                    start_graph_frame = frame + reg_start_frame;
                end
                if (frame > end_graph_frame)
                    end_graph_frame = frame + reg_start_frame;
                end
                label = val.Nodes(iNode,2);
                frame_str = num2str(frame + reg_start_frame,'%05.3d') ;
                label_str = num2str(label,'%05.3d');
                node_str = strcat(frame_str,'_',label_str);
                gg1 = gg1.addnode(node_str);
            end


    
            for iEdge=1:nEdges
                frame = val.Edges(iEdge,1,1) + reg_start_frame;
                label = val.Edges(iEdge,1,2);
                frame_str = num2str(frame,'%05.3d');
                label_str = num2str(label,'%05.3d');
                node_str1 = strcat(frame_str,'_',label_str);
                frame = val.Edges(iEdge,2,1)  + reg_start_frame;
                label = val.Edges(iEdge,2,2);
                frame_str = num2str(frame,'%05.3d');
                label_str = num2str(label,'%05.3d');
                node_str2 = strcat(frame_str,'_',label_str);
                gg1 = gg1.addedge(node_str1,node_str2);
            end

            % remove duplicate edges
            gg = graph;
            nNodes = size(gg1.Nodes);
            for iNode = 1:nNodes
                gg = gg.addnode(gg1.Nodes{iNode,1});
            end
            %newgg.Nodes = gg.Nodes;
            nEdges = size(gg1.Edges,1);
            for iEdge=1:2:nEdges
                node1 = gg1.Edges{iEdge,1}(1,1);
                node2 = gg1.Edges{iEdge,1}(1,2);
                gg = gg.addedge(node1,node2);
            end
            nEdges = size(gg.Edges,1);

            % plot the graph
            figure;
            H = plot(gg,'layout','layered');

            nEdges = size(gg.Edges)
            nNodes = size(gg.Nodes)
            
            n = size(sim_name,2);
            sim_name_out = sim_name(1:(n-5));
            disp(sim_name_out);

            % save the plot of the graph
            savefig(strcat(sim_dir,sim_name_out,'.fig'));
            % save as mat file
            save(strcat(sim_dir,strcat(sim_name_out),'.mat'),"gg");

            % combine with previous graph
            full_sim_graph = CombineGraphs(full_sim_graph,gg);

        end
    end
    disp('get next file')
end

% plot full graph and save as mat and fig
figure;
H = plot(full_sim_graph,'layout','layered');  

% add the split labels
start_frame = 0;
d =degree(full_sim_graph);
[s,t] = findedge(full_sim_graph); % all source and targets for graph
index = 1;
% find every node with only one link (it is a 'start' node) - label it
nNodes = size(full_sim_graph.Nodes);
for iNode = 1:nNodes
    bfound = false;
    deg = d(iNode);
    iframe = str2double( full_sim_graph.Nodes{iNode,1}{1,1}(1:3) );
    if  (iframe == 0)
        labelnode(H,iNode,full_sim_graph.Nodes{iNode,1});
        bfound = true;
    end
    % if start of graph has a splits (from 4 to 8)
    if (deg == 2) & (iframe == start_frame)
        labelnode(H,iNode,full_sim_graph.Nodes{iNode,1});
        bfound = true;
    end
        
    if (deg > 2)
        labelnode(H,iNode,full_sim_graph.Nodes{iNode,1});
        bfound = true;
    end
    
end 

nEdges = size(gg.Edges)
nNodes = size(gg.Nodes)
% save the fig
savefig(strcat(sim_dir,'Full_Sim_Graph.fig'));
% save the graph - also as json
save(strcat(sim_dir,'Full_Sim_Graph.mat'),'full_sim_graph');
 % now make the json version
jH = jsonencode(full_sim_graph);
json_FileName = strcat('Full_Sim_Graph.json');
fid = fopen(strcat(sim_dir,json_FileName),'w');
fprintf(fid, jH);
fclose(fid);

end
