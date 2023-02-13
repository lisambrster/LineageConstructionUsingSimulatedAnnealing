
function [] = CheckSplits(config_name)

% output visualizations to check lineage constructions
% outputs are Output_Path/LineageVisualizations
% output matlab figs to check if lineage is correct

% output point clouds of nuclei for each of the two frames involved in splits with labels/links of matches
%   Output_Path/LineageVisualizations/Point_Clouds_FrameXX.fig (XX is the frame number)
% output centroids of nuclei for each of the two frames involved in splits with labels/links of matches
%   Output_Path/LineageVisualizations/Just_Points_FrameXX.fig
% where Output_Path is specified in configuration

    config_path = '../';
config_path = strcat(config_path,config_name);

addpath(genpath('../CPD2/core'));
addpath(genpath('../CPD2/data'));
addpath(genpath('../YAMLMatlab_0.4.3'));
addpath(genpath('../klb_io'));
addpath(genpath('../common'));

config_opts = ReadYaml(fullfile(config_path,'config.yaml'));

output_path = strcat(config_opts.output_dir,'/LineageVisualizations/');
[status, msg, msgID] = mkdir(output_path);
reg_start_frame  =  config_opts.register_begin_frame;

sim_graph_file = fullfile(config_opts.output_dir,'SimGraphs/Full_Sim_Graph.mat');
gg = load(sim_graph_file);
gg = gg.full_sim_graph;

% plot the graph
figure;
H = plot(gg,'layout','layered');

% find minimum frame
start_graph_frame = 1000;
nNodes = size(gg.Nodes,1);
for iNode = 1:nNodes
        node_str = gg.Nodes{iNode,1};
        node_frame = str2num(node_str{1,1}(1:3));
        if (node_frame <  start_graph_frame)
            start_graph_frame = node_frame
        end
end

%% get all frames with division based on graph
nSplits = 0;
d =degree(gg);
last_split_frame = -1;
nNodes = size(gg.Nodes,1);
for iNode = 1:nNodes
    deg = d(iNode);
    node_str = gg.Nodes{iNode,1};
    node_frame = str2num(node_str{1,1}(1:3));
    if (deg > 2)  | ((node_frame == start_graph_frame) & (deg == 2))
        mother_str = gg.Nodes{iNode,1};
        iframe = str2num(mother_str{1,1}(1:3));
        if (iframe ~= last_split_frame)
            nSplits = nSplits + 1;
            frames_with_splits(nSplits) = iframe;
            last_split_frame = iframe;
        end
    end
end
nSplitFrames = size(frames_with_splits,2);
%% just temp for GT graph
%%frames_with_splits = frames_with_splits(reg_start_frame:nSplitFrames);

trans_path = fullfile(config_opts.output_dir, ...
    strcat(config_opts.register_file_name_prefix,'_transforms.mat'));
% frames that need to be re-registered
register_again = [];%[10, 13];%[16];
% define threshold for registration (if necessary)
sigma_thres = 30;

% TREE SETUP
% if ground truth tree (for testing) is available, we can use it to clean labels
tree_path = '';
next_graph_file = '';
% should we consider only labels present in ground truth tree above?
clean_data = false;

% DAUGHTER VOLUME THRESHOLD SETUP
vol_thres = 1000;

% PLOTTING SETUP
plot_all = false; %% true doesn't work - don't get links on point clouds
% Set plot sizes 
plot_width = 600;
plot_height = 400;

% IMAGE INDICES
firstTime = config_opts.register_begin_frame;
lastTime =  config_opts.register_end_frame-1;
% to consider overall
which_number_vect = 1:lastTime+1;
% to use for tracking

inds_to_track = firstTime:lastTime;
%-----END_OF_MAIN_SETUP-----

% DISTANCE SETUP
% I did not find this useful so set it to very liberal values
%minimal distance between mother and daughters
dist_thres = 0;
%maximal distance between mother and daughters' centroid
dist_cent_thres = 100;


% Initialize empty graph and cell array for storing registration
G_based_on_nn = graph;
valid_time_indices = which_number_vect;
store_registration = cell((length(valid_time_indices)-1), 1);

for itest = 1: nSplitFrames

    time_index_index = frames_with_splits(itest);
    time_index = time_index_index;
    if ~plot_all
        close all;
    end

    % store next in series
    time_index_plus_1 = valid_time_indices(time_index_index+1);
    
    % store combined image for both.
    combined_image1 = read_embryo_frame(config_opts.data_path, config_opts.name_of_embryo, ...
        config_opts.suffix_for_embryo, ...
        config_opts.suffix_for_embryo_alternative, ...
        time_index+reg_start_frame);

    combined_image2 = read_embryo_frame(config_opts.data_path, config_opts.name_of_embryo, ...
        config_opts.suffix_for_embryo, ...
        config_opts.suffix_for_embryo_alternative, ...
        time_index_plus_1+reg_start_frame);
  
    % STORE MESHGRID
    [X, Y, Z] = meshgrid(1:size(combined_image1, 2), 1:size(combined_image1, 1), 1:size(combined_image1, 3));
    
    % FRACTION OF POINTS (DOWNSAMPLING)
    fraction_of_selected_points =  1/10;  % slow to run at full scale - but make full res points and xform?
    find1 = find(combined_image1(:)~=0);  % this is the indices into combined_image1 to get indices into (X,Y,Z) to the full set of point
    number_of_points = length(find1);
    
    % why random points - why not just subsample by 10 ?
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    find1 = find1(p);
    
    ptCloud1 = [X(find1), Y(find1), Z(find1)] - [mean(X(find1)), mean(Y(find1)), mean(Z(find1))];
    %
    find2 = find(combined_image2(:)~=0);
    number_of_points = length(find2);
    
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    find2 = find2(p);
    
    ptCloud2 = [X(find2), Y(find2), Z(find2)] - [mean(X(find2)), mean(Y(find2)), mean(Z(find2))];
    ptCloud2 = pointCloud(ptCloud2);

    % Calculate nuclear volumes
    volumes1 = regionprops3(combined_image1, 'Volume').Volume;
    volumes2 = regionprops3(combined_image2, 'Volume').Volume;

  
    % register the two frames using the registration transform
    transforms = load(trans_path);
    Transform = transforms.store_registration{time_index_index  + 1 +reg_start_frame, 1};  %%% LB: added + 1
    R = Transform.Rotation;
    t = Transform.Translation;
    [M, D]=size(ptCloud2.Location);
    %Transform.Y=(ptCloud2.Location-repmat(t(1,:), [M,1]))*R';
     Transform.Y=(ptCloud2.Location + repmat(t(1,:), [M,1]))*R';
     [M, D]=size(ptCloud1);
     newX = (ptCloud1)*R - repmat(t',[M,1]);


    f = figure; f.Position = [50 50 plot_width plot_height]; hold all; 
    title('After registering X to Y.'); 
    %cpd_plot_iter(ptCloud1, Transform.Y);
    %% LB swap
    cpd_plot_iter(newX, ptCloud2.Location);
    hold all;
   % End of registration
    
    % Begin tracking
    %alpha_shape_for_each_label1, alpha_shape_for_each_label2] = compute_matches_based_on_point_clouds_CPD(Transform.Y,ptCloud1,...
    [iou_matrix, M, corresponding_ious_for_matches, ...
            cell_labels_I_care_about1, cell_labels_I_care_about2, ...
            center_point_for_each_label1, center_point_for_each_label2, ...
            match_based_on_nearest_neighbors, ~, ~, ...
            alpha_shape_for_each_label1, alpha_shape_for_each_label2] = compute_matches_based_on_point_clouds_CPD(ptCloud2.Location,newX,...
            combined_image1,combined_image2,find1,find2);
% 
%     disp('Cell_labels_I_care_about1');
%     disp(cell_labels_I_care_about1);
%     disp('Cell_labels_I_care_about2');
%     disp(cell_labels_I_care_about2);

    store_iou_table{time_index_index, 1} = iou_matrix;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    color_vec = [];
    color_map_setting = [1 0 0; 0 0 1];
    sample_graph_sim = graph;
    sample_graph_gt = graph;
    missing_nodes  = cell2table(cell(0, 3), 'VariableNames', {'xpos' 'ypos' 'zpos'}); 
    for iind = 1:length(cell_labels_I_care_about1)
        this_label = cell_labels_I_care_about1(iind);
        
        % store node props table... so that node can be added with volume
        node_id = {[num2str(time_index,'%05.3d'),'_', ...
            num2str(this_label,'%05.3d')]};
        NodePropsTable = table(node_id, ...
            center_point_for_each_label1(iind, 1), ...
            center_point_for_each_label1(iind, 2), ...
            center_point_for_each_label1(iind, 3), ...
            'VariableNames',{'Name' 'xpos' 'ypos' 'zpos'});
        
        sample_graph_sim = addnode(sample_graph_sim, NodePropsTable);
        sample_graph_gt = addnode(sample_graph_gt, NodePropsTable);
     end

    for iind = 1:length(cell_labels_I_care_about2)
        this_label = cell_labels_I_care_about2(iind);
        
        % store node props table... so that node can be added with volume
        node_id = {[num2str(time_index_plus_1,'%05.3d'),'_', ...
            num2str(this_label,'%05.3d')]};
        NodePropsTable = table(node_id, ...
            center_point_for_each_label2(iind, 1), ...
            center_point_for_each_label2(iind, 2), ...
            center_point_for_each_label2(iind, 3), ...
            'VariableNames',{'Name' 'xpos' 'ypos' 'zpos'});
        
        sample_graph_sim = addnode(sample_graph_sim, NodePropsTable);
        sample_graph_gt = addnode(sample_graph_gt, NodePropsTable);

    end
    
    nNuclei1 = size(cell_labels_I_care_about1,1);
    for kk=1:nNuclei1
        color_vec = [color_vec; 1];
    end
    nNuclei2 = size(cell_labels_I_care_about2,1);
    for kk=1:nNuclei2
        color_vec = [color_vec; 2];
    end
    
    % add edges based on simulated annealing graph
    bSimGraph = true;  %% when both are drawn some stuff gets remembered..
    if (bSimGraph)
        nEdges = size(gg.Edges,1);
        for iEdge = 1:nEdges
            time1 = str2num(gg.Edges{iEdge,1}{1}(1:3));
            label1 = str2num(gg.Edges{iEdge,1}{1}(5:7));
            ind = find(cell_labels_I_care_about1 == label1);
            n1 = size(ind,1);            
            time2 = str2num(gg.Edges{iEdge,1}{2}(1:3));
            label2 = str2num(gg.Edges{iEdge,1}{2}(5:7));
            ind = find(cell_labels_I_care_about2 == label2);
            n2 = size(ind,1);
%             if (n1 == 0) | (n2  == 0)
%                 disp('size match issue');
%             end
            if (time1 == (time_index)) && (time2 == (time_index_plus_1 )) && (n1 > 0) && (n2 > 0)
                sample_graph_sim = addedge(sample_graph_sim,[num2str(time_index,'%05.3d'),'_', ...
                     num2str(label1,'%05.3d')],...
                    [num2str(time_index_plus_1,'%05.3d'),'_', ...
                    num2str(label2,'%05.3d')]);    
            end
        end    
    
        % put labels and links on point cloud plot
        hold all; plot(sample_graph_sim, 'XData', sample_graph_sim.Nodes.xpos, ...
            'YData', sample_graph_sim.Nodes.ypos, 'ZData', sample_graph_sim.Nodes.zpos, ...
            'NodeCData', color_vec, 'NodeFontSize', config_opts.marker_font_size, ...
            'EdgeColor', 'k', 'LineWidth', 2.0, 'MarkerSize', config_opts.marker_size, ...
            'Interpreter','none');
        hold all;
        colormap(color_map_setting);
        cloud_name = sprintf('Point_Clouds_Frame%0.5d.fig',time_index);
        savefig(strcat(output_path,cloud_name));

        % make new figure with just centroids and links
        f2 = figure; f2.Position = [600 50 plot_width plot_height];
        h2 = plot(sample_graph_sim, 'XData', sample_graph_sim.Nodes.xpos, 'YData', ...
            sample_graph_sim.Nodes.ypos, 'ZData', sample_graph_sim.Nodes.zpos, ...
            'NodeCData', color_vec, 'NodeFontSize', config_opts.marker_font_size, ...
            'EdgeColor', 'k', 'LineWidth', 2.0,'NodeLabel', sample_graph_sim.Nodes.Name, ...
            'MarkerSize', config_opts.marker_size, 'Interpreter','none');
    
        hold all;
        colormap(color_map_setting);
        % save point graph
        point_name = sprintf('Just_Points_Frame%0.5d.fig',time_index);
        savefig(strcat(output_path,point_name));
    end
    
    % add edges based on graph G
    bGroundTruthGraph =false;
    if (bGroundTruthGraph)
        nEdges = size(g.Edges);
        for iEdge = 1:nEdges
            e1 = g.Edges{iEdge,1}{1};
            time1 = str2num( e1(1:3));
            label1 = str2num(e1(5:7));
            e2 = g.Edges{iEdge,1}{2};
            time2 = str2num( e2(1:3));
            label2 = str2num(e2(5:7));
            if (time1 == time_index) && (time2 == time_index_plus_1)
                sample_graph_gt = addedge(sample_graph_gt,[num2str(time_index,'%05.3d'),'_', ...
                     num2str(label1,'%05.3d')],...
                    [num2str(time_index_plus_1,'%05.3d'),'_', ...
                    num2str(label2,'%05.3d')]);    
            end
        end
    end
    
    % visualization for checking if everything is correct
    bLBOff = bGroundTruthGraph;
    if (bLBOff)
        % put labels and links on point cloud plot
        hold all; plot(sample_graph_gt, 'XData', sample_graph_gt.Nodes.xpos, ...
            'YData', sample_graph_gt.Nodes.ypos, 'ZData', sample_graph_gt.Nodes.zpos, ...
            'NodeCData', color_vec, 'NodeFontSize', config_opts.marker_font_size, ...
            'EdgeColor', 'k', 'LineWidth', 2.0, 'MarkerSize', config_opts.marker_size, ...
            'Interpreter','none');
        colormap(color_map_setting);
        % make new figure with just centroids and links
        f2 = figure; f2.Position = [600 50 plot_width plot_height];
        h2 = plot(sample_graph_gt, 'XData', sample_graph_gt.Nodes.xpos, 'YData', ...
            sample_graph_gt.Nodes.ypos, 'ZData', sample_graph_gt.Nodes.zpos, ...
            'NodeCData', color_vec, 'NodeFontSize', config_opts.marker_font_size, ...
            'EdgeColor', 'k', 'LineWidth', 2.0,'NodeLabel', sample_graph_gt.Nodes.Name, ...
            'MarkerSize', config_opts.marker_size, 'Interpreter','none');
    
        hold all;
        colormap(color_map_setting);
    end

    disp('time index');
    disp(time_index);
    %pause;
    
end


end

