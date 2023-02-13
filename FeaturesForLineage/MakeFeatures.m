


% make outputs for simulated annealing

function [] = MoreFeatures( config_name )

% given  registration precomputed
% find actual number of nuclei per frame
% output centroid, volume, IoU, solidity, etc

% register frames for centroids


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear;
%close all;
%config_name = '210809_stack3'
config_path = '../';
config_path = strcat(config_path,config_name)

bVerbose =false; % show plots

addpath(genpath('../CPD2/core'));
addpath(genpath('../CPD2/data'));
addpath(genpath('../YAMLMatlab_0.4.3'));
addpath(genpath('../klb_io'));
addpath(genpath('../common'));

config_opts = ReadYaml(fullfile(config_path,'config.yaml'));

out_path = config_opts.output_dir;

image_path = config_opts.image_path;

data_path = config_opts.data_path;
name_of_embryo = config_opts.name_of_embryo;
suffix_for_embryo = config_opts.suffix_for_embryo;
suffix_for_embryo_alternative = config_opts.suffix_for_embryo_alternative;
read_type = 0; %config_opts.read_type;

firstTime = config_opts.register_begin_frame;
lastTime = config_opts.register_end_frame - 1;

RegistrationFileName =  fullfile(config_opts.output_dir, ...
    strcat(config_opts.register_file_name_prefix,'_transforms.mat'));
transforms = load(RegistrationFileName);

figure;
hold on;
firstTime1 = firstTime; % note - usually zero - need to be careful if changing in config..

for i=firstTime:(lastTime-1)
    % s and t get +1 because might start at 0
    s(i+1) = transforms.store_registration{i - firstTime1 + 1,1}.minSigma;
    t(i+1) = transforms.store_registration{i - firstTime1 + 1,1};
end

plot(firstTime:lastTime-1,s(firstTime+1:lastTime),'LineWidth',4,'Color','b');

pito1over3 = (pi)^(1/3);
% Voxel size before making isotropic
pixel_size_xy_um = 0.208; % um
pixel_size_z_um = 2.0; % um
% Voxel size after making isotropic
xyz_res = 0.8320;
% Volume of isotropic voxel
voxel_vol = xyz_res^3;

% Which image indices to run over...
which_number_vect = 1:200;
valid_time_indices = which_number_vect;
store_registration = cell((length(valid_time_indices)-1), 1);
sigma2_vect_saved = zeros((length(valid_time_indices)-1), 1);
sigma2tests = zeros((length(valid_time_indices)-1),1);  % these are the 100 tests for one pair of images

%% Note - compute these after registration
nuclei_centroids = zeros((length(valid_time_indices)-1),250,3); % for each time and nucleus
nuclei_centroids_next = zeros((length(valid_time_indices)-1),250,3); % for each adjacent time step and nucleus
nuclei_volumes = zeros((length(valid_time_indices)-1),250,1);
nuclei_boxes = zeros((length(valid_time_indices)-1),250,6);
nuclei_solidities = zeros((length(valid_time_indices)-1),250,1);
nuclei_IoUs = zeros((length(valid_time_indices)-1),250,1);
IoU_tables = cell((length(valid_time_indices)-1), 1);
label_ids = zeros((length(valid_time_indices)-1),250,1);
nuclei_sphericity = zeros((length(valid_time_indices)-1),250,1);
nuclei_ratio = zeros((length(valid_time_indices)-1),250,1);
mean_intensity = zeros((length(valid_time_indices)-1),250,1);
mean_intensity20 = zeros((length(valid_time_indices)-1),250,1);
std_intensity = zeros((length(valid_time_indices)-1),250,1);



for time_index_index = firstTime:lastTime  
     
    % store this time index
    %time_index = valid_time_indices(time_index_index)
    
    % store next in series
    %time_index_plus_1 = valid_time_indices(time_index_index+1);
    
    % store combined image for both.
    % read in segmented images and rescale 
    if (read_type == 1)
        combined_image1 = read_embryo_frame_FSEQ(data_path, name_of_embryo, ...
            suffix_for_embryo, ...
            suffix_for_embryo_alternative, ...
            time_index_index);
    else        
        combined_image1 = read_embryo_frame(data_path, name_of_embryo, ...
            suffix_for_embryo, ...
            suffix_for_embryo_alternative, ...
            time_index_index);
    end
    % note just undoing what happens in read
    %combined_image1 = permute(combined_image1, [2 1 3]);

    % read in segmented images and rescale 
    if (read_type == 1)
        combined_image2 = read_embryo_frame_FSEQ(data_path, name_of_embryo, ...
            suffix_for_embryo, ...
            suffix_for_embryo_alternative, ...
            time_index_index + 1);
    else        
        combined_image2 = read_embryo_frame(data_path, name_of_embryo, ...
            suffix_for_embryo, ...
            suffix_for_embryo_alternative, ...
            time_index_index + 1);
    end
    % note just undoing what happens in read
    %combined_image2 = permute(combined_image2, [2 1 3]);

    % STORE MESHGRID
    [X, Y, Z] = meshgrid(1:size(combined_image1, 2), 1:size(combined_image1, 1), 1:size(combined_image1, 3));
    
    % FRACTION OF POINTS (DOWNSAMPLING)
    fraction_of_selected_points =  1/10;  % slow to run at full scale - but make full res points and xform?
    nNuclei1 = length(unique(combined_image1)) -1

    for iNuc = 1:nNuclei1
        find_nuc = find(combined_image1(:) == iNuc);
    end
    find1 = find(combined_image1(:)~=0);  % this is the indices into combined_image1 to get indices into (X,Y,Z) to the full set of point
    number_of_points = length(find1);
    
    % why random points - why not just subsample by 10 ?
    rng(1);
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    find1 = find1(p);
    

    meanX1 = mean(X(find1));
    meanY1 = mean(Y(find1));
    meanZ1 = mean(Z(find1));
    ptCloud1 = [X(find1), Y(find1), Z(find1)] - [meanX1, meanY1, meanZ1];

    [X, Y, Z] = meshgrid(1:size(combined_image2, 2), 1:size(combined_image2, 1), 1:size(combined_image2, 3));

    nNuclei2 = length(unique(combined_image2)) -1
    find2 = find(combined_image2(:)~=0);
    number_of_points = length(find2);
    
    rng(1);
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    find2 = find2(p);
    
  
    meanX2 = mean(X(find2));
    meanY2 = mean(Y(find2));
    meanZ2 = mean(Z(find2));
    
    ptCloud2 = [X(find2), Y(find2), Z(find2)] - [meanX2, meanY2, meanZ2];

    Y = ptCloud2;
    X = ptCloud1;
    [Mx, D] = size(X);
    [My, D] = size(Y);
    
        
    if bVerbose
        figure; hold all; title('X registered to Y'); cpd_plot_iter(X,Y);
    end
    
    bFullReg = true;
    if bFullReg
        % perform the transformation iteratively (x1 -> x2 -> ... -> xn)
        if time_index_index == lastTime
            newX = X;
        else
            newX = X;
            for iFrame = time_index_index: lastTime -1 % last frame stays the same
                ThisTransform = t(iFrame+1);
                this_rotation = ThisTransform.Rotation;
                this_translation = ThisTransform.Translation;
                newX = newX*this_rotation - repmat(this_translation',[Mx,1]);  %%%%%%%% for pairs not newX
            end
        end

        % perform the transformation iteratively (x1 -> x2 -> ... -> xn-1) --
        % for Y
        if time_index_index >= lastTime
            newY = Y;
        else
            newY = Y;
            for iFrame = time_index_index+1: lastTime  -1% last frame stays the same
                ThisTransform = t(iFrame+1);
                this_rotation = ThisTransform.Rotation;
                this_translation = ThisTransform.Translation;
                newY = newY*this_rotation - repmat(this_translation',[My,1]);  %%%%%%%% for pairs not newX
            end
        end
    else % just perform registration between these two frames (more accurate)
         if time_index_index == lastTime
            newX = X;
        else
            newX = X;
            iFrame = time_index_index;
            ThisTransform = t(iFrame+1);
            this_rotation = ThisTransform.Rotation;
            this_translation = ThisTransform.Translation;
            newX = newX*this_rotation - repmat(this_translation',[Mx,1]);
            
            newY = Y;
         end
        
    end
    
    
    
    if bVerbose
        figure; hold all; title('X registered to Y'); cpd_plot_iter(newX,newY);
    end
   
    % Calculate nuclear volumes and solidity
    %%  the label images are sub sampled
    volumes1 = regionprops3(combined_image1, 'Volume').Volume;
    volumes2 = regionprops3(combined_image2, 'Volume').Volume;
    solidity1 = regionprops3(combined_image1,'Solidity').Solidity;
    solidity2 = regionprops3(combined_image2,'Solidity').Solidity;
    %[ulf_x ulf_y ulf_z width_x width_y width_z]
    box1 = regionprops3(combined_image1,'BoundingBox').BoundingBox;
    box2 = regionprops3(combined_image2,'BoundingBox').BoundingBox;

    principal_axes = regionprops3(combined_image1,"PrincipalAxisLength").PrincipalAxisLength;
    surface_area = regionprops3(combined_image1,'SurfaceArea').SurfaceArea;
    for in = 1: size(volumes1,1)
        sphericity(in,1) = (pito1over3 * ((6 *volumes1(in,1)) ^ (2/3)))/surface_area(in,1);
        aspectratio(in,1) = principal_axes(in,1)/principal_axes(in,3);
    end
    %% intensity based metrics
    if  endsWith(image_path, 'klb')
        image_full_path = sprintf(image_path, time_index_index,time_index_index)
        raw_image= readKLBstack(image_full_path); % x,y,z
        % why here since on label read done twice?
        %raw_image = permute(raw_image, [2 1 3]);
    elseif endsWith(image_path,'tif')
        image_full_path = sprintf(image_path, time_index_index)
        A = imread(image_full_path,1);
        tiff_info = imfinfo(image_full_path);
        % combine all tiff stacks into 1 3D image.
        raw_image = zeros(size(A,1), size(A,2), size(tiff_info, 1));
        for j = 1:size(tiff_info, 1)
            A = imread(image_full_path,j);
            raw_image(:,:,j) = A(:,:,1);
        end
        % do not use for F32/F44
        %raw_image = permute(raw_image, [2 1 3]);
    end


    % downsample to match labels and swap axes
    resXY = 0.208;
    resZ = 2.0;
    reduceRatio = 1/4;
    raw_image = isotropicSample_nearest(double(raw_image), resXY, resZ, reduceRatio);
    raw_image = permute(raw_image,[2 1 3]);

    disp(size(raw_image));
    disp(size(combined_image1));
    [d,h,w] = size(raw_image);
    [d1,h1,w1] = size(combined_image1);
    combined_image1 = combined_image1(1:d,1:h,1:w);
    disp(size(combined_image1));

    % find mean intensity of each label
    %MeanIntensity = regionprops3(combined_image1, raw_image, 'MeanIntensity')
    %MinIntensity = regionprops3(combined_image1, raw_image, 'MinIntensity');
    %MaxIntensity = regionprops3(combined_image1, raw_image, 'MaxIntensity');
    %Voxels = regionprops3(combined_image1, raw_image,"VoxelValues");
    %for iNuc = 1: nNuclei1
    %    StdIntensity(iNuc) = single(std(Voxels(iNuc,1).VoxelValues{:,1}));
    %end
    %StdIntensity = single(StdIntensity);
    %disp('std')
    %disp(class(StdIntensity))
    %disp(size(StdIntensity))
    %disp(StdIntensity)
    % original meanI and stdI - but get top 20%
    nlabels = length(unique(combined_image1)) - 1;
    labels = unique(combined_image1);
    ind = find(combined_image1 == 0);
    meanIntensity0 = mean(raw_image(ind));

    nbins = 1000;
    for ilabel = 0:nlabels
        this_label = labels(ilabel + 1);
        ind = find(combined_image1 == this_label);
        [counts, centers] = hist(raw_image(ind),nbins);
        % get total of all counts - same as size of ind
        total_counts = size(ind,1);
        % find last bin before including 80%
        cum_count = 0;
        last_bin = nbins;
        ibin = 1;
        while  (last_bin == nbins)
            cum_count = cum_count + counts(ibin);
            if (cum_count/total_counts > 0.8)
                last_bin = ibin-1;
            end
            ibin = ibin + 1;
        end
        % get mean of all the intensities after last_bin
        meanI = 0;
        meanCount = 0;
        for ibin = last_bin:nbins
            meanCount = meanCount + counts(ibin);
            meanI = meanI + counts(ibin)*centers(ibin);
        end
        if (ilabel ~= 0)
            meanIntensity20(this_label) = meanI / meanCount;
            meanIntensity(this_label) = mean(raw_image(ind)) - meanIntensity0;
            stdIntensity(this_label) = std(raw_image(ind));
        else
            meanI0 = meanI / meanCount;
        end
    end
    %disp(meanIntensity20)
    %disp(meanIntensity)
    %disp(stdIntensity)

    % object ordered by labels (what if missing label) - then NaNs

    [volumes_for_each_label1, volumes_for_each_label2, iou_matrix, M, corresponding_ious_for_matches, ...
            cell_labels_I_care_about1, cell_labels_I_care_about2, ...
            center_point_for_each_label1, center_point_for_each_label2, match_based_on_nearest_neighbors, alpha_shape_for_each_label1, alpha_shape_for_each_label2 ] = compute_matches_based_on_point_clouds_CPD_more( newY,newX, ...
            combined_image1,combined_image2,find1,find2);
     store_matches{time_index_index-firstTime1+1, 1} = M;
     store_iou_table{time_index_index-firstTime1+1, 1} = iou_matrix;
     IoU_tables{time_index_index-firstTime1+1, 1} = iou_matrix;

     disp('centroids ')
     disp(center_point_for_each_label1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % make the graph..
   bMakeGraph = false;
   if bMakeGraph
    [nn nd]=kNearestNeighbors(center_point_for_each_label1, center_point_for_each_label2,min(3,length(center_point_for_each_label2)));
    if length(nn(:,1))~=length(unique(nn(:,1))) % Reject duplicate nearest neighbors
        dup=find_duplicates(nn(:,1));
        for lvd=1:size(dup,1)
            [ic,ia,ib]=intersect(nn(dup(lvd).ind,2),setdiff(1:size(nn,1),nn(:,1)));
            if ~isempty(ia)
                nn(dup(lvd).ind(ia),1)=nn(dup(lvd).ind(ia),2);
                nd(dup(lvd).ind(ia),1)=nd(dup(lvd).ind(ia),2);
                loi=dup(lvd).ind(setdiff(1:size(dup(lvd).ind,1),ia)); % treat triple and more entries
                if length(loi)>1
                    [mv,mi]=min(nd(loi,1));
                    loi1=setdiff(1:length(loi),mi);
                    nn(loi(loi1),1)=NaN;
                end
            else
                [mv mi]=min(sum(nd(dup(lvd).ind,1),2));
                loi=setdiff(dup(lvd).ind,dup(lvd).ind(mi));
                nn(loi,1)=NaN;
            end
        end
    end
    
    nn=nn(:,1);
    nd=nd(:,1);
    
    sample_graph = graph;
    for iind = 1:length(cell_labels_I_care_about1)
        this_label = cell_labels_I_care_about1(iind);
        
        % store node props table... so that node can be added with volume
        NodePropsTable = table({[num2str(time_index,'%05.3d'),'_', num2str(this_label,'%05.3d')]}, center_point_for_each_label1(iind, 1), center_point_for_each_label1(iind, 2), center_point_for_each_label1(iind, 3), ...
            'VariableNames',{'Name' 'xpos' 'ypos' 'zpos'});
        
        sample_graph = addnode(sample_graph, NodePropsTable);
    end
    
    for iind = 1:length(cell_labels_I_care_about2)
        this_label = cell_labels_I_care_about2(iind);
        
        % store node props table... so that node can be added with volume
        NodePropsTable = table({[num2str(time_index_plus_1,'%05.3d'),'_', num2str(this_label,'%05.3d')]}, center_point_for_each_label2(iind, 1), center_point_for_each_label2(iind, 2), center_point_for_each_label2(iind, 3), ...
            'VariableNames',{'Name' 'xpos' 'ypos' 'zpos'});
        
        sample_graph = addnode(sample_graph, NodePropsTable);
    end
    
    for point_index = 1:length(nn)
        
        if (~isnan(nn(point_index)))
            % make directed edges (in time) between matches + store iou for the match as a graph weight
            sample_graph = addedge(sample_graph, [num2str(time_index,'%05.3d'),'_', num2str(cell_labels_I_care_about1(nn(point_index)),'%05.3d')],...
                [num2str(time_index_plus_1,'%05.3d'),'_', num2str(cell_labels_I_care_about2(point_index),'%05.3d')]);
        end
%         if (~isnan(nn(point_index)))
%             
%             % make directed edges (in time) between matches + store iou for the match as a graph weight
%             G_based_on_nn = addedge(G_based_on_nn, [num2str(time_index,'%05.3d'),'_', num2str(cell_labels_I_care_about1(nn(point_index)),'%05.3d')],...
%                 [num2str(time_index_plus_1,'%05.3d'),'_', num2str(cell_labels_I_care_about2(point_index),'%05.3d')]);
%             
%         end
        
    end
    % visualization for checking if everything is correct - 3d plot of
    % edges and nodes
    hold all; plot(sample_graph, 'XData', sample_graph.Nodes.xpos, 'YData', sample_graph.Nodes.ypos, 'ZData', sample_graph.Nodes.zpos, 'EdgeColor', 'k', 'LineWidth', 2.0);
    figure; h1 = plot(sample_graph, 'XData', sample_graph.Nodes.xpos, 'YData', sample_graph.Nodes.ypos, 'ZData', sample_graph.Nodes.zpos, 'EdgeColor', 'k', 'LineWidth', 2.0,'NodeLabel',sample_graph.Nodes.Name);
    disp(time_index);
        
    % loop through all matches
    % LB: M is nMatches x 2 (label at time t, label at time t+1)  

    for i = 1:length(M)
        find_non_zero_please = find(iou_matrix(M(i,1),:));
        if (length(find_non_zero_please) > 1)  % more than one iou match
            figure; hold all; cpd_plot_iter(newX, ptCloud2.Location);
            hold all; plot(alpha_shape_for_each_label1{M(i,1),1},'FaceColor','red','FaceAlpha',0.5);
            for j = 1:length(find_non_zero_please)
                if (find_non_zero_please(j) == M(i,2)) % this is the best match?
                    hold all; plot(alpha_shape_for_each_label2{M(i,2),1},'FaceColor','green','FaceAlpha',0.5);
                else
                    hold all; plot(alpha_shape_for_each_label2{find_non_zero_please(j),1},'FaceColor','black','FaceAlpha',0.5);
                end
                
            end
            
            title([num2str(corresponding_ious_for_matches(i)),';',num2str(i)]);
        end
    end
    

    %% this is good place for breakpoint to check that all nuclei are matched
    figure;
    plot(G_based_on_nn,'layout','layered')
    disp('time index');
    disp(time_index);
   end % make graph
    
    % compute the centroids/volumes/IoU (post registration) - in time pairs
    % output centroids/volumes/IoU pairs (nuclei_centroids,
    % nuclei_centroids_next, nuclei_volumes, nuclei_IoUs)
    % Note: this is isotropic, downsampled ..
    for iind = 1:length(cell_labels_I_care_about1)
        ilabel = cell_labels_I_care_about1(iind);
        if (ilabel ~= 999)
            %labels_id(time_index_index-firstTime1+1,ilabel) = cell_labels_I_care_about1
            %nuclei_volumes(time_index_index,ilabel,1) = volumes_for_each_label1(iind);  
            nuclei_volumes(time_index_index-firstTime1+1,ilabel,1) = volumes1(ilabel); 
            nuclei_solidities(time_index_index-firstTime1+1,ilabel,1) = solidity1(ilabel);
            nuclei_boxes(time_index_index-firstTime1+1,ilabel,:) = box1(ilabel,:);
            nuclei_ratio(time_index_index-firstTime1+1,ilabel,:) = aspectratio(ilabel,:);
            nuclei_sphericity(time_index_index-firstTime1+1,ilabel,:) = sphericity(ilabel,:);
            mean_intensity(time_index_index-firstTime1+1,ilabel,1) = meanIntensity(ilabel);
            mean_intensity20(time_index_index-firstTime1+1,ilabel,1) = meanIntensity20(ilabel);
            std_intensity(time_index_index-firstTime1+1,ilabel,1) = stdIntensity(ilabel);
            %% do these need to be registered??
            nuclei_centroids(time_index_index-firstTime1+1,ilabel,:) = center_point_for_each_label1(iind,:);
            %  if there is a match, get the IoU (**** output all IoUs - the
            %  whole table
            nMatches = size(M,1);
            nuclei_IoUs(time_index_index-firstTime1+1,ilabel) = 0;
            for iMatch=1:nMatches
                if (M(iMatch,1) == ilabel)
                    nuclei_IoUs(time_index_index-firstTime1+1,ilabel) = corresponding_ious_for_matches(iMatch);
                end
            end
        end
    end
    
    for iind = 1:length(cell_labels_I_care_about2)
        ilabel = cell_labels_I_care_about2(iind);
        %nuclei_volumes(time_index_index,ilabel,1) = volumes_for_each_label1(ilabel);  
        %% do these need to be registered??
        if (ilabel ~= 999)
            nuclei_centroids_next(time_index_index-firstTime1+1,ilabel,:) = center_point_for_each_label2(iind,:);
        end
    end
        

    %pause;
    close all;
end

% Save vector of transformations...
%save('transforms.mat', 'store_registration');
%save('matches.mat','store_matches');
%save('/Users/lbrown/Documents/PosfaiLab/Visualizations/Gata6NanogJan22/iou_table.mat','store_iou_table');
%save('graph.mat','G_based_on_nn');
info.centroids = nuclei_centroids;
info.centroids_next = nuclei_centroids_next;
info.volumes = nuclei_volumes;
info.boxes = nuclei_boxes;
info.solidities = nuclei_solidities;
info.IoU_tables = IoU_tables;
info.aspectratios = nuclei_ratio;
info.sphericity = nuclei_sphericity;
info.mean_intensity = mean_intensity;
info.mean_intensity20 = mean_intensity20;
info.std_intensity = std_intensity;
%info.IoUs = nuclei_IoUs;
%info_file = 'info_withIoUtable.mat';
%info_full_file = strcat(config_path,info_file);
%save(info_full_file,'info');

jinfo = jsonencode(info);
fid = fopen(strcat(out_path,'Features.json'),'w');
fprintf(fid, jinfo);
fclose(fid);

end

