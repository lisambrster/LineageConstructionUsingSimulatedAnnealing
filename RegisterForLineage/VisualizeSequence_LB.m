
% for sequence of label images
% using sequence of registration transforms
% put all label images in same reference frame

function [] = VisualizeSequence(config_name)

config_path = '../';
config_path = strcat(config_path,config_name)

%% %%%%% NO CHNAGES BELOW %%%%%%%
addpath(genpath('../CPD2/core'));
addpath(genpath('../CPD2/data'));

addpath(genpath('../YAMLMatlab_0.4.3'));
addpath(genpath('../klb_io'));
addpath(genpath('../common'));

config_opts = ReadYaml(fullfile(config_path,'config.yaml'));

firstTime = config_opts.register_begin_frame;
lastTime =  config_opts.register_end_frame-1;

RegistrationFileName = fullfile(config_opts.output_dir, ...
    strcat(config_opts.register_file_name_prefix,'_transforms.mat'));
transforms = load(RegistrationFileName);

for i=firstTime:lastTime
    s(i + 1) = transforms.store_registration{i-firstTime + 1,1}.minSigma;
     t(i + 1) = transforms.store_registration{i - firstTime + 1,1};
end
gcf1 = figure;
plot(firstTime:lastTime-1,s(firstTime+1:lastTime),'LineWidth',4,'Color','b');
xlabel('Frame');
ylabel('Registration Sigma ');
% save the plot to output dir
saveas(gcf1,strcat(config_opts.output_dir ,'RegistrationSigma.png'));


bMakeVideo = true;
%figure; % use this to look at pairs
if bMakeVideo
    % create the video writer with 1 fps
    writerObj = VideoWriter(strcat(config_opts.output_dir,'RegisteredPoints.avi'));
    writerObj.FrameRate = 10;
    % set the frames per second
    % open the video writer
    open(writerObj);
end

% Voxel size before making isotropic
pixel_size_xy_um = 0.208; % um
pixel_size_z_um = 2.9; % um
% Voxel size after making isotropic
xyz_res = 0.8320;
% Volume of isotropic voxel
voxel_vol = xyz_res^3;

% Initialize empty graph and cell array for storing registration
% Which image indices to run over...
which_number_vect = 1:lastTime;
valid_time_indices = which_number_vect;

for time_index_index = firstTime:lastTime-1
     
    % store this time index
    time_index = time_index_index;
    
    % store next in series
    time_index_plus_1 = valid_time_indices(time_index_index+1);
    
    % store combined image for both.
    combined_image1 = read_embryo_frame(config_opts.data_path, ...
            config_opts.name_of_embryo, ...
            config_opts.suffix_for_embryo, ...
            config_opts.suffix_for_embryo_alternative, ...
            time_index);
    %combined_image1 = permute(combined_image1, [2 1 3]);
    combined_image2 = read_embryo_frame(config_opts.data_path, ...
            config_opts.name_of_embryo, ...
            config_opts.suffix_for_embryo, ...
            config_opts.suffix_for_embryo_alternative, ...
            time_index_plus_1);
    %combined_image2 = permute(combined_image2, [2 1 3]);
    
    % STORE MESHGRID
    [X, Y, Z] = meshgrid(1:size(combined_image1, 2), 1:size(combined_image1, 1), 1:size(combined_image1, 3));
    
    % FRACTION OF POINTS (DOWNSAMPLING)
    fraction_of_selected_points =  1/10;  % slow to run at full scale - but make full res points and xform?
    find1 = find(combined_image1(:)~=0); 
    number_of_points = length(find1);
        
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    find1 = find1(p);
    
    ptCloud1 = [X(find1), Y(find1), Z(find1)] - [mean(X(find1)), mean(Y(find1)), mean(Z(find1))];
   
    [X, Y, Z] = meshgrid(1:size(combined_image2, 2), 1:size(combined_image2, 1), 1:size(combined_image2, 3));

    find2 = find(combined_image2(:)~=0);
    number_of_points = length(find2);
    
    p = randperm(number_of_points,round(number_of_points * fraction_of_selected_points));
    find2 = find2(p);
    
    ptCloud2 = [X(find2), Y(find2), Z(find2)] - [mean(X(find2)), mean(Y(find2)), mean(Z(find2))];
    ptCloud2 = pointCloud(ptCloud2);
    
    X = ptCloud1;
    
    bPairs = false;
    if bPairs
            if time_index_index == (lastTime-1 )
                newX = X;
            else
                newX = X;
                Transform = transforms.store_registration{time_index_index,1};
                R = Transform.Rotation;
                t = Transform.Translation;
                [M, D]=size(ptCloud2.Location);
                Transform.Y=(ptCloud2.Location + repmat(t(1,:), [M,1]))*R';
                [M, D]=size(ptCloud1);
                newX = X*R - repmat(t,[M,1]);
            end
    else
        % perform the transformation iteratively (x1 -> x2 -> ... -> xn)

        if time_index_index ==(lastTime-1 )
            newX = X;
        else
            newX = X;
            for iFrame = time_index_index: lastTime - 2 % last frame stays the same
                ThisTransform = transforms.store_registration{iFrame-firstTime+1,1};
                this_rotation = ThisTransform.Rotation;
                this_translation = ThisTransform.Translation;
                [M, D]=size(newX);
                newX = newX*this_rotation -  repmat(this_translation',[M,1]);  %%%%%%%% for pairs not newX
            end
            if time_index_index == (firstTime)
                prev1X = newX;
                prev2X = newX;
            end
        end
    end
    % for no registration
    %newX = X;
    
    none = [];
%     figure; hold all;
%     title_str = strcat({'Registering: '},string(time_index_index),{' to '},string(time_index_plus_1));
%     title(title_str); 
%     view(45,45);
    
    figure;
    hold all;
    title_str = strcat ( string(time_index_index), '.jpg');
    title(title_str); 
    xlim([-60,80]);
    ylim([-60,80]);
    zlim([-60,80]);
    view(45,45);

    if bPairs  % video of each pair
        if (time_index_index == lastTime)
            cpd_plot_iter(X,X);
        else
            %cpd_plot_iter(newX,Transform.Y);
            cpd_plot_iter(newX, ptCloud2.Location);
        end
    else
        labels = unique(combined_image1);
        nNuclei = size(labels,1) - 1 
        for ilabel=1:nNuclei
            ilabel = labels(ilabel);
            if (ilabel ~= 0)
                ind = find(newX == ilabel);
            end
        end
     
         hold on;
         %plot3(prev2X(:,1),prev2X(:,2),prev2X(:,3),'g.');
         plot3(prev1X(:,1),prev1X(:,2),prev1X(:,3),'c.');   
         plot3(newX(:,1),newX(:,2),newX(:,3),'b.');
         
         prev2X = prev1X;
         prev1X = newX;
        % get cdx2 val
        %ind = find(combined_image1
        % this is to just plot one at a time per image for making the video
%         if (mod(time_index_index,2) == 0) & (time_index_index <1000)
%             plot3(newX(:,1),newX(:,2),newX(:,3),'b.');
%         elseif (time_index_index < 10000)
%             plot3(newX(:,1),newX(:,2),newX(:,3),'b.');  % set to red to look at pairs
%         end
    end

    F(time_index_index) = getframe;
    if (bMakeVideo)
        writeVideo(writerObj, F(time_index_index));
    end
    %pause;
end

close all;
close(writerObj);
end
