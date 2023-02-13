function combined_image = read_embryo_frame(data_path, name_of_embryo, ...
                suffix_for_embryo, ...
                suffix_for_embryo_alt, ...
                time_index)
%% Reads a TIF or KLB image and does some sampling to handle isotropy
%% Inputs: 
%%  name_of_embryo: The prefix name of the file name upto the 5 digit time_index
%%  suffix_for_embryo: The suffix of the file name after the time index 
%%  suffix_for_embryo_alternative: The suffix of the file name after the time index if the file has been hand corrected
%%  time_index: The time index

%addpath('/mnt/home/lbrown/Registration/RegisterSequenceMatlab/Abhishek/lineage_track/klb_io/');
%addpath('/mnt/home/lbrown/Registration/RegisterSequenceMatlab/Lever/leverUtilities/src/MATLAB/+MicroscopeData/+KLB/')
% used David's instructions but didn't need to make the mex (was in the make??)
addpath('/mnt/home/lbrown/Registration/RegisterSequenceMatlab/Abhishek/lineage_track/KLBsource/keller-lab-block-filetype/matlabWrapper/klb_mexa64/')
    emb_name = fullfile(data_path, strcat(name_of_embryo,num2str(time_index,'%05.5d'),suffix_for_embryo_alt));
    if ~isfile(emb_name)
        emb_name = fullfile(data_path, strcat(name_of_embryo,num2str(time_index,'%05.5d'),suffix_for_embryo));
        if ~isfile(emb_name)
            disp(emb_name)
            disp('(orig read) file does not exist, passing to next iteration');
            combined_image = zeros();
            return;
        end
    end
    if endsWith(suffix_for_embryo, 'klb')
        try
            combined_image= readKLBstack(emb_name);
            %% LB Note - this only for KLB -- not in Abhishek's code
            %combined_image = permute(combined_image, [2 1 3]);
        catch
            disp(emb_name);
            error('KLB could not be read. Try using TIF images. quitting...');
        end

%           %% temp code to convert to tiff
%           nc = size(emb_name,2);
%           tiff_name = strcat(emb_name(1:(nc-3)),'tiff');
%           disp(tiff_name)
%           combined_image = uint8(combined_image ); % lowest 8-bits
%           nslices = size(combined_image,3);
%            % First slice:
%            imwrite(im2uint8(combined_image(:,:,1)),tiff_name);
%            for islice = 2:nslices
%              imwrite(im2uint8(combined_image(:,:,islice)),tiff_name,'WriteMode','append');
%            end 
          
        %% need extrasamples or 3D ??
%         t = Tiff(tiff_name, 'w');
%         tagstruct.ImageLength = size(combined_image, 1);
%         tagstruct.ImageWidth = size(combined_image, 2);
%         tagstruct.Compression = Tiff.Compression.None;
%         tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
%         %if size(combined_image,3) == 1
%         %tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
%         %else
%         tagstruct.Photometric = Tiff.Photometric.RGB;
%         %end
%         tagstruct.BitsPerSample = 8;
%         nslices = size(combined_image,3);
%         tagstruct.SamplesPerPixel = nslices;
%         tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
%         t.setTag(tagstruct);
%         write(t,squeeze(im2uint8(combined_image(:,:,:))));
% 
%         %% write all the slices
%         %for islice=1:nslices
%         %    write(t,squeeze(im2uint8(combined_image(:,:,islice))));
%         %end
%         t.close();

    elseif endsWith(suffix_for_embryo, 'tif')||endsWith(suffix_for_embryo, 'tiff')
        A = imread(emb_name,1);
        tiff_info = imfinfo(emb_name);
        % combine all tiff stacks into 1 3D image.
        combined_image = zeros(size(A,1), size(A,2), size(tiff_info, 1));
        for j = 1:size(tiff_info, 1)
            A = imread(emb_name,j);
            combined_image(:,:,j) = A(:,:,1);
        end
    else
        error('Filename should end with tif, tiff or klb');
    end
    combined_image = permute(combined_image, [2 1 3]);
    resXY = 0.208;
    resZ = 2.0;
    reduceRatio = 1/4;
    combined_image = isotropicSample_nearest(double(combined_image), resXY, resZ, reduceRatio);
end