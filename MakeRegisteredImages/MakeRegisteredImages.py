# register raw images
# register label images (from stardist and/or correction tool)
# and make maximum intensity projection video

# inputs - from config files
# raw image path and format
# label image path and format
# registration transforms (json) 
# output path (by default label path)


import numpy as np
import os
import math
from scipy.ndimage import affine_transform
from scipy.ndimage import zoom
import tifffile as tiff
import glob
import argparse
import json
import pyklb
import cv2
from GeometricInterpolation import InterpolateTubes
import h5py
import yaml

########################################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_name', type=str, default='config.yaml',
                        help='the name of the config file in CurrentConfig directory')
    parser.add_argument('-l', '--label_path', type=str,  default='.',
                        help='the path of the 3D label images')
    parser.add_argument('-i', '--image_path', type=str,  default='.',
                        help='the path of the 3D raw images')
    parser.add_argument('-t', '--transform_filename', type=str,  default='transforms.json',
                        help='the filename of the json transform file')
    parser.add_argument('-o', '--output_path', type=str,
                        help='the output path - for registered raw and label images')
    parser.add_argument('-s', '--start', type=int,  default=None,
                        help='start of sequence - integer')
    parser.add_argument('-e', '--end', type=int,  default=None,
                        help='end of sequence - integer')
    args = parser.parse_args()

config_path = '/mnt/home/lbrown/LineageConstruction/CurrentConfig/';

config_name = args.config_name

print('Config ',config_name)
with open(os.path.join(config_path,config_name), 'r') as file:
    config_opts = yaml.safe_load(file)

# additional argument - for on/off labels, membrane (with path) etc
#img_type = 'membrane'
img_type = 'nuclei'

label_path = config_opts["data_path"]
print(label_path)


name_of_embryo = config_opts['name_of_embryo']
suffix_for_embryo = config_opts['suffix_for_embryo']
suffix_for_embryo_alternative = config_opts['suffix_for_embryo_alternative']
# for F Sequences
#label_path = label_path + '/' + '%d' + name_of_embryo + suffix_for_embryo
numstr = '%05d'
label_path = label_path + '/' + name_of_embryo + numstr + suffix_for_embryo
print('label_path ', label_path)

out_path = config_opts["output_dir"]

if not os.path.exists(out_path):
    os.mkdir(out_path)
print('out_path ',out_path)

image_path = config_opts['image_path']

start_frame = config_opts['register_begin_frame']
end_frame = config_opts['register_end_frame']
print('start, end ',start_frame,end_frame)
transform_path = out_path #os.path.dirname(label_path)
transform_name = config_opts['register_file_name_prefix'] + '_transforms.json'

# read in registration transforms
embryo_path = transform_path
fid = open(os.path.join(embryo_path,transform_name),'r')
transforms = json.load(fid)

new_image_path = os.path.join(out_path,'registered_images')
if not os.path.exists(new_image_path):
    os.mkdir(new_image_path)
new_label_path = os.path.join(out_path,'registered_label_images')
if not os.path.exists(new_label_path):
    os.mkdir(new_label_path)
MIP_path = os.path.join(out_path,'MIP-Frames')
if not os.path.exists(MIP_path):
    os.mkdir(MIP_path)

bLabels = True # set to false if only want to register raw images
nPercents = image_path.count('%')
#crop_offset = [0,200,700] # z,y,x
#crop_box = [:64,200:1000,700:1700]

# for 220827_stack1 ONLY - REMOVE THIS!!!!
#start_frame = 95

for iframe in range(start_frame,end_frame): #
    frame_str = '%d' % iframe
    if (nPercents == 1):
        image_fullname = image_path % (iframe)
    else:
        image_fullname = image_path % (iframe,iframe)
    print(image_fullname)
    suffix = image_fullname[-3:]
    if (suffix == 'klb'):
        img = pyklb.readfull(image_fullname)
    elif (suffix == '.h5'):
        f = h5py.File(image_fullname, 'r+')
        # Print all root level object names (aka keys)
        # these can be group or dataset names
        #print("Keys: %s" % f.keys())
        # get first object name/key; may or may NOT be a group
        a_group_key = list(f.keys())[0]
        img = f[a_group_key][()]  # returns as a numpy array
    else :
        img = tiff.imread(image_fullname)
    #img = img[:64,200:1000,700:1700]
    d,h,w = img.shape
    print('img shape ',img.shape)
    if bLabels:
        # read label image - 100_masks_mi_0001.tif
        label_fullname = label_path % iframe
        # if doesn't exist - try
        #lux_suffix_for_embryo_alternative = '.lux_SegmentationCorrected.tif';
        #suffix_for_embryo_alternative = '.SegmentationCorrected.tif';
        alt_suffix = '.lux_SegmentationCorrected.klb'
        print(label_fullname)
        if (os.path.exists(label_fullname)):
            label_fullname = label_fullname
        else:
            label_fullname = label_fullname[:-14] + alt_suffix ## was -10
        print(label_fullname)
        suffix = label_fullname[-3:]
        if (suffix == 'klb'):
            label_img = pyklb.readfull(label_fullname,1)
            hdr = pyklb.readheader(label_fullname)
            #print(hdr)
        elif (os.path.exists(label_fullname)):
            label_img = tiff.imread(label_fullname)
        else:
            lux_part = label_fullname[-8:-4]
            #print(lux_part)
            if (lux_part == '.lux'):
                label_fullname = label_fullname[:-8] + lux_suffix_for_embryo_alternative
            else:
                label_fullname = label_fullname[:-8] + lux_suffix_for_embryo_alternative
            label_img = tiff.imread(label_fullname)
        #label_img = np.ascontiguousarray(label_img)
        #print('label img shape ',label_img.shape)
        #print('label type ',type(label_img[0,0,0]))
        #list_nuclei = np.unique(label_img)
        #print('number of labels',len(list_nuclei))
        #print(list_nuclei)
        #ind = np.where(label_img == 161)
        #print(len(ind[0]))


        #label_img = label_img[:64, 200:1000, 700:1500]

        # before downsampling  by 4 -- label image is fully upsampled
        new_label_img = InterpolateTubes(label_img, 10)
        #print('label image shape after interpolation ',new_label_img.shape)
        #print('number of labels ', len(np.unique(new_label_img)))

    print('max, min before shift ', np.amax(img), np.amin(img))
    img = img >> 4
    # added 4/2023
    print('max, min before capping ',np.amax(img),np.amin(img))
    ind = np.where(img > 255)
    img[ind] = 255
    print('max, min after capping ', np.amax(img), np.amin(img))
    img = img.astype(dtype=np.uint8)
    print('max, min, mean after conversion ',np.amax(img),np.amin(img),np.mean(img))
    #exit()
    # down sample x/y by 4  -> each pixel is 0.832 apart
    # up sample z by 2.0/(0.832) -> each pixel is 0.832 apart
    # for nuclei
    dest_size = ( 2.0/0.832, .25, .25 )
    # for membrane (already (1.0 x 1.0 x 1.0)
    #dest_size = (1.0 / 0.832, 1.0 / 0.832, 1.0 / 0.832)
    new_img = zoom(img, dest_size, order=3) # up to 5 - higher order spline interpolation
    #print('image shape after zoom ',new_img.shape)

    if bLabels:
        new_label_img = np.swapaxes(new_label_img,0,2)    

    new_img = np.swapaxes(new_img, 0, 2)
    #print('image shape after swap ',new_img.shape)
    w,h,d = new_img.shape
 
    if (iframe == (end_frame-1)):
        # no change to original
        Rev_Rotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        Rev_Offset = [0,0,0]
        reg_img = new_img
        if bLabels:
            reg_label_img = new_label_img
    else:
        # compute cumulative registration transform
        Cumulative_Rotation = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        RevR1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        RevT1 = [0, 0, 0]
        for iframe_seq in range(iframe, (end_frame-2)):
            # get current mean
            C1 = transforms[iframe_seq - start_frame]['Centroids1']  # frames are zero based start at start_frame
            meanX = C1[0]
            meanY = C1[1]
            meanZ = C1[2]
            C1 = transforms[iframe_seq - start_frame+1]['Centroids1']  # frames are off by one (zero-based)
            meanX1 = C1[0]
            meanY1 = C1[1]
            meanZ1 = C1[2]
            #print('iframe_setq C1 ', iframe_seq, meanX, meanY, meanZ)
            #print('iframe_setq C2 ', iframe_seq,meanX1,meanY1,meanZ1)
            Mean_Offset = np.array([meanX1, meanY1, meanZ1])
            Rotation =  transforms[iframe_seq - start_frame]['Rotation']
            Translation = -np.array(transforms[iframe_seq - start_frame]['Translation']) ## was [0]  # this is 1x3
            #print('translation ',Translation)
            Cumulative_Offset = np.matmul(Rotation, Mean_Offset + Translation)
            # needs to be reverse transform
            RevR2 = Rotation
            RevT2 = -Cumulative_Offset + [meanX, meanY, meanZ] #[meanX, meanY, meanZ] + DiffC

            # combine this with previous
            Rev_Rotation = np.matmul(RevR1, RevR2)
            Rev_Offset = np.matmul(RevR1, RevT2) + RevT1
            # update for next iteration
            RevR1 = Rev_Rotation
            RevT1 = Rev_Offset

        #print('Final Offset ',Rev_Offset)
        #print('Final Roation ',Rev_Rotation)
        reg_img = affine_transform(new_img,Rev_Rotation,Rev_Offset,order=3) # higher order for higher order spline interpolation
        if bLabels:
            sf = 4 # label image is 4x bigger than raw image
            lab_Offset = [value * sf for value in Rev_Offset]
            #reg_label_img = affine_transform(new_label_img, Rev_Rotation, lab_Offset, order = 0, mode = 'nearest')

    # if saving the image here - reverse x/z
    reg_img_swap = np.swapaxes(reg_img,0,2)
    if bLabels:
        dest_scale = [0.25, 0.25, 0.25]
        reg_label = zoom(new_label_img, dest_scale, order=0, mode='nearest')
        reg_label = affine_transform(reg_label, Rev_Rotation, Rev_Offset, order=0, mode='nearest')
        reg_label_swap = np.swapaxes(reg_label, 0, 2)
        # make it binary
        #ind = np.where(reg_label_swap > 0)
        #reg_label_swap[ind] = 255
        # make MIP frame
        #max_proj_img = np.max(reg_label_swap, axis=1)
        #cv2.imwrite(os.path.join(MIP_path, 'MIP_LABEL_Y_' + frame_str + '.jpg'), max_proj_img)

    #reg_img_swap = np.ascontiguousarray(reg_img_swap)
    #pyklb.writefull(reg_img_swap,os.path.join(new_image_path,'image_reg_' + frame_str + '.klb'))
    tiff.imwrite(os.path.join(new_image_path,img_type + '_reg8_' + frame_str + '.tif'),reg_img_swap)
    if bLabels:
        #reg_label_swap = np.ascontiguousarray(reg_label_swap)
        #pyklb.writefull(reg_label_swap, os.path.join(new_label_path, 'label_reg_' + frame_str + '.klb'))
        tiff.imwrite(os.path.join(new_label_path, 'label_reg8_' + frame_str + '.tif'), reg_label_swap)

    # now make max intensity projection image
    max_proj_img = np.max(reg_img, axis=1)
    cv2.imwrite(os.path.join(MIP_path,'MIP_Y8' + img_type + '_' + frame_str + '.jpg'), max_proj_img)



