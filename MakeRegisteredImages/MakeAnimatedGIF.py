import glob
from PIL import Image
import os
import argparse
import yaml

# make animated gif from MIP-Frames

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config_name', type=str, default='config.yaml',
                        help='the name of the config file in CurrentConfig directory')
    args = parser.parse_args()

    config_path = '/mnt/home/lbrown/LineageConstruction/CurrentConfig/';
    config_name = args.config_name
    print('Config ', config_name)
    with open(os.path.join(config_path, config_name), 'r') as file:
        config_opts = yaml.safe_load(file)

    out_path = config_opts["output_dir"]
    frame_folder = os.path.join(out_path, 'MIP-Frames')

    embryo_base = 'MIP_Y8'
    img_type = 'nuclei_'
    #img_type = 'membrane_'
    frames = []
    start = config_opts['register_begin_frame']
    end = config_opts['register_end_frame']
    print('start, end ', start, end)
    for iframe in range(start,end):
        range_str = str(start) + '_' + str(end)
        frame_str = '%05d' % iframe
        frame_str = '%d' % iframe
        embryo_file = embryo_base + img_type + frame_str + '.jpg'
        fname = os.path.join(frame_folder,embryo_file)
        print(fname)
        img = Image.open(fname)
        frames.append (img)

    frame_one = frames[0]
    frame_one.save(os.path.join(frame_folder,'RegisteredMIP' + img_type + range_str + '.gif'), format="GIF", append_images=frames,
                   save_all=True, duration=60000) # loop=1)  # ~2 sec/frame ?