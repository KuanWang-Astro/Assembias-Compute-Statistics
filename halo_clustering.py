import argparse

parser = argparse.ArgumentParser(description='Compute wp(rp) for halos')
parser.add_argument('--td',type=int,default=0,dest='td',help='time delay in seconds') 
parser.add_argument('--Lbox',type=int,required=True,dest='Lbox')
parser.add_argument('--simname',required=True,dest='simname')
parser.add_argument('--version',default='halotools_v0p4',required=True,dest='version')
parser.add_argument('--redshift',type=float,default=0.,dest='redshift')
parser.add_argument('--halofinder',default='rockstar',dest='halofinder')
parser.add_argument('--outfile',required=True,dest='outfile')
parser.add_argument('--Vmax',type=float,default=0.,dest='Vmax')
args = parser.parse_args()

import time
time.sleep(args.td)

import collections
import gc
import numpy as np
from concurrent.futures import ProcessPoolExecutor as Pool

from datetime import datetime

from halotools.sim_manager import CachedHaloCatalog

from HOD_models import decorated_hod_model
from HOD_models import standard_hod_model

from halotools.empirical_models import MockFactory

from halotools.mock_observables import return_xyz_formatted_array
from halotools.mock_observables import counts_in_cylinders
from halotools.mock_observables import void_prob_func
from halotools.mock_observables import wp
from halotools.utils import randomly_downsample_data
from halotools.utils import crossmatch

##########################################################
output_names = ('wprp','null')


Lbox = args.Lbox

pi_max = 60
r_wp = np.logspace(-1, np.log10(Lbox)-1, 20)
##wp

##########################################################

def calc_all_observables():
    
    gc.collect()
    
    output = []
        
    if args.Vmax!=0:
        mask_Vmax = halocat.halo_table['halo_vmax']>args.Vmax
        
    pos_halos = return_xyz_formatted_array(halocat.halo_table['halo_x'][mask_Vmax],\
                                           halocat.halo_table['halo_y'][mask_Vmax],\
                                           halocat.halo_table['halo_z'][mask_Vmax],period=Lbox)
    
    # wprp
    output.append(wp(pos_halos, r_wp, pi_max, period=Lbox))
    
    output.append(np.arange(1))

    
    return output


############################################################
consuelo20_box_list = ['0_4001','0_4002','0_4003','0_4004','0_4020','0_4026','0_4027','0_4028','0_4029','0_4030',\
            '0_4032','0_4033','0_4034','0_4035','0_4036','0_4037','0_4038','0_4039','0_4040']


def main(output_fname):
    
    output_dict = collections.defaultdict(list)
    
    global halocat
    
    if args.simname=='consuelo20' and args.version=='all':
        for box in consuelo20_box_list:
            halocat = CachedHaloCatalog(simname = args.simname, version_name = box,redshift = args.redshift, \
                                        halo_finder = args.halofinder)
            output_data = calc_all_observables()
            for name, data in zip(output_names, output_data):
                output_dict[name].append(data)
            print box
    else:
        halocat = CachedHaloCatalog(simname = args.simname, version_name = args.version,redshift = args.redshift, \
                                            halo_finder = args.halofinder)
        output_data = calc_all_observables()
        for name, data in zip(output_names, output_data):
            output_dict[name].append(data)

    for name in output_names:
        output_dict[name] = np.array(output_dict[name])

    np.savez(output_fname, **output_dict)


if __name__ == '__main__':
    main(args.outfile)
    f = open(args.outfile+'_log','w')
    for arg in vars(args):
        f.write(str(arg)+':'+str(getattr(args, arg))+'\n')

    f.close()


