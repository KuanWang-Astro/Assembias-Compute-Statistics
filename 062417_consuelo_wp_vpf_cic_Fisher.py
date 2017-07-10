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

##########################################################

param_names = ('alpha','logM1','sigma_logM','logM0','logMmin','mean_occupation_centrals_assembias_param1','mean_occupation_satellites_assembias_param1')
output_names = ('ngals','Pcic','vpf','wprp','param')

##########################################################

Lbox = 420

proj_search_radius = 2.0         ##a cylinder of radius 2 Mpc/h
cylinder_half_length = 10.0      ##half-length 10 Mpc/h

##cic

r_vpf = np.logspace(0, 1.3, 20)
num_spheres = int(1e5)
##vpf

pi_max = 60
r_wp = np.logspace(-1, np.log10(Lbox)-1, 20)
##wp

##########################################################

def calc_all_observables(param):

    model.param_dict.update(dict(zip(param_names, param)))    ##update model.param_dict with pairs (param_names:params)

    try:
        model.mock.populate()
    except:
        model.populate_mock(halocat)
    
    gc.collect()
    
    output = []


    pos_gals_d = return_xyz_formatted_array(*(model.mock.galaxy_table[ax] for ax in 'xyz'), \
            velocity=model.mock.galaxy_table['vz'], velocity_distortion_dimension='z',\
                                          period=model.mock.Lbox)             ##redshift space distorted
    pos_gals_d = np.array(pos_gals_d,dtype=float)
    
    # ngals
    output.append(model.mock.galaxy_table['x'].size)
    
    # Pcic
    output.append(np.bincount(counts_in_cylinders(pos_gals_d, pos_gals_d, proj_search_radius, \
            cylinder_half_length), minlength=100)[1:100]/float(model.mock.galaxy_table['x'].size))

    
    # vpf
    output.append(void_prob_func(pos_gals_d, r_vpf, num_spheres, period=model.mock.Lbox))
    
    # wprp
    output.append(wp(pos_gals_d, r_wp, pi_max, period=model.mock.Lbox))
    
    # parameter set
    output.append(param)
    
    return output


############################################################
box_list = ['0_4001','0_4002','0_4003','0_4004','0_4020','0_4026','0_4027','0_4028','0_4029','0_4030',\
            '0_4032','0_4033','0_4034','0_4035','0_4036','0_4037','0_4038','0_4039','0_4040']


def main(model_gen_func, params_fname, params_usecols, output_fname):
    global model
    model = model_gen_func()

    nparams = 1540
    params = np.loadtxt(params_fname, usecols=params_usecols)
    params_median = np.median(params, axis=0)
    params = params_median*np.ones((nparams,len(params_usecols)))
    nreal = 20
    for i in range(7):
        for j in range(11):
            for k in range(nreal):
                params[i*11*nreal+j*nreal+k,i] = params_median[i]-0.05+0.01*j
                
    output_dict = collections.defaultdict(list)
    nproc = 55
    
    global halocat
    
    with Pool(nproc) as pool:
        for box in box_list:
            halocat = CachedHaloCatalog(simname = 'consuelo20', version_name = box,redshift = 0, \
                            halo_finder = 'rockstar')
            model.populate_mock(halocat)
            for i, output_data in enumerate(pool.map(calc_all_observables, params)):
                if i%55 == 54:
                    print i
                    print str(datetime.now())
                for name, data in zip(output_names, output_data):
                    output_dict[name].append(data)
            print box
    
    for name in output_names:
        output_dict[name] = np.array(output_dict[name])

    np.savez(output_fname, **output_dict)


if __name__ == '__main__':
    main(decorated_hod_model, '../ABMCMCfiles/corr1_wp20.0.abfit.covar.chain', range(7), '062417_consuelo_Fisher_w')
    print 'w_1_20 done'


