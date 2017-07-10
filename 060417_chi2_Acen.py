import collections
import gc
import numpy as np
from concurrent.futures import ProcessPoolExecutor as Pool

from datetime import datetime

from halotools.sim_manager import CachedHaloCatalog
from halotools.empirical_models import PrebuiltHodModelFactory
from halotools.empirical_models import AssembiasZheng07Cens
from halotools.empirical_models import TrivialPhaseSpace
from halotools.empirical_models import AssembiasZheng07Sats
from halotools.empirical_models import NFWPhaseSpace
from halotools.empirical_models import HodModelFactory
from halotools.empirical_models import PrebuiltHodModelFactory
from halotools.empirical_models import MockFactory
from halotools.mock_observables import return_xyz_formatted_array
from halotools.mock_observables import counts_in_cylinders
from halotools.mock_observables import delta_sigma
from halotools.mock_observables import void_prob_func
from halotools.mock_observables import wp
from halotools.utils import randomly_downsample_data

halocat = CachedHaloCatalog(simname = 'bolplanck', version_name = 'halotools_v0p4',redshift = 0, \
                            halo_finder = 'rockstar')


def decorated_hod_model():
    cen_occ_model = AssembiasZheng07Cens(prim_haloprop_key='halo_mvir', sec_haloprop_key='halo_nfw_conc')
    cen_prof_model = TrivialPhaseSpace()
    sat_occ_model = AssembiasZheng07Sats(prim_haloprop_key='halo_mvir', sec_haloprop_key='halo_nfw_conc')
    sat_prof_model = NFWPhaseSpace()
    return HodModelFactory(centrals_occupation=cen_occ_model, centrals_profile=cen_prof_model, satellites_occupation=sat_occ_model, satellites_profile=sat_prof_model)


def standard_hod_model():
    return PrebuiltHodModelFactory('zheng07', threshold=-20)


param_names = ('alpha','logM1','sigma_logM','logM0','logMmin','mean_occupation_centrals_assembias_param1','mean_occupation_satellites_assembias_param1')
output_names = ('ngals','Pcic','ggl','vpf','wprp','Acen')

##########################################################

proj_search_radius = 2.0         ##a cylinder of radius 2 Mpc/h
cylinder_half_length = 10.0      ##half-length 10 Mpc/h

##cic

rp_bins_ggl = np.logspace(-1, 1.9, 30)
num_ptcls_to_use = int(1e5)
particle_masses = np.zeros(num_ptcls_to_use)+halocat.particle_mass
total_num_ptcls_in_snapshot = halocat.num_ptcl_per_dim**3
downsampling_factor = total_num_ptcls_in_snapshot/float(num_ptcls_to_use)
##ggl

rbins = np.logspace(0, 1.9, 20)
num_spheres = int(1e5)
##vpf

pi_max = 60
rp_bins = np.logspace(-1, 1.9, 30)
rp_bin_centers = (rp_bins[:1] + rp_bins[1:])/2.
##wp

##########################################################

pos_part = return_xyz_formatted_array(*(halocat.ptcl_table[ax] for ax in 'xyz'), period=halocat.Lbox)
    
pos_part = randomly_downsample_data(pos_part, num_ptcls_to_use)

def calc_all_observables(param):

    model.param_dict.update(dict(zip(param_names, param)))    ##update model.param_dict with pairs (param_names:params)

    try:
        model.mock.populate()
    except:
        model.populate_mock(halocat)
    
    gc.collect()
    
    output = []

    pos_gals_u = return_xyz_formatted_array(*(model.mock.galaxy_table[ax] for ax in 'xyz'), 
                                          period=model.mock.Lbox)             ##redshift space undistorted
    pos_gals_u = np.array(pos_gals_u,dtype=float)
    pos_gals_d = return_xyz_formatted_array(*(model.mock.galaxy_table[ax] for ax in 'xyz'), \
            velocity=model.mock.galaxy_table['vz'], velocity_distortion_dimension='z',\
                                          period=model.mock.Lbox)             ##redshift space distorted
    pos_gals_d = np.array(pos_gals_d,dtype=float)
    
    # ngals
    output.append(model.mock.galaxy_table['x'].size)
    
    # Pcic
    output.append(np.bincount(counts_in_cylinders(pos_gals_d, pos_gals_d, proj_search_radius, \
                                                  cylinder_half_length), minlength=100)[1:100]/float(model.mock.galaxy_table['x'].size))

    # delta sigma
    output.append(delta_sigma(pos_gals_u, pos_part, particle_masses=particle_masses, downsampling_factor=downsampling_factor, rp_bins=rp_bins_ggl, period=model.mock.Lbox)[1])

    
    # vpf
    output.append(void_prob_func(pos_gals_d, rbins, num_spheres, period=model.mock.Lbox))
    
    # wprp
    output.append(wp(pos_gals_d, rp_bins, pi_max, period=model.mock.Lbox))
    
    # Acen
    output.append(param[5])
    
    return output


## median of corr1
median_wo = np.array((1.14385007593,13.2858403826,0.348464903173,11.3075027005,11.9718570282))
median_w = np.array((1.02654941214,13.1879106954,0.8781897069,12.1031391855,12.2692942798,0.91596941012,0.0258608345476))

def main(model_gen_func, Acen_zp, params_usecols, output_fname):
    global model
    model = model_gen_func()

    nparams = 4200

    params = median_w*np.ones((nparams,7))  ##take medians for other parameters than Acen and Asat
    params[:,6] = np.zeros(nparams)    ##set Asat=0

    for i in range(21):
        params[i*200:i*200+200,5] = Acen_zp+0.01*i-0.1


    output_dict = collections.defaultdict(list)
    nproc = 55
    with Pool(nproc) as pool:
        for i, output_data in enumerate(pool.map(calc_all_observables, params)):
            if i%55 == 54:
                print i, 'of', nparams
                print str(datetime.now())
            for name, data in zip(output_names, output_data):
                output_dict[name].append(data)
    
    for name in output_names:
	print name
        output_dict[name] = np.array(output_dict[name])

    np.savez(output_fname, **output_dict)


if __name__ == '__main__':
    main(decorated_hod_model, 0.8, range(7), '060417_p08')
    print '0.8 done'
    main(decorated_hod_model, 0.0, range(7), '060417_00')
    print '0.0 done'
    main(decorated_hod_model, -0.8, range(7), '060417_m08')
    print '-0.8 done'


