import numpy as np
from halotools.sim_manager import UserSuppliedHaloCatalog as USHC

redshift = 0.0
Lbox = 420
particle_mass = 0.187e10

hlist4038 = np.loadtxt('/home/yymao/consuelo/4038/hlist_1.00000.list')

print 'loaded'

hlist = hlist4038
h_scale_factor = hlist[:,0]
h_id = hlist[:,1].astype(np.int64)
h_pid = hlist[:,5]
h_upid = hlist[:,6]
h_mvir = hlist[:,9]
h_rvir = hlist[:,11]/1000.
h_rs = hlist[:,12]/1000.
h_scale_factor_last_mm = hlist[:,15]
h_vmax = hlist[:,16]
h_x = hlist[:,17]
h_y = hlist[:,18]
h_z = hlist[:,19]
h_vx = hlist[:,20]
h_vy = hlist[:,21]
h_vz = hlist[:,22]
h_jx = hlist[:,27]
h_jy = hlist[:,28]
h_jz = hlist[:,29]
h_spin = hlist[:,30]
h_macc = hlist[:,23]
h_mpeak = hlist[:,24]
h_vacc = hlist[:,25]
h_vpeak = hlist[:,26]

h_vrms = hlist[:,13]

h_nfw_conc = h_rvir/h_rs
print 'columns'


halo_catalog = USHC(redshift = redshift, Lbox = Lbox, particle_mass = particle_mass, number_of_particles = 1400**3,\
                    omega_m = 0.25, omega_L = 0.75, omega_b = 0.04, h = 0.7, sigma8 = 0.8, ns = 1.0,\
                    halo_scale_factor = h_scale_factor, halo_id = h_id, halo_pid = h_pid, halo_upid = h_upid,\
                    halo_mvir = h_mvir, halo_rvir = h_rvir, halo_rs = h_rs,\
                    halo_scale_factor_last_mm = h_scale_factor_last_mm, halo_vmax = h_vmax,\
                    halo_x = h_x, halo_y = h_y, halo_z = h_z, halo_vx = h_vx, halo_vy = h_vy, halo_vz = h_vz,\
                    halo_jx = h_jx, halo_jy = h_jy, halo_jz = h_jz, halo_spin = h_spin,\
                    halo_macc = h_macc, halo_mpeak = h_mpeak, halo_vacc = h_vacc, halo_vpeak = h_vpeak,\
                    halo_vrms = h_vrms, halo_nfw_conc = h_nfw_conc)
print 'catalog'

halo_catalog.add_halocat_to_cache('/home/kuw8/.astropy/cache/halotools/halo_catalogs/consuelo20/rockstar/hlist_1.00000.list.0_4038.hdf5', 'consuelo20', 'rockstar', '0_4038', '#')
print 'cached'

