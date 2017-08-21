import numpy as np
from halotools.sim_manager import UserSuppliedHaloCatalog as USHC

redshift = 0.0
Lbox = 250
particle_mass = 1.44e8

hlist = np.loadtxt('/home/yymao/chinchilla/c250-2048/hlist_1.00000.list')

print 'loaded'

h_scale_factor = hlist[:,0]
h_id = hlist[:,1].astype(np.int64)
h_pid = hlist[:,5]
h_upid = hlist[:,6]
h_mvir = hlist[:,10]
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
h_jx = hlist[:,23]
h_jy = hlist[:,24]
h_jz = hlist[:,25]
h_spin = hlist[:,26]
h_rs_klypin = hlist[:,34]
h_mvir_all = hlist[:,35]
h_m200b = hlist[:,36]
h_m200c = hlist[:,37]
h_m500c = hlist[:,38]
h_m2500c = hlist[:,39]
h_xoff = hlist[:,40]
h_voff = hlist[:,41]
h_spin_bullock = hlist[:,42]
h_b_to_a = hlist[:,43]
h_c_to_a = hlist[:,44]
h_Ax = hlist[:,45]
h_Ay = hlist[:,46]
h_Az = hlist[:,47]
h_b_to_a_500c = hlist[:,48]
h_c_to_a_500c = hlist[:,49]
h_Ax_500c = hlist[:,50]
h_Ay_500c = hlist[:,51]
h_Az_500c = hlist[:,52]
h_t_by_u = hlist[:,53]
h_macc = hlist[:,54]
h_mpeak = hlist[:,55]
h_vacc = hlist[:,56]
h_vpeak = hlist[:,57]
h_halfmass_scale_factor = hlist[:,58]
h_dmvir_dt_inst = hlist[:,59]
h_dmvir_dt_100myr = hlist[:,60]
h_dmvir_dt_tdyn = hlist[:,61]

h_vrms = hlist[:,13]

h_nfw_conc = h_rvir/h_rs
print 'columns'


halo_catalog = USHC(redshift = redshift, Lbox = Lbox, particle_mass = particle_mass, number_of_particles = 2048**3,\
                    omega_m = 0.286, omega_L = 0.714, omega_b = 0., h = 0.7,\
                    halo_scale_factor = h_scale_factor, halo_id = h_id, halo_pid = h_pid, halo_upid = h_upid,\
                    halo_mvir = h_mvir, halo_rvir = h_rvir, halo_rs = h_rs,\
                    halo_scale_factor_last_mm = h_scale_factor_last_mm, halo_vmax = h_vmax,\
                    halo_x = h_x, halo_y = h_y, halo_z = h_z, halo_vx = h_vx, halo_vy = h_vy, halo_vz = h_vz,\
                    halo_jx = h_jx, halo_jy = h_jy, halo_jz = h_jz, halo_spin = h_spin,\
                    halo_rs_klypin = h_rs_klypin,\
                    halo_mvir_all = h_mvir_all, halo_m200b = h_m200b, halo_m200c = h_m200c,\
                    halo_m500c = h_m500c, halo_m2500c = h_m2500c, halo_xoff = h_xoff, halo_voff = h_voff,\
                    halo_spin_bullock = h_spin_bullock, halo_b_to_a = h_b_to_a, halo_c_to_a = h_c_to_a,\
                    halo_Ax = h_Ax, halo_Ay = h_Ay, halo_Az = h_Az, halo_b_to_a_500c = h_b_to_a_500c,\
                    halo_c_to_a_500c = h_c_to_a_500c, halo_Ax_500c = h_Ax_500c,\
                    halo_Ay_500c = h_Ay_500c, halo_Az_500c = h_Az_500c, halo_t_by_u = h_t_by_u,\
                    halo_macc = h_macc, halo_mpeak = h_mpeak, halo_vacc = h_vacc, halo_vpeak = h_vpeak,\
                    halo_halfmass_scale_factor = h_halfmass_scale_factor,\
                    halo_dmvir_dt_inst = h_dmvir_dt_inst, halo_dmvir_dt_100myr = h_dmvir_dt_100myr,\
                    halo_dmvir_dt_tdyn = h_dmvir_dt_tdyn, halo_vrms = h_vrms, halo_nfw_conc = h_nfw_conc)
print 'catalog'

halo_catalog.add_halocat_to_cache('/home/kuw8/.astropy/cache/halotools/halo_catalogs/chinchilla/rockstar/hlist_1.00000.list.250-2048.hdf5', 'chinchilla', 'rockstar', '250-2048', '#')
print 'cached'

