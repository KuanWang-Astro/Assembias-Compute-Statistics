import numpy as np
from halotools.sim_manager import UserSuppliedHaloCatalog as USHC

redshift = 0.0
Lbox = 500
particle_mass = 1.03395e10

diemer500 = np.loadtxt('/home/kuw8/Assembly Bias Project/diemer_antonio/z0.0_L0500.catalog')

print 'loaded'

hlist = diemer500

h_id = hlist[:,0].astype(np.int64)
h_descid = hlist[:,1].astype(np.int64)
h_mvir = hlist[:,2]
h_vmax = hlist[:,3]
h_vrms = hlist[:,4]
h_rvir = hlist[:,5]/1000.
h_rs = hlist[:,6]/1000.
h_x = hlist[:,8]
h_y = hlist[:,9]
h_z = hlist[:,10]
h_vx = hlist[:,11]
h_vy = hlist[:,12]
h_vz = hlist[:,13]
h_jx = hlist[:,14]
h_jy = hlist[:,15]
h_jz = hlist[:,16]
h_spin = hlist[:,17]
h_rs_klypin = hlist[:,18]
h_m200b_all = hlist[:,19]
h_m200b = hlist[:,20]
h_m200c = hlist[:,21]
h_m500c = hlist[:,22]
h_m2500c = hlist[:,23]
h_xoff = hlist[:,24]
h_voff =hlist[:,25]
h_spin_bullock = hlist[:,26]
h_b_to_a = hlist[:,27]
h_c_to_a = hlist[:,28]
h_Ax = hlist[:,29]
h_Ay = hlist[:,30]
h_Az = hlist[:,31]
h_b_to_a_500c = hlist[:,32]
h_c_to_a_500c = hlist[:,33]
h_Ax_500c = hlist[:,34]
h_Ay_500c = hlist[:,35]
h_Az_500c = hlist[:,36]
h_t_by_u = hlist[:,37]
h_halfmass_radius = hlist[:40]/1000.
h_pid = hlist[:,41].astype(np.int64)
h_nfw_conc = h_rvir/h_rs

h_upid = np.copy(h_pid)
for i in range(len(h_pid)):
    j = i
    while(h_pid[j]!=-1):
        j = np.where(h_id==h_pid[j])
        h_upid[i] = h_id[j]

        
print 'columns'


halo_catalog = USHC(redshift = redshift, Lbox = Lbox, particle_mass = particle_mass,\
                    omega_m = 0.32, omega_L = 0.68, h = 0.67,\
                    halo_id = h_id, halo_descid = h_descid, halo_mvir = h_mvir, halo_vmax = h_vmax,\
                    halo_vrms = h_vrms, halo_rvir = h_rvir, halo_rs = h_rs,\
                    halo_x = h_x, halo_y = h_y, halo_z = h_z, halo_vx = h_vx, halo_vy = h_vy, halo_vz = h_vz,\
                    halo_jx = h_jx, halo_jy = h_jy, halo_jz = h_jz, halo_spin = h_spin,\
                    halo_rs_klypin = h_rs_klypin, halo_m200b_all = h_m200b_all,\
                    halo_m200b = h_m200b, halo_m200c = h_m200c, halo_m500c = h_m500c,\
                    halo_m2500c = h_m2500c, halo_xoff = h_xoff, halo_voff = h_voff,\
                    halo_spin_bullock = h_spin_bullock, halo_b_to_a = h_b_to_a,\
                    halo_c_to_a = h_c_to_a, halo_Ax = h_Ax, halo_Ay = h_Ay, halo_Az = h_Az,\
                    halo_b_to_a_500c = h_b_to_a_500c, halo_c_to_a_500c = h_c_to_a_500c,\
                    halo_Ax_500c = h_Ax_500c, halo_Ay_500c = h_Ay_500c, halo_Az_500c = h_Az_500c,\
                    halo_t_by_u = h_t_by_u, halo_halfmass_radius = h_halfmass_radius,\
                    halo_pid = h_pid, halo_nfw_conc = h_nfw_conc, halo_upid = h_upid)
print 'catalog'

halo_catalog.add_halocat_to_cache('/home/kuw8/.astropy/cache/halotools/halo_catalogs/diemerL0500/rockstar/z0.0_L0500.catalog.antonio.hdf5', 'diemerL0500', 'rockstar', 'antonio', '#')
print 'cached'

