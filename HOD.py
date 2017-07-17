from concurrent.futures import ProcessPoolExecutor as Pool

import numpy as np
import halotools
from halotools.sim_manager import CachedHaloCatalog
import math
import random
from datetime import datetime

########

halocat = CachedHaloCatalog(simname = 'bolplanck', version_name = 'halotools_v0p4',redshift = 0, \
                            halo_finder = 'rockstar')

N_halo = halocat.halo_table['halo_mvir'].size
conc_split = np.median(halocat.halo_table['halo_nfw_conc'])

def poisson_distribution(k,l):
    return np.power(l,k)*np.exp(-l)/math.factorial(k)

########

def standard_hod_cen_moment1(Mvir,logMmin,sigma_logM):
    return 0.5*(1+math.erf((np.log10(Mvir)-logMmin)/sigma_logM))

def standard_hod_sat_moment1(Mvir,logM0,logM1,alpha):
    if Mvir>np.power(10,logM0):
        return np.power((Mvir-np.power(10,logM0))/np.power(10,logM1),alpha)
    else:
        return 0.0

def concentration_split(conc):
    if conc<conc_split:
        return -1
    else:
        return 1
    
def decorated_hod_cen_moment1(Mvir,logMmin,sigma_logM,Acen,conc):
    standard = 0.5*(1+math.erf((np.log10(Mvir)-logMmin)/sigma_logM))
    if standard>0.5:
        return standard+concentration_split(conc)*Acen*(1.0-standard)
    else:
        return standard*(1.0+concentration_split(conc)*Acen)
    
def decorated_hod_sat_moment1(Mvir,logM0,logM1,alpha,Asat,conc):
    if Mvir>np.power(10,logM0):
        return (1.0+Asat*concentration_split(conc))*np.power((Mvir-np.power(10,logM0))/np.power(10,logM1),alpha)
    else:
        return 0.0

###########


def calculate_HOD(pset):
    P_n = np.zeros(10)
        
    for n in range(9):
        prob = 0
        for q in range(N_halo):
            Mvir = halocat.halo_table['halo_mvir'][q]
            conc = halocat.halo_table['halo_nfw_conc'][q]
            if pset.size==5:
                prob += standard_hod_cen_moment1(Mvir,pset[4],pset[2])*\
                    poisson_distribution(n,standard_hod_sat_moment1(Mvir,pset[3],pset[1],pset[0]))
            else:
                prob += decorated_hod_cen_moment1(Mvir,pset[4],pset[2],pset[5],conc)*\
                    poisson_distribution(n,decorated_hod_sat_moment1(Mvir,pset[3],pset[1],pset[0],pset[6],conc))
        print n
        P_n[n+1] = prob/N_halo
    prob0 = 0
    for q in range(N_halo):
        Mvir = halocat.halo_table['halo_mvir'][q]
        conc = halocat.halo_table['halo_nfw_conc'][q]
        if pset.size==5:
            prob0 += 1-standard_hod_cen_moment1(Mvir,pset[4],pset[2])
        else:
            prob0 += 1-decorated_hod_cen_moment1(Mvir,pset[4],pset[2],pset[5],conc)
    P_n[0] = prob0/N_halo

    
    return P_n



def main(AB_bool,input_filename,output_filename,nparams):
    f1 = open(output_filename,'w')
    usecols = range(5+2*AB_bool)
    params = np.loadtxt(input_filename, usecols=usecols)
    params = params[np.random.choice(len(params), nparams)]
    nproc = 55
    with Pool(nproc) as pool:
        for i, ret in enumerate(pool.map(calculate_HOD, params)):
            ret.tofile(f1)
            if i%55 == 54:
 	        print i, 'of', nparams
 	        print str(datetime.now())
    f1.close()

if __name__ == '__main__':
#    main(False,'../ABMCMCfiles/corr1_wp20.0.covar.chain','HOD_1_20_wo',1000)
    main(True,'../ABMCMCfiles/corr1_wp20.0.abfit.covar.chain','HOD_1_20_w',1000)
