import matplotlib
#matplotlib.use('Agg')
from pylab import *
import os

#data_dir = '/aurora_nobackup/eva/romerowo/nutau_outputs/20161115/1.5km_ice'
#tag = '1.5km_ice'
#data_dir = '/aurora_nobackup/eva/romerowo/nutau_outputs/20161115/3.7km_ocn'
#tag = '3.7km_ocn'
#data_dir = '/aurora_nobackup/eva/romerowo/nutau_outputs/20161128/1.5km_ice_lowCS'
#tag = '1.5km_ice_lowCS'
#data_dir = '/aurora_nobackup/eva/romerowo/nutau_outputs/20161128/3.7km_ocn_lowCS'
#tag = '3.7km_ocn_lowCS'

tag = '0.0km_ice_midCS_stdEL'
#tag = '0.0km_ice_regCS_lowEL'
#tag = '0.0km_ice_lowCS_regEL'
#tag = '0.0km_ice_lowCS_lowEL'
#tag = '1.0km_ice_regCS_regEL'
#tag = '1.0km_ice_regCS_lowEL'
#tag = '1.0km_ice_lowCS_regEL'
#tag = '1.0km_ice_lowCS_lowEL'
#tag = '2.0km_ice_regCS_regEL'
#tag = '2.0km_ice_regCS_lowEL'
#tag = '2.0km_ice_lowCS_regEL'
#tag = '2.0km_ice_lowCS_lowEL'
#tag = '3.0km_ice_regCS_regEL'
#tag = '3.0km_ice_regCS_lowEL'
#tag = '3.0km_ice_lowCS_regEL'
#tag = '3.0km_ice_lowCS_lowEL'
#tag = '4.0km_ice_regCS_regEL'
#tag = '4.0km_ice_regCS_lowEL'
#tag = '4.0km_ice_lowCS_regEL'
#tag = '4.0km_ice_lowCS_lowEL'

def read_emerging(filename):
    lc = 0
    num_CC = []
    num_NC = []
    num_decays = []
    num_particles = []
    energy = []
    for line in file(filename):
        if(lc!=0 and 'END' not in line):
            num_CC.append(int(line.split()[0]))
            num_NC.append(int(line.split()[1]))
            num_decays.append(int(line.split()[2]))
            num_particles.append(int(line.split()[3]))
            energy.append(float(line.split()[4]))
            #print line
        lc+=1
    return np.array(num_CC), np.array(num_NC), np.array(num_decays), np.array(num_particles), np.array(energy)

tag_list = []
for thickness in ['0.0', '1.0', '2.0', '3.0', '4.0']:
  for CS in ['low','mid', 'upp']:
    #for EL in ['std', 'low']:
    for EL in ['std', 'low']:
        tag_list.append('%skm_ice_%sCS_%sEL'%(thickness, CS, EL))

print tag_list
for tag in tag_list:
    data_dir = '/aurora_nobackup/eva/romerowo/nutau_outputs/20170409/%s'%tag
    count = 0
    count_true=0
    ang_array = np.arange(85.0, 90.1, 0.1)
    ang_array = np.concatenate([ np.arange(90.,95.,0.1) , np.arange(95.,180.,1.) ])
    th_exit_array = 90.-ang_array
    e_array = np.array([1e+15, 3e+15, 1e+16, 3e+16, 1e+17, 3e+17, 1e+18, 3e+18, 1e+19, 3e+19, 1e+20, 3e+20, 1e+21])

    missing_count = 0
    outdir = './LUTs/%s'%tag
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    for e in e_array:
        e_s = '%1.0e'%e
        #e_s = e_s.replace('+','')
        #print e_s
        data_array  = []
        for ang in ang_array:
            fnm = data_dir + '/Emerging_tau_info_reg_%s_%1.1f_%s.dat'%(e_s, ang, tag)
            print fnm
            if not os.path.exists(fnm):
                print fnm, os.path.exists(fnm)
                missing_count += 1
                data_array.append([])
            if os.path.exists(fnm):
                num_CC, num_NC, num_decays, num_particles, energy = read_emerging(fnm)
                print fnm, len(energy)
                data_array.append(energy)
        np.savez('%s/LUT_%s_eV.npz'%(outdir,e_s), data_array = data_array, th_exit_array = th_exit_array)
    print 'missing_count', missing_count

    '''
    fnm = 'Emerging_tau_info_reg_1e19_94.1.dat'
    print 'reading file'
    num_CC, num_NC, num_decays, num_particles, energy = read_emerging(fnm)
    '''

