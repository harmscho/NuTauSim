from pylab import *
from matplotlib.colors import LogNorm
from matplotlib import gridspec
rcParams['figure.facecolor'] = 'white'
rcParams['font.size']=18
def load_LUT(LUT_fnm):
    f = np.load(LUT_fnm)

    th_array = f['th_exit_array']
    th_em_array = -th_array
    data_array = f['data_array']
    P_exit = np.zeros(len(th_array))
    mean_num_CC     = f['mean_num_CC']
    mean_num_NC     = f['mean_num_NC']
    mean_num_decays = f['mean_num_decays']
    rms_num_CC      = f['rms_num_CC']
    rms_num_NC      = f['rms_num_NC']
    rms_num_decays  = f['rms_num_decays']
    for k in range(0,len(th_array)):
        P_exit[k] = float(len(data_array[k]))/1.e6
    return th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays

def get_mean_sigma(sorted_array, conf=0.68):
    N = len(sorted_array)
    k_conf = int(float(N) *(conf))
    k_init = range(0,N-k_conf)
    k_stop = range(k_conf, len(sorted_array))
    #print N, k_conf, len(k_init), len(k_stop)
    intervals = sorted_array[k_stop]-sorted_array[k_init]
    k1_conf = np.argmin(intervals)
    return sorted_array[k_init[k1_conf]], sorted_array[k_stop[k1_conf]]
    #return 0.5*(sorted_array[k1_68]+sorted_array[k1_68+k_stop[0]]), 0.5*(sorted_array[k1_68+k_stop[0]]- sorted_array[k1_68])


#plot_type = 'energy'
#plot_type = 'thickness'
plot_type = 'models'
#plot_type = 'regen'

energy_list = ['1e+16', '3e+16', '1e+17', '3e+17', '1e+18', '3e+18', '1e+19', '3e+19','1e+20', '3e+20', '1e+21']
#energy_list = ['1e+15', '3e+15', '1e+16', '3e+16', '1e+17', '3e+17', '1e+18', '3e+18', '1e+19', '3e+19','1e+20', '3e+20', '1e+21']
ice_thickness_list = ['0.0', '1.0', '2.0', '3.0', '4.0']

if(plot_type =='energy'):
    colors = cm.hot(np.linspace(0, 1, int(2.*float(len(energy_list)))))
    #figure(1, figsize=(8,10))
    figure(1, figsize=(8,6))
    figure(2, figsize=(8,10))
    figure(3, figsize=(8,10))
    figure(4, figsize=(10,10))
    gs = gridspec.GridSpec(5, 2, width_ratios=[1, 3]) 
    cc = 0
    ice_thickness = '4.0'
    ax2_array = [] 
    count2=0
    count3=0
    for energy in energy_list[::-1]:
        th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays = load_LUT('/home/romerowo/nutau_sim/LUTs/%skm_ice_midCS_stdEL/LUT_%s_eV.npz'%(ice_thickness,energy))
        #print P_exit
        figure(1)
        ax=subplot(111)
        ax.set_xscale('log')
        semilogy(th_em_array[P_exit>0.], P_exit[P_exit>0.], '-', lw=2, color=colors[cc], label='%s eV'%energy)
        new_ticks = np.array([0.1, 0.3, 1., 3., 10.,  30.])
        xticks(new_ticks, new_ticks, fontsize=16)
        yticks(fontsize=15)
        xlim(0.,50.)
        ylim(1.e-6,1.)
        legend(loc=1, fontsize=14, borderpad=0.1, borderaxespad=0, labelspacing=0.1)
        grid(True, which='both')
        xlabel('Emergence Angle, degrees')
        ylabel('Probability of Tau Exit')

        '''
        ax=subplot(211)
        semilogy(th_em_array[P_exit>0.], P_exit[P_exit>0.], '-', lw=2, color=colors[cc], label='%s eV'%energy)
        yticks(fontsize=15)
        xlim(0.,50.)
        ylim(1.e-6,1.)
        legend(loc=1, fontsize=14, borderpad=0.1, borderaxespad=0, labelspacing=0.1)
        grid(True, which='both')
        xlabel('Emergence Angle, degrees')
        ylabel('Probability of Tau Exit')

        subplot(212)
        semilogy(th_em_array[P_exit>0.], P_exit[P_exit>0.], '-', lw=2, color=colors[cc])
        xlim(0.,5.)
        ylim(1.e-5,1.)
        xlabel('Emergence Angle, degrees')
        ylabel('Probability of Tau Exit')
        grid(True, which='both')
        '''
        # 2dhistogram version
        if '3' not in energy and float(energy)>3.e16:

            #subplot(7,2,count3*2+2)
            img = np.zeros((51,len(arange(14.,21.1,0.1))-1))
            th_em_grid, E_tau_grid = np.meshgrid(th_em_array[range(0,51)], arange(14.,21.05,0.1))
            for k in range(0,51):
                #print th_em_array[k]
                hist, bin_edges = np.histogram(data_array[k], bins=arange(14.,21.1,0.1))
                img[k,:] = hist/1.e6
                #semilogy(bin_edges[1:], hist+1.e-3)
                #ylim(0.5,1.e5)
            img = np.rot90(img)
            img = np.flipud(img)
            #print th_em_grid.shape, E_tau_grid.shape, img.shape
            figure(3)
            '''
            ax3=subplot(gs[count3*2])
            ax3.set_yscale('log')
            mu = np.zeros(51)
            eL = np.zeros(51)
            eH = np.zeros(51)
            for k in range(0,51):
                #errorbar([th_em_array[k]], [np.mean(data_array[k])], yerr=[np.std(data_array[k])], fmt='.', color='k')
                if not isnan(np.mean(10.**data_array[k])):
                    print np.mean(10.**data_array[k]), np.std(10.**data_array[k]), data_array[k]
                    std_err = np.std(10.**data_array[k])
                    ss = np.sort(10.**data_array[k])
                    cs = np.cumsum(np.ones(len(ss)))/float(len(ss))
                    ss1, ss2 =  get_mean_sigma(ss, conf=0.68)
                    #figure()
                    #plot(ss,cs, '.')
                    #plot([ss1, ss1], [0.,1.], 'r--')
                    #plot([ss2, ss2], [0.,1.], 'r--')
                    #grid(True)
                    #if(np.random.randint(0,5)==0): show()
                    m_ss = np.mean([ss1, ss2])
                    mu[k] = np.mean( get_mean_sigma(ss, conf=0.05))
                    eL[k] = ss1
                    eH[k] = ss2
            plot(th_em_array[1:51], eL[1:], 'b-')
            plot(th_em_array[1:51], eH[1:], 'b-')
            plot(th_em_array[1:51], mu[1:], 'k-')
            fill_between(th_em_array[1:51], eL[1:], eH[1:], alpha=0.3)
            ylim(1.e14, 1.e21)
            '''
            figure(2)
            ax2=subplot(gs[count3*2])
            pcolormesh(th_em_grid, E_tau_grid, img, norm=LogNorm(1.e-7,0.05), cmap='magma')
            xticks(fontsize=14)
            yticks(fontsize=14)
            grid(True)
            ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
            if(count3==4): xlabel('Emergence Angle, deg', fontsize=14)

            index_array = np.concatenate([range(0,51)[::5], range(51,51+84)])
            th_em_grid, E_tau_grid = np.meshgrid(th_em_array[index_array], arange(14.,21.05,0.1))
            #print len(index_array)
            img = np.zeros((len(index_array),len(arange(14.,21.1,0.1))-1))
            #print img.shape
            for k in range(0,len(index_array)):
                #print th_em_array[k]
                hist, bin_edges = np.histogram(data_array[index_array[k]], bins=arange(14.,21.1,0.1))
                img[k,:] = hist/1.e6
                #semilogy(bin_edges[1:], hist+1.e-3)
                #ylim(0.5,1.e5)
            img = np.rot90(img)
            img = np.flipud(img)
            figure(3)
            #ax4=subplot(gs[count3*2+1])
            ax4=subplot(5,1,count3+1)
            #ax4.set_yscale('log')
            ax4.set_xscale('log')
            mu = np.zeros(len(th_em_array))
            eL68 = np.zeros(len(th_em_array))
            eH68 = np.zeros(len(th_em_array))
            eL95 = np.zeros(len(th_em_array))
            eH95 = np.zeros(len(th_em_array))
            eL99 = np.zeros(len(th_em_array))
            eH99 = np.zeros(len(th_em_array))
            for k in range(0,len(index_array)):
                print th_em_array[k], th_em_array[k]%1.0
                #errorbar([th_em_array[k]], [np.mean(data_array[k])], yerr=[np.std(data_array[k])], fmt='.', color='k')
                #if len(data_array[k])>3 and (abs(th_em_array[k]%1.0)<0.09 or abs(th_em_array[k]%1.0)>0.99 or abs(th_em_array[k]-0.1)<0.01):
                if len(data_array[k])>3:
                    print '%1.2e %1.2e'%(np.mean(10.**data_array[k]), np.std(10.**data_array[k])), len(data_array[k])
                    std_err = np.std(10.**data_array[k])
                    ss = np.sort(10.**data_array[k])
                    cs = np.cumsum(np.ones(len(ss)))/float(len(ss))
                    #plot(ss,cs, '.')
                    #ss1, ss2 =  get_mean_sigma(ss, conf=0.95)
                    #m_ss = np.mean([ss1, ss2])
                    #show()
                    #eL[k] = ss1
                    #eH[k] = ss2
                    eL68[k], eH68[k] = np.log10(get_mean_sigma(ss, conf=0.68))
                    eL95[k], eH95[k] = np.log10(get_mean_sigma(ss, conf=0.95))
                    eL99[k], eH99[k] = np.log10(get_mean_sigma(ss, conf=0.99))
                    mu[k] = np.log10(np.mean( get_mean_sigma(ss, conf=0.20)))
                    ylim(14.1, 21.)
            cut=eL68!=0.
            #plot(th_em_array[cut], eL68[cut], 'b-')
            #plot(th_em_array[cut], eH68[cut], 'b-')
            plot(th_em_array[cut], mu[cut], 'r-', lw=2)
            fill_between(th_em_array[cut], eL68[cut], eH68[cut], facecolor = 'k', alpha=0.5)
            fill_between(th_em_array[cut], eL95[cut], eH95[cut], facecolor = 'k', alpha=0.45)
            #fill_between(th_em_array[cut], eL99[cut], eH99[cut], alpha=0.3)
            xlim(0.1,50.)
            grid(True, which='both')
            ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
            new_ticks = np.array([0.1, 0.3, 1., 3., 10.,  30.,]) 
            xticks(new_ticks, new_ticks, fontsize=16)
            yticks(fontsize=15)
            xlabel('Emergence Angle, degrees')
            ax4.text(0.975, 0.95, r'$\log_{10}(E_{\nu} / eV) = %1.0f$'%np.log10(float(energy)),
                verticalalignment='top', horizontalalignment='right',
                transform=ax4.transAxes,
                color='black', fontsize=18)

            figure(2)
            ax=subplot(gs[count3*2+1])
            pcolormesh(th_em_grid, E_tau_grid, img, norm=LogNorm(1.e-7,0.05), cmap='magma')
            xticks(fontsize=14)
            yticks(fontsize=14)
            grid(True)
            cbar = colorbar()
            cbar.ax.tick_params(labelsize=14) 
            ax.text(0.95, 0.95, r'$\log_{10}(E_{\nu} / eV) = %1.0f$'%np.log10(float(energy)),
                verticalalignment='top', horizontalalignment='right',
                transform=ax.transAxes,
                color='black', fontsize=18)

            if(count3==4): xlabel('Emergence Angle, deg', fontsize=16)
            suptitle('Exiting Tau Lepton Energy Spectra', fontsize=20)
            count3+=1      

            figure(4)
            ax_CC = subplot(3,1,1)
            ax_CC.set_yscale('log')
            ax_CC.set_xscale('log')
            cut = mean_num_CC >= 1.
            plot(th_em_array[cut], mean_num_CC[cut], color = colors[cc], lw=2, label='%s eV'%energy)
            leg = legend(loc=3, fontsize=14, title=r'$\nu_{\tau}$ Energy') #, borderpad=0.1, borderaxespad=0, labelspacing=0.1)
            plt.setp(leg.get_title(),fontsize=16)
            title('Mean Number of Interactions')
            ylim(1., 10.)
            xlim(0.1, 50.)
            grid(True, which='both')
            xticks(new_ticks, new_ticks)
            ylabel('Mean Number')

            ax_NC = subplot(3,1,2)
            ax_NC.set_yscale('log')
            ax_NC.set_xscale('log')
            grid(True, which='both')
            cut = mean_num_NC > 1.e-5
            plot(th_em_array[cut], mean_num_NC[cut], color = colors[cc], lw=2)
            xticks(new_ticks, new_ticks)
            xlim(0.1, 50.)
            ylim(1.e-3, 10.)
            ylabel('Mean Number')

            ax_DK = subplot(3,1,3)
            ax_DK.set_yscale('log')
            ax_DK.set_xscale('log')
            cut = mean_num_decays > 1.e-5
            plot(th_em_array[cut], mean_num_decays[cut], color = colors[cc], lw=2) 
            grid(True, which='both')
            ylim(1.e-3, 10.)
            xlim(0.1, 50.)
            xticks(new_ticks, new_ticks)
            ylabel('Mean Number')
            
        cc+=1

    figure(4)
    ax_CC.text(0.05, 0.95, 'CC Interactions',
        verticalalignment='top', horizontalalignment='left',
        transform=ax_CC.transAxes,
        color='black', fontsize=20)
    ylabel('Mean Number')
    ax_NC.text(0.05, 0.95, 'NC Interaction',
        verticalalignment='top', horizontalalignment='left',
        transform=ax_NC.transAxes,
        color='black', fontsize=20)
    ylabel('Mean Number')
    ax_DK.text(0.05, 0.95, 'Tau Decays',
        verticalalignment='top', horizontalalignment='left',
        transform=ax_DK.transAxes,
        color='black', fontsize=20)
    xlabel('Emergence Angle, deg')
    subplots_adjust(top=0.945, bottom=0.075, right=0.95, hspace=0.2, wspace=0.125)
    savefig('Tau_Interactions.pdf')


    figure(2)
    subplots_adjust(top=0.945, bottom=0.075, right=0.95, hspace=0.2, wspace=0.125)
    savefig('Tau_Energy_Spectra.pdf')
    figure(3)
    subplots_adjust(top=0.915, bottom=0.075, right=0.95, hspace=0.25, wspace=0.125)
    suptitle('Tau Lepton Exiting Energies for Various\nNeutrino Energies and Emergence Angles')
    savefig('Tau_Energy_Distrib.pdf')

    figure(1)
    #subplot(211)
    subplot(111)
    title('Probability of Tau Lepton Exiting for Various\nNeutrino Energies and Emergence Angles')
    #ax.text(0.1, 0.95, 'Ice Thickness: %s km\nCross-Section Model: Middle\nEnergy Loss Model: ALLM'%ice_thickness,
    #    verticalalignment='top', horizontalalignment='left',
    #    transform=ax.transAxes,
    #    color='black', fontsize=18)
    ax.text(0.05, 0.05, 'Ice Thickness: %s km\nCross-Section Model: Middle\nEnergy Loss Model: ALLM'%ice_thickness,
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=18)
    subplots_adjust(top=0.875, bottom=0.1, right=0.95)
    savefig('Tau_Exit_Prob.pdf')


    show()


###########################

if(plot_type =='thickness'):
    colors = cm.hot(np.linspace(0, 1, int(2.*float(len(ice_thickness_list)))))
    figure(10, figsize=(8,6))
    #fig = figure(11, figsize=(8,2.75))
    fig = figure(11, figsize=(8,7))
    figure(12, figsize=(8,6.5))

    cc = 0
    energy = '1e+20'
    for ice_thickness in ice_thickness_list[::-1]:
        print 'ice_thickness', ice_thickness
        figure(10)
        th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays = load_LUT('/home/romerowo/nutau_sim/LUTs/%skm_ice_midCS_stdEL/LUT_%s_eV.npz'%(ice_thickness,energy))

        #subplot(211)
        #semilogy(th_em_array, P_exit, '-', color=colors[cc], label='%s km'%ice_thickness)
        #xlim(0.,90.)
        #ylim(1.e-7,1.)
        lgnd = legend(loc=1, fontsize=18, borderpad=0.1, borderaxespad=0, labelspacing=0.1, title='Ice Thickness')
        grid(True, which='both')
        ax=subplot(111)
        ax.set_xscale('log')
        semilogy(th_em_array[P_exit>0.], P_exit[P_exit>0.], '-', lw=2, color=colors[cc], label='%s km'%ice_thickness)
        lgnd = legend(loc=1, fontsize=18, borderpad=0.1, borderaxespad=0, labelspacing=0.1, title='Ice Thickness')
        new_ticks = np.array([0.1, 0.3, 1., 3., 10.,  30.,]) 
        xticks(new_ticks, new_ticks)#, fontsize=16)
        xlim(0.,50.)

        ylim(1.e-6,1.)
        xlabel('Emergence Angle, degrees')
        ylabel('Probability of Tau Exit')
        grid(True, which='both')
        xlabel('Emergence Angle, degrees')

        '''
        figure(11)
        ax2=subplot(1,5,5-cc)
        img = np.zeros((51,len(arange(14.,21.1,0.1))-1))
        th_em_grid, E_tau_grid = np.meshgrid(th_em_array[range(0,51)], arange(14.,21.05,0.1))
        for k in range(0,51):
            #print th_em_array[k]
            hist, bin_edges = np.histogram(data_array[k], bins=arange(14.,21.1,0.1))
            img[k,:] = hist/1.e6
            #semilogy(bin_edges[1:], hist+1.e-3)
            #ylim(0.5,1.e5)
        img = np.rot90(img)
        img = np.flipud(img)
        print th_em_grid.shape, E_tau_grid.shape, img.shape
        pcolormesh(th_em_grid, E_tau_grid, img, norm=LogNorm(1.e-7,0.05), cmap='magma')
        xticks(fontsize=13)
        yticks(fontsize=13)
        ylim(15.,20.)

        title('D=%s km'%ice_thickness, fontsize=14)
        grid(True)
        if(cc==4): ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
        if(cc==0): 
            cbar_ax = fig.add_axes([0.93, 0.225, 0.01, 0.465])
            cbar=colorbar(cax = cbar_ax)
            cbar.ax.tick_params(labelsize=14) 
        if(cc==2): xlabel('Emergence Angle, degrees', fontsize=16)
        '''
        figure(11)
        if(cc!=1):
          #ax2=subplot(1,5,5-cc)
          if(cc==0): ax2=subplot(2,2,4)
          if(cc==2): ax2=subplot(2,2,3)
          if(cc==3): ax2=subplot(2,2,2)
          if(cc==4): ax2=subplot(2,2,1)

          img = np.zeros((51,len(arange(14.,21.1,0.1))-1))
          th_em_grid, E_tau_grid = np.meshgrid(th_em_array[range(0,51)], arange(14.,21.05,0.1))
          for k in range(0,51):
            #print th_em_array[k]
            hist, bin_edges = np.histogram(data_array[k], bins=arange(14.,21.1,0.1))
            img[k,:] = hist/1.e6
            #semilogy(bin_edges[1:], hist+1.e-3)
            #ylim(0.5,1.e5)
          img = np.rot90(img)
          img = np.flipud(img)
          print th_em_grid.shape, E_tau_grid.shape, img.shape
          pcolormesh(th_em_grid, E_tau_grid, img, norm=LogNorm(1.e-7,0.05), cmap='magma')
          xticks(fontsize=15)
          yticks(fontsize=15)
          ylim(15.,20.)

          title('D=%s km'%ice_thickness, fontsize=16)
          ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
          xlabel('Emergence Angle, degrees', fontsize=16)
          grid(True)

          if(cc==0): 
            cbar_ax = fig.add_axes([0.91, 0.3, 0.015, 0.4])
            cbar=colorbar(cax = cbar_ax)
            cbar.ax.tick_params(labelsize=17) 

        if(cc!=1):
            figure(12)
            if(cc==0): ax4=subplot(2,2,4)
            if(cc==2): ax4=subplot(2,2,3)
            if(cc==3): ax4=subplot(2,2,2)
            if(cc==4): ax4=subplot(2,2,1)
            ax4.set_xscale('log')
            mu = np.zeros(len(th_em_array))
            eL68 = np.zeros(len(th_em_array))
            eH68 = np.zeros(len(th_em_array))
            eL95 = np.zeros(len(th_em_array))
            eH95 = np.zeros(len(th_em_array))
            eL99 = np.zeros(len(th_em_array))
            eH99 = np.zeros(len(th_em_array))
            for k in range(0,len(th_em_array)):
                print th_em_array[k], th_em_array[k]%1.0
                #errorbar([th_em_array[k]], [np.mean(data_array[k])], yerr=[np.std(data_array[k])], fmt='.', color='k')
                #if len(data_array[k])>3 and (abs(th_em_array[k]%1.0)<0.09 or abs(th_em_array[k]%1.0)>0.99 or abs(th_em_array[k]-0.1)<0.01):
                if len(data_array[k])>3:
                    print '%1.2e %1.2e'%(np.mean(10.**data_array[k]), np.std(10.**data_array[k])), len(data_array[k])
                    std_err = np.std(10.**data_array[k])
                    ss = np.sort(10.**data_array[k])
                    cs = np.cumsum(np.ones(len(ss)))/float(len(ss))
                    #plot(ss,cs, '.')
                    #ss1, ss2 =  get_mean_sigma(ss, conf=0.95)
                    #m_ss = np.mean([ss1, ss2])
                    #show()
                    #eL[k] = ss1
                    #eH[k] = ss2
                    eL68[k], eH68[k] = np.log10(get_mean_sigma(ss, conf=0.68))
                    eL95[k], eH95[k] = np.log10(get_mean_sigma(ss, conf=0.95))
                    eL99[k], eH99[k] = np.log10(get_mean_sigma(ss, conf=0.99))
                    mu[k] = np.log10(np.mean( get_mean_sigma(ss, conf=0.20)))
                    ylim(14.1, 21.)
            cut=eL68!=0.
            #plot(th_em_array[cut], eL68[cut], 'b-')
            #plot(th_em_array[cut], eH68[cut], 'b-')
            plot(th_em_array[cut], mu[cut], 'r-', lw=2)
            fill_between(th_em_array[cut], eL68[cut], eH68[cut], facecolor = 'k', alpha=0.5)
            fill_between(th_em_array[cut], eL95[cut], eH95[cut], facecolor = 'k', alpha=0.45)
            #fill_between(th_em_array[cut], eL99[cut], eH99[cut], alpha=0.3)
            xlim(0.1,50.)
            grid(True, which='both')
            new_ticks = np.array([0.1, 0.3, 1., 3., 10.,  30.,]) 
            xticks(new_ticks, new_ticks, fontsize=16)
            yticks(fontsize=15)
            if(cc==2 or cc==4): ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
            if(cc==0 or cc==2): xlabel('Emergence Angle, degrees')
            ax4.text(0.975, 0.95, r'Ice Thickness %s km'%ice_thickness,
                verticalalignment='top', horizontalalignment='right',
                transform=ax4.transAxes,
                color='black', fontsize=18)

        cc+=1

    figure(12)
    subplots_adjust(top=0.875, bottom=0.1, left=0.1,right=0.975, hspace=0.175, wspace=0.15)
    suptitle('Tau Lepton Exiting Energies for Various\nIce Layer Thicknesses', fontsize=20)
    savefig('Tau_Energy_Distrib_Depth.pdf')

    figure(11)
    suptitle('Exiting Tau Lepton Energy Spectra for\nVarious Ice Thicknesses')
    #subplots_adjust(top=0.70, left=0.1, right = 0.925, bottom =0.225, wspace=0.3)
    subplots_adjust(hspace=0.4, wspace=0.3, top = 0.86, left=0.1, right=0.89)

    savefig('Tau_Energy_Spectra_Depth.pdf')
    figure(10)
    title('Probability of Tau Lepton Exiting for\nVarious Ice Thicknesses', fontsize=21)
    subplots_adjust(top=0.875)
    ax.text(0.05, 0.05, 'Neutrino Energy: %s eV\nCross-Section Model: Middle\nEnergy Loss Model: ALLM'%energy,
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=18)
    subplots_adjust(bottom=0.125, right=0.95)
    savefig('Tau_Exit_Prob_Depth.pdf')

    show()
###########################
if(plot_type =='models'):
    linestyles= ['-','--']
    figure(100, figsize=(8,6))
    fig=figure(101, figsize=(8,9.5))
    figure(102, figsize=(8,9.5))

    colors = cm.hot(np.linspace(0, 1, int(2.*6.)))
    cc = 0
    energy = '1e+20'
    ice_thickness = '4.0'
    for sigma_model in ['upp', 'mid','low']:
        for e_loss_model in ['low','std']:
            figure(100)
            th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays  = load_LUT('/home/romerowo/nutau_sim/LUTs/%skm_ice_%sCS_%sEL/LUT_%s_eV.npz'%(ice_thickness, sigma_model, e_loss_model,energy))

            grid(True, which='both')
            ax1=subplot(111)
            ax1.set_xscale('log')
            if sigma_model=='upp': CS_m = 'Upper'
            if sigma_model=='mid': CS_m = 'Middle'
            if sigma_model=='low': CS_m = 'Lower'
            if e_loss_model=='std': EL_m = 'ALLM'
            if e_loss_model=='low': EL_m = 'ASW'
            lab_text = '%s, %s'%(CS_m, EL_m)
            semilogy(th_em_array[P_exit>0.], P_exit[P_exit>0.], '-', lw=2, color=colors[cc], label=lab_text)
            lgnd = legend(loc=1, fontsize=16, borderpad=0.1, borderaxespad=0, labelspacing=0.1, title=r'$\sigma_{\nu N}$, $dE_{\tau}/dx$ model')
            plt.setp(lgnd.get_title(),fontsize=16)
            ylim(1.e-6,1.)
            new_ticks = np.array([0.1, 0.3, 1., 3., 10.,  30.,]) 
            xticks(new_ticks, new_ticks)#, fontsize=16)
            xlim(0.,50.)
            xlabel('Emergence Angle, degrees')
            ylabel('Probability of Tau Exit')
            grid(True, which='both')
            '''
            subplot(212)
            semilogy(th_em_array, P_exit, '-', lw=2, color=colors[cc], label=lab_text)
            #lgnd = legend(loc=1, fontsize=18, borderpad=0.1, borderaxespad=0, labelspacing=0.1, title=r'$\sigma_{\nu N}$, $dE_{\tau}/dx$ model')
            xlim(0.,5.)
            ylim(1.e-3,1.)
            xlabel('Emergence Angle, degrees')
            ylabel('Probability of Tau Exit')
            grid(True, which='both')
            xlabel('Emergence Angle, degrees')
            '''

            figure(101)
            print cc
            ax2=subplot(3,2,cc+1)
            img = np.zeros((51,len(arange(14.,21.1,0.1))-1))
            th_em_grid, E_tau_grid = np.meshgrid(th_em_array[range(0,51)], arange(14.,21.05,0.1))
            for k in range(0,51):
	            #print th_em_array[k]
	            hist, bin_edges = np.histogram(data_array[k], bins=arange(14.,21.1,0.1))
	            img[k,:] = hist/1.e6
	            #semilogy(bin_edges[1:], hist+1.e-3)
	            #ylim(0.5,1.e5)
            img = np.rot90(img)
            img = np.flipud(img)
            print th_em_grid.shape, E_tau_grid.shape, img.shape
            pcolormesh(th_em_grid, E_tau_grid, img, norm=LogNorm(1.e-7,0.05), cmap='magma')
            xticks(fontsize=15)
            yticks(fontsize=15)
            ylim(15.,20.)

            title('%s, %s'%(CS_m, EL_m), fontsize=16)
            grid(True)
            ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=17)
            xlabel('Emergence Angle, degrees', fontsize=15)
            #if(cc==0 or cc==2): ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
                        #if(cc==2 or cc==3): xlabel('Emergence Angle, degrees', fontsize=16)
            if(cc==3):
	            cbar_ax = fig.add_axes([0.91, 0.3, 0.015, 0.4])
	            cbar=colorbar(cax = cbar_ax)
	            cbar.ax.tick_params(labelsize=17)


            figure(102)
            ax2=subplot(3,2,cc+1)
            ax2.set_xscale('log')
            mu = np.zeros(len(th_em_array))
            eL68 = np.zeros(len(th_em_array))
            eH68 = np.zeros(len(th_em_array))
            eL95 = np.zeros(len(th_em_array))
            eH95 = np.zeros(len(th_em_array))
            eL99 = np.zeros(len(th_em_array))
            eH99 = np.zeros(len(th_em_array))
            for k in range(0,len(th_em_array)):
                print th_em_array[k], th_em_array[k]%1.0
                #errorbar([th_em_array[k]], [np.mean(data_array[k])], yerr=[np.std(data_array[k])], fmt='.', color='k')
                #if len(data_array[k])>3 and (abs(th_em_array[k]%1.0)<0.09 or abs(th_em_array[k]%1.0)>0.99 or abs(th_em_array[k]-0.1)<0.01):
                if len(data_array[k])>3:
                    print '%1.2e %1.2e'%(np.mean(10.**data_array[k]), np.std(10.**data_array[k])), len(data_array[k])
                    std_err = np.std(10.**data_array[k])
                    ss = np.sort(10.**data_array[k])
                    cs = np.cumsum(np.ones(len(ss)))/float(len(ss))
                    #plot(ss,cs, '.')
                    #ss1, ss2 =  get_mean_sigma(ss, conf=0.95)
                    #m_ss = np.mean([ss1, ss2])
                    #show()
                    #eL[k] = ss1
                    #eH[k] = ss2
                    eL68[k], eH68[k] = np.log10(get_mean_sigma(ss, conf=0.68))
                    eL95[k], eH95[k] = np.log10(get_mean_sigma(ss, conf=0.95))
                    eL99[k], eH99[k] = np.log10(get_mean_sigma(ss, conf=0.99))
                    mu[k] = np.log10(np.mean( get_mean_sigma(ss, conf=0.20)))
                    ylim(14.1, 21.)
            cut=eL68!=0.
            #plot(th_em_array[cut], eL68[cut], 'b-')
            #plot(th_em_array[cut], eH68[cut], 'b-')
            plot(th_em_array[cut], mu[cut], 'r-', lw=2)
            fill_between(th_em_array[cut], eL68[cut], eH68[cut], facecolor = 'k', alpha=0.5)
            fill_between(th_em_array[cut], eL95[cut], eH95[cut], facecolor = 'k', alpha=0.45)
            #fill_between(th_em_array[cut], eL99[cut], eH99[cut], alpha=0.3)
            xlim(0.1,50.)
            grid(True, which='both')
            new_ticks = np.array([0.1, 0.3, 1., 3., 10.,  30.,]) 
            xticks(new_ticks, new_ticks, fontsize=16)
            yticks(fontsize=15)
            if(cc%2==0): ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
            if(cc==4 or cc==5): xlabel('Emergence Angle, degrees')
            ax2.text(0.05, 0.975, '%s, %s'%(CS_m, EL_m),
                verticalalignment='top', horizontalalignment='left',
                transform=ax2.transAxes,
                color='black', fontsize=18)
            #title('%s, %s'%(CS_m, EL_m), fontsize=16)
            grid(True)
            #ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=17)
            #xlabel('Emergence Angle, degrees', fontsize=15)
            #if(cc==0 or cc==2): ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
            #if(cc==2 or cc==3): xlabel('Emergence Angle, degrees', fontsize=16)

            cc+=1
    figure(102)
    suptitle('Exiting Tau Lepton Energies for Various\nCross-section and Energy Loss Models', fontsize=21)
    subplots_adjust(hspace=0.175, wspace=0.175, bottom=0.075, top = 0.9, left=0.1, right=0.975)
    savefig('Tau_Energy_Models.pdf')
    figure(101)
    suptitle('Exiting Tau Lepton Spectra for Various\nCross-section and Energy Loss Models', fontsize=21)
    subplots_adjust(hspace=0.4, wspace=0.4, bottom=0.075, top = 0.875, left=0.1, right=0.89)
    savefig('Tau_Energy_Spectra_Models.pdf')
    figure(100)
    #subplot(211)
    title('Probability of Tau Lepton Exiting for Various\nCross-section and Energy Loss Models')
    subplots_adjust(top=0.875, bottom=0.1, hspace=0.25, right=0.95)
    ax1.text(0.05, 0.05, 'Neutrino Energy: %s eV\nIce Thickness: %s km'%(energy, ice_thickness),
        verticalalignment='bottom', horizontalalignment='left',	
        transform=ax1.transAxes,
        color='black', fontsize=18)
    savefig('Tau_Exit_Prob_Models.pdf')

    show()	

###########################
if(plot_type =='regen'):
    linestyles= ['-','--']
    figure(100, figsize=(8,6))
    fig=figure(101, figsize=(8,4))
    figure(102, figsize=(8,4))

    colors = cm.hot(np.linspace(0, 1, int(2.*6.)))
    cc = 0
    #energy = '1e+20'
    ice_thickness = '4.0'
    sigma_model = 'mid'
    e_loss_model = 'std'
    for energy in ['1e+20', '1e+19', '1e+18']:
        for regen_model in ['', '_no_regen']:
            figure(100)
            th_em_array, P_exit, data_array, mean_num_CC, mean_num_NC, mean_num_decays = load_LUT('/home/romerowo/nutau_sim/LUTs/%skm_ice_%sCS_%sEL%s/LUT_%s_eV.npz'%(ice_thickness, sigma_model, e_loss_model, regen_model, energy))

            linestyle='-'
            if(regen_model!=''): linestyle='--'
            grid(True, which='both')
            ax1=subplot(111)
            ax1.set_xscale('log')
            if sigma_model=='upp': CS_m = 'Upper'
            if sigma_model=='mid': CS_m = 'Middle'
            if sigma_model=='low': CS_m = 'Lower'
            if e_loss_model=='std': EL_m = 'ALLM'
            if e_loss_model=='low': EL_m = 'ASW'
            yesno = 'yes'
            if(regen_model != ''): yesno='no'
            lab_text = '%1.0f, %s'%(np.log10(float(energy)), yesno)
            semilogy(th_em_array[P_exit>0.], P_exit[P_exit>0.], linestyle, lw=2, color=colors[cc], label=lab_text)
            lgnd = legend(loc=1, fontsize=16, borderpad=0.1, borderaxespad=0, labelspacing=0.1, title=r'$\log_{10}(E_{\nu}/eV)$, Regen')
            plt.setp(lgnd.get_title(),fontsize=16)

            ylim(1.e-6,1.)
            new_ticks = np.array([0.1, 0.3, 1., 3., 10.,  30.,]) 
            xticks(new_ticks, new_ticks)#, fontsize=16)

            xlim(0.,50.)
            xlabel('Emergence Angle, degrees')
            ylabel('Probability of Tau Exit')
            grid(True, which='both')
            '''
            subplot(212)
            semilogy(th_em_array, P_exit, linestyle, lw=2, color=colors[cc], label=lab_text)
            #lgnd = legend(loc=1, fontsize=18, borderpad=0.1, borderaxespad=0, labelspacing=0.1, title=r'$\sigma_{\nu N}$, $dE_{\tau}/dx$ model')
            xlim(0.,10.)
            ylim(1.e-4,1.)
            xlabel('Emergence Angle, degrees')
            ylabel('Probability of Tau Exit')
            grid(True, which='both')
            xlabel('Emergence Angle, degrees')
            '''
            if(cc<2): 
                figure(101)
                print cc
                ax2=subplot(1,2,cc+1)
                img = np.zeros((51,len(arange(14.,21.1,0.1))-1))
                th_em_grid, E_tau_grid = np.meshgrid(th_em_array[range(0,51)], arange(14.,21.05,0.1))
                for k in range(0,51):
                    #print th_em_array[k]
                    hist, bin_edges = np.histogram(data_array[k], bins=arange(14.,21.1,0.1))
                    img[k,:] = hist/1.e6
                    #semilogy(bin_edges[1:], hist+1.e-3)
                    #ylim(0.5,1.e5)
                img = np.rot90(img)
                img = np.flipud(img)
                print th_em_grid.shape, E_tau_grid.shape, img.shape
                pcolormesh(th_em_grid, E_tau_grid, img, norm=LogNorm(1.e-7,0.05), cmap='magma')
                xticks(fontsize=15)
                yticks(fontsize=15)
                ylim(15.,20.)

                #title('%s, %s'%(CS_m, EL_m), fontsize=16)
                title(r'$\log_{10}(E_{\nu}/eV)=$%1.0f, %s'%(np.log10(float(energy)), yesno.replace('yes','Regen').replace('no','No Regen')), fontsize=16)

                grid(True)
                ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=17)
                xlabel('Emergence Angle, degrees', fontsize=15)
                #if(cc==0 or cc==2): ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
                #if(cc==2 or cc==3): xlabel('Emergence Angle, degrees', fontsize=16)
                if(cc==1):
                    cbar_ax = fig.add_axes([0.91, 0.15, 0.015, 0.625])
                    cbar=colorbar(cax = cbar_ax)
                    cbar.ax.tick_params(labelsize=17)

                figure(102)
                print cc
                ax3=subplot(1,2,cc+1)
                ax3.set_xscale('log')
                mu = np.zeros(len(th_em_array))
                eL68 = np.zeros(len(th_em_array))
                eH68 = np.zeros(len(th_em_array))
                eL95 = np.zeros(len(th_em_array))
                eH95 = np.zeros(len(th_em_array))
                eL99 = np.zeros(len(th_em_array))
                eH99 = np.zeros(len(th_em_array))
                for k in range(0,len(th_em_array)):
                    print th_em_array[k], th_em_array[k]%1.0
                    #errorbar([th_em_array[k]], [np.mean(data_array[k])], yerr=[np.std(data_array[k])], fmt='.', color='k')
                    #if len(data_array[k])>3 and (abs(th_em_array[k]%1.0)<0.09 or abs(th_em_array[k]%1.0)>0.99 or abs(th_em_array[k]-0.1)<0.01):
                    if len(data_array[k])>3:
                        print '%1.2e %1.2e'%(np.mean(10.**data_array[k]), np.std(10.**data_array[k])), len(data_array[k])
                        std_err = np.std(10.**data_array[k])
                        ss = np.sort(10.**data_array[k])
                        cs = np.cumsum(np.ones(len(ss)))/float(len(ss))
                        #plot(ss,cs, '.')
                        #ss1, ss2 =  get_mean_sigma(ss, conf=0.95)
                        #m_ss = np.mean([ss1, ss2])
                        #show()
                        #eL[k] = ss1
                        #eH[k] = ss2
                        eL68[k], eH68[k] = np.log10(get_mean_sigma(ss, conf=0.68))
                        eL95[k], eH95[k] = np.log10(get_mean_sigma(ss, conf=0.95))
                        eL99[k], eH99[k] = np.log10(get_mean_sigma(ss, conf=0.99))
                        mu[k] = np.log10(np.mean( get_mean_sigma(ss, conf=0.20)))
                        ylim(14.1, 21.)
                cut=eL68!=0.
                #plot(th_em_array[cut], eL68[cut], 'b-')
                #plot(th_em_array[cut], eH68[cut], 'b-')
                plot(th_em_array[cut], mu[cut], 'r-', lw=2)
                fill_between(th_em_array[cut], eL68[cut], eH68[cut], facecolor = 'k', alpha=0.5)
                fill_between(th_em_array[cut], eL95[cut], eH95[cut], facecolor = 'k', alpha=0.45)
                #fill_between(th_em_array[cut], eL99[cut], eH99[cut], alpha=0.3)
                xlim(0.1,50.)
                grid(True, which='both')
                new_ticks = np.array([0.1, 0.3, 1., 3., 10.,  30.,]) 
                xticks(new_ticks, new_ticks, fontsize=16)
                yticks(fontsize=15)
                if(cc%2==0): ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
                xlabel('Emergence Angle, degrees')
                

                ax3.text(0.05, 0.95, '%s'%(yesno.replace('yes','With Regeneration').replace('no','No Regeneration')),
                    verticalalignment='top', horizontalalignment='left',
                    transform=ax3.transAxes,
                    color='black', fontsize=18)
                #title('%s, %s'%(CS_m, EL_m), fontsize=16)
                grid(True)
                #ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=17)
                #xlabel('Emergence Angle, degrees', fontsize=15)
                #if(cc==0 or cc==2): ylabel(r'log$_{10}$($E_{\tau}$ / eV)', fontsize=16)
                #if(cc==2 or cc==3): xlabel('Emergence Angle, degrees', fontsize=16)

            cc+=1
    figure(102)
    suptitle('Exiting Tau Lepton Energy\nWith and Without Tau Regeneration', fontsize=19)
    subplots_adjust(hspace=0.5, wspace=0.15, bottom=0.175, top = 0.8, left=0.1, right=0.975)
    savefig('Tau_Energy_Distrib_Regen.pdf')

    figure(101)
    suptitle('Exiting Tau Lepton Spectra\nWith and Without Tau Regeneration', fontsize=19)
    subplots_adjust(hspace=0.5, wspace=0.3, bottom=0.15, top = 0.785, left=0.1, right=0.89)
    savefig('Tau_Energy_Regen.pdf')
    figure(100)
    subplot(111)
    title('Probability of Tau Lepton Exiting\n With and Without Tau Regeneration')
    subplots_adjust(top=0.85, bottom=0.1, hspace=0.25, right=0.95)
    ax1.text(0.05, 0.05, 'Ice Thickness: %s km'%(ice_thickness),
        verticalalignment='bottom', horizontalalignment='left',	
        transform=ax1.transAxes,
        color='black', fontsize=18)
    savefig('Tau_Exit_Regen.pdf')

    show()	
