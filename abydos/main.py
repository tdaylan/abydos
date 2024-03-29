import tdpy.util
import ephesos.util
import pexo.main
from tdpy.util import summgene

import os, sys
import emcee

import allesfitter.postprocessing.plot_char

import matplotlib.pyplot as plt
from allesfitter import config

import astropy

import numpy as np

import tessmaps
#from tessmaps.src import get_time_on_silicon as gts

import numpy.random
from scipy import stats
import os, datetime

import emcee


def plot_kosm(gdat):
    
    '''
    Demonstrates Kolmogorov-Smirnov test statistics
    '''
     
    path = gdat.pathdata + 'plt.png'
    
    if not os.path.exists(path):
        numbsamp = 1000
        numbpara = 4
        indxpara = np.arange(numbpara)
        
        shft = 1.5 * np.random.randn(numbpara)
        shrk = 0.5 * np.random.randn(numbpara) + 1.
        
        chan = np.random.randn(numbpara * numbsamp).reshape((numbsamp, numbpara))
        for k in indxpara:
            chan[:, k] = shft[k] + np.random.randn(numbsamp) / shrk[k]
        
        strgtitl = ''
        figr, axis = plt.subplots()
        for k in indxpara:
            plt.hist(chan[:, k], alpha=0.4, lw=10, histtype='stepfilled')
            if k != 0:
                kosm, pval = stats.ks_2samp(chan[:, 0], chan[:, k])
                strgtitl += ' %.5g' % pval
        axis.set_title(strgtitl)
        plt.tight_layout()
        figr.savefig(path)
        plt.close()


def main():
   
    '''
    Pipeline to analyze previously known planets detected by TESS
    '''
    
    # construct global object
    gdat = tdpy.util.gdatstrt()
    
    # set the numpy random seed
    np.random.seed(0)
    
    gdat.strgtimestmp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    print('Previously known characterization pipeline for TESS started at %s...' % gdat.strgtimestmp)
        
    gdat.factmjme = 317.907
    gdat.factmsmj = 1048.
    gdat.factrjre = 11.2
    gdat.factrsrj = 9.95
    
    # paths
    gdat.pathbase = os.environ['ABYDOS_DATA_PATH'] + '/'
    pathimag = gdat.pathbase + 'imag/'
    pathdata = gdat.pathbase + 'data/'
    os.system('mkdir -p %s' % pathimag)
    os.system('mkdir -p %s' % pathdata)
    
    gdat.liststrgmodl = ['orbt']
    gdat.listlablmodl = ['Orbital']
    gdat.liststrgstar = ['WASP-77A', 'WASP-43'] 
    gdat.numbstar = len(gdat.liststrgstar)
    gdat.indxstar = np.arange(gdat.numbstar)
    gdat.listlablstar = gdat.liststrgstar

    print('Will analyze data on the following stars:')
    for strgstar in gdat.liststrgstar:
        print(strgstar)
    print

    # plot ESMs
    dictexarcomp = pexo.main.retr_exarcomp()
    
    esmm = ephesos.util.retr_esmm(dictexarcomp['tmptplan'], dictexarcomp['tmptstar'], dictexarcomp['radiplan'], dictexarcomp['radistar'], \
                                                                                                        dictexarcomp['kmagstar'])
    
    indxplansort = np.argsort(esmm)
    
    figr, axis = plt.subplots(figsize=(12, 6))
    axis.hist(esmm[np.where(np.isfinite(esmm))], 100)
    axis.set_xlabel('ESM')
    axis.set_ylabel(r'N')
    axis.set_yscale('log')
    plt.tight_layout()
    path = pathimag + 'histesmm.pdf'
    print('Writing to %s...' % path)
    plt.savefig(path)
    plt.close()
    
    ### TRAPPIST photometry
    #trap = np.loadtxt(pathdata + 'WASP-62/original_data/TRAPPIST_z_WASP62.txt', skiprows=1)
    #trap = trap[:, :3]
    #trap[:, 0] += 2.45e6
    #np.savetxt(pathdata + 'WASP-62/TRAPPIST.csv', trap, delimiter=',')
    ### EULER photometry
    #eulr = np.loadtxt(pathdata + 'WASP-62/original_data/09-01-12_Wasp-62_20111124_Euler_RG.dat')
    #eulr[:, 0] += 2.4e6
    #np.savetxt(pathdata + 'WASP-62/EULER.csv', eulr, delimiter=',')
    ### WASP photometry
    #wasp = np.loadtxt(pathdata + 'WASP-62/original_data/1SWASPJ054833.59-635918.3_WASP62B.lc', skiprows=109)
    #wasp[:, 0] += 2.45e6
    #indxsort = np.argsort(wasp[:, 0])
    #wasp = wasp[indxsort, :]
    #wasp[:, 1] = 10**(-wasp[:, 1] / 2.5)
    #wasp[:, 2] = wasp[:, 1] * wasp[:, 2] / 1.09
    #for k in range(len(wasp)):
    #    if k != len(wasp) - 1 and wasp[k, 0] == wasp[k+1, 0]:
    #        wasp[k+1, 0] += 1e-8
    #np.savetxt(pathdata + 'WASP-62/WASP.csv', wasp, delimiter=',')
    ### CORALIE RV
    #cora = np.loadtxt(pathdata + 'WASP-62/original_data/Wasp62_CORALIE.dat')
    #cora[:, 0] += 2.45e6
    #cora[:, 1] -= np.mean(cora[:, 1])
    #np.savetxt(pathdata + 'WASP-62/CORALIE.csv', cora, delimiter=',')
    
    
    for k in gdat.indxstar:
        strgexar = gdat.liststrgstar[k]

        #for strgruns in ['wout', 'with']:
        #    dictsett = {}
        #    dictpara = {}
        #    dictsett['fast_fit'] = ['True']
        #    dictsett['mcmc_total_steps'] = ['2000']
        #    dictsett['mcmc_burn_steps'] = ['0']
        #    dictpara['b_epoch'] = ['%.13g' % listepoc[k][0], '1,uniform 0 1e12,$P_b$,$\mathrm{d}$']
        #    dictpara['b_period'] = ['%.13g' % listperi[k][0], '1,uniform 0 1e12,$P_b$,$\mathrm{d}$']
        #    if strgruns == 'wout':
        #        liststrgdata = ['TRAPPIST', 'EULER', 'WASP']
        #        dictsett['inst_phot'] = ['TRAPPIST EULER WASP']
        #    if strgruns == 'with':
        #        liststrgdata = ['SPOC', 'TRAPPIST', 'EULER', 'WASP']
        #        dictsett['inst_phot'] = ['TESS TRAPPIST EULER WASP']
    
        # run the pipeline
        #pexo(strgexar=strgexar)

    # tables
    ## planets
    path = pathdata + 'tablinit.tex'
    objtfile = open(path, 'w')
    objtfile.write('& WASP-62b & WASP-100b & WASP-119b & WASP-126\n')
    objtfile.write('Sectors & 1-4, 6-12 & 1-4, 1-12 & 1-4,7,11 & 1-12\n')
    for k, strgplan in enumerate(gdat.liststrgstar):
        objtfile.write('')
    objtfile.write('\n')
    objtfile.close()

    
    ## violin plots
    gdat.pathbaseesmm = gdat.pathbase + 'ESM_Top10/'
    
    gdat.liststrgdata = ['woutTESS', 'alldata']
    gdat.listlabldata = ['w/o TESS', 'w/ TESS']
    allesfitter.postprocessing.plot_char.plot_char(gdat.pathbaseesmm, gdat.liststrgstar, gdat.liststrgdata, gdat.liststrgmodl, \
                                                        gdat.listlablstar, gdat.listlabldata, gdat.listlablmodl, pathimag)
    
    gdat.liststrgdata = ['woutTESS', 'alldata', 'TESS']
    gdat.lisrlabldata = ['w/o TESS', 'w/ TESS', 'o TESS']
    allesfitter.postprocessing.plot_char.plot_char(gdat.pathbaseesmm, gdat.liststrgstar, gdat.liststrgdata, gdat.liststrgmodl, \
                                                        gdat.listlablstar, gdat.listlabldata, gdat.listlablmodl, pathimag)
    
main()
    
