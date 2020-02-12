import tdpy.util
import tesstarg.util
import pexo.main
from tdpy.util import summgene

import os, sys
import emcee

import allesfitter.postprocessing.plot_viol

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
        
    # paths
    pathbase = os.environ['KNWNTESS_DATA_PATH'] + '/'
    pathimag = pathbase + 'imag/'
    pathdata = pathbase + 'data/'
    os.system('mkdir -p %s' % pathimag)
    os.system('mkdir -p %s' % pathdata)
    
    gdat.liststrgstar = ['WASP-77A', 'WASP-43'] 
    gdat.numbstar = len(gdat.liststrgstar)
    gdat.indxstar = np.arange(gdat.numbstar)
    #listtici = [ 149603524, 38846515, 388104525, 25155310]
    
    print('Will analyze data on the following stars:')
    for strgstar in gdat.liststrgstar:
        print(strgstar)
    print

    # plot ESMs
    dictexarcomp = pexo.main.retr_exarcomp()
    
    esmm = tesstarg.util.retr_esmm(dictexarcomp['tmptplanequb'], dictexarcomp['tmptstar'], dictexarcomp['radiplan'], dictexarcomp['radistar'], \
                                                                                                        dictexarcomp['kmag'])
    
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
    pathbaseesmm = pathbase + 'ESM_Top10/'
    gdat.liststrgruns = ['woutTESS', 'alldata']
    gdat.lisrlablruns = ['w/o TESS', 'w/ TESS']
    allesfitter.postprocessing.plot_viol.plot_viol(pathbaseesmm, gdat.liststrgstar, gdat.liststrgruns, gdat.lisrlablruns, pathimag)
    
    gdat.liststrgruns = ['woutTESS', 'alldata', 'TESS']
    gdat.lisrlablruns = ['w/o TESS', 'w/ TESS', 'o TESS']
    allesfitter.postprocessing.plot_viol.plot_viol(pathbaseesmm, gdat.liststrgstar, gdat.liststrgruns, gdat.lisrlablruns, pathimag)
    

    ## eccentricity and TTV tests
    numbrunstest = 2
    numbtest = 2
    levi = np.empty((numbrunstest, gdat.numbstar, numbtest))
    indxrunstest = np.arange(numbrunstest)
    indxtest = np.arange(numbrunstest)
    
    for a in indxtest:
        for b in indxrunstest:
            for k, strgstar in enumerate(gdat.liststrgstar):
                
                # make sure the run is complete

                # get the edivences
                levi[b, k, a] = np.random.rand()
    
        deltlevi = levi[1, :, :] - levi[0, :, :]
        figr, axis = plt.subplots(figsize=(12, 6))
        axis.plot(gdat.indxstar, deltlevi)
        axis.set_xlabel(gdat.liststrgstar)
        axis.set_ylabel(r'$\Delta \log$ Z')
        plt.tight_layout()
        if a == 0:
            strg = 'ecce'
        if a == 1:
            strg = 'ttvr'
        path = pathimag + strg + 'levi.pdf'
        print('Writing to %s...' % path)
        plt.savefig(path)
        plt.close()
    
    # plot p values
    ## threshold p value to conclude significant difference between posteriors with and without TESS
    pvalthrs = 1e-6
    
    ## calculate the KS test statistic between the posteriors
    numbparacomp = len(lablparacomp[u])
    pval = np.empty(numbparacomp)
    for j in range(numbparacomp):
        kosm, pval[j] = scipy.stats.ks_2samp([gdat.indxrunsfrst[u]][:, j], chan[gdat.indxrunsseco[u]][:, j])
        kosm, pval[j] = scipy.stats.ks_2samp(chan[gdat.indxrunsfrst[u]][:, j], chan[gdat.indxrunsseco[u]][:, j])
    
    ## find the list of parameters whose posterior with and without TESS are unlikely to be drawn from the same distribution
    figr, axis = plt.subplots(figsize=(12, 5))
    indxparacomp = np.arange(numbparacomp)
    axis.plot(indxparacomp, pval, ls='', marker='o')
    indxparagood = np.where(pval < pvalthrs)[0]
    if indxparagood.size > 0:
        axis.plot(indxparacomp[indxparagood], pval[indxparagood], ls='', marker='o', color='r')
    axis.set_yscale('log')
    axis.set_xticks(indxparacomp)
    axis.set_xticklabels(lablparacomp[u])
    axis.axhline(pvalthrs, ls='--', color='black', alpha=0.3)
    plt.tight_layout()
    path = gdat.pathimag + 'kosm_com%d.%s' % (u, gdat.strgplotextn)
    print('Writing to %s...' % path)
    figr.savefig(path)
    plt.close()
    
main()
    
