import tesstarg.util
import os, sys
import numpy as np
from astroquery.mast import Catalogs
import emcee
import allesfitter.postprocessing.plot_viol
import matplotlib.pyplot as plt
from allesfitter import config
from tdpy.util import summgene
import astropy
import ttvr.util

import numpy as np

import tessmaps
#from tessmaps.src import get_time_on_silicon as gts

import matplotlib.pyplot as plt
import numpy.random
from scipy import stats
import os, datetime
import allesfitter

import emcee
import h5py
from astropy.coordinates import SkyCoord
from allesfitter import config


print 'Known TESS planet pipeline started.'

if len(sys.argv) == 1:
    raise Exception('One argument is expected.')

#pathdata = os.environ['KNWNTESS_DATA_PATH'] + '/'
#liststrgstar = ['WASP-121', 'WASP-62', 'WASP-100', 'WASP-119', 'WASP-126']
#listtici = [ 149603524, 38846515, 388104525, 25155310]
pathdata = os.environ['DATA_PATH'] + '/wasp121b/'
liststrgstar = ['WASP-121']
listtici = ['22529346']
listepoc = [ \
            [2455855.39195, 0.00027], \
            [2456272.3395, 0.0009], \
            [2456537.547, 0.002], \
            [2456890.3191, 0.0006], \
           ]
listperi = [ \
            [4.411953, 0.000003], \
            [2.849375, 0.000008], \
            [2.49979, 0.00001], \
            [3.28880, 0.00001], \
           ]

if sys.argv[1] != 'plot':
    indxruns = liststrgstar.index(sys.argv[1])
    liststrgstar = [liststrgstar[indxruns]]
    listtici = [listtici[indxruns]]
    listepoc = [listepoc[indxruns]]
    listperi = [listperi[indxruns]]

print 'liststrgstar'
print liststrgstar
print 'listepoc'
print listepoc
print 'listperi'
print listperi

print 'Preparing data sets other than TESS...'
# WASP-62
## TRAPPIST photometry
trap = np.loadtxt(pathdata + 'WASP-62/original_data/TRAPPIST_z_WASP62.txt', skiprows=1)
trap = trap[:, :3]
trap[:, 0] += 2.45e6
np.savetxt(pathdata + 'WASP-62/TRAPPIST.csv', trap, delimiter=',')
## EULER photometry
eulr = np.loadtxt(pathdata + 'WASP-62/original_data/09-01-12_Wasp-62_20111124_Euler_RG.dat')
eulr[:, 0] += 2.4e6
np.savetxt(pathdata + 'WASP-62/EULER.csv', eulr, delimiter=',')
## WASP photometry
wasp = np.loadtxt(pathdata + 'WASP-62/original_data/1SWASPJ054833.59-635918.3_WASP62B.lc', skiprows=109)
wasp[:, 0] += 2.45e6
indxsort = np.argsort(wasp[:, 0])
wasp = wasp[indxsort, :]
wasp[:, 1] = 10**(-wasp[:, 1] / 2.5)
wasp[:, 2] = wasp[:, 1] * wasp[:, 2] / 1.09
for k in range(len(wasp)):
    if k != len(wasp) - 1 and wasp[k, 0] == wasp[k+1, 0]:
        wasp[k+1, 0] += 1e-8
np.savetxt(pathdata + 'WASP-62/WASP.csv', wasp, delimiter=',')
## CORALIE RV
cora = np.loadtxt(pathdata + 'WASP-62/original_data/Wasp62_CORALIE.dat')
cora[:, 0] += 2.45e6
cora[:, 1] -= np.mean(cora[:, 1])
np.savetxt(pathdata + 'WASP-62/CORALIE.csv', cora, delimiter=',')


boolstar = False
for k, strgplan in enumerate(liststrgstar):
    if sys.argv[1] == strgplan:
        boolstar = True

if boolstar:
    numbtici = len(listtici)
    for k, tici in enumerate(listtici):
        dictsett = {}
        dictsett['fast_fit'] = ['True']
        dictsett['mcmc_total_steps'] = ['2000']
        dictsett['mcmc_burn_steps'] = ['0']
        dictpara = {}
        dictpara['b_epoch'] = ['%.13g' % listepoc[k][0], '1,uniform 0 1e12,$P_b$,$\mathrm{d}$']
        dictpara['b_period'] = ['%.13g' % listperi[k][0], '1,uniform 0 1e12,$P_b$,$\mathrm{d}$']
        for n in indxruns:
            if n == 0:
                dictpara['b_fc'] = ['0,1,uniform 0 1,$P_b$,$\mathrm{d}$']
                dictpara['b_fc'] = ['0,1,uniform 0 1,$P_b$,$\mathrm{d}$']
            #tesstarg.util.init_alle(tici, dictpara=dictpara, dictsett=dictsett, strgtarg=liststrgstar[k], pathfold=pathdata, strgalleextn='only')

elif sys.argv[1] == 'plot':

    pathdata = os.environ['KNWNTESS_DATA_PATH'] + '/'
    pathpost = pathdata + 'post/'
    os.system('mkdir -p %s' % pathpost)
    
    # tables
    ## planets
    path = pathdata + 'exar.csv'
    objtfile = open(path)
    path = pathdata + 'tablinit.tex'
    objtfile = open(path, 'w')
    objtfile.write('& WASP-62b & WASP-100b & WASP-119b & WASP-126\n')
    objtfile.write('Sectors & 1-4, 6-12 & 1-4, 1-12 & 1-4,7,11 & 1-12\n')
    for k, strgplan in enumerate(liststrgstar):
        objtfile.write('')
    objtfile.write('\n')
    objtfile.close()
    
    ## parameter constraints of WASP-62b
    for strgextn in ['wout', 'with', 'only', 'phas']:
        for strgfile in ['mcmc_fit_b', 'mcmc_corner']:
            cmnd = 'cp %s/WASP-62/allesfit_%s/results/%s.pdf %s/%s_%s.pdf' % (pathdata, strgextn, strgfile, pathpost, strgfile, strgextn)
            print cmnd
            os.system(cmnd)

    ## violin plots
    liststrgrtyp = ['wout', 'with', 'only']
    listpath = []
    for strgrtyp in liststrgrtyp:
        listpath.append(pathdata + 'WASP-62/allesfit_%s/' % strgrtyp)
    
    for strgrtyp in ['wout', 'with']:
        dictsett = {}
        dictpara = {}
        dictsett['fast_fit'] = ['True']
        dictsett['mcmc_total_steps'] = ['2000']
        dictsett['mcmc_burn_steps'] = ['0']
        dictpara['b_epoch'] = ['%.13g' % listepoc[k][0], '1,uniform 0 1e12,$P_b$,$\mathrm{d}$']
        dictpara['b_period'] = ['%.13g' % listperi[k][0], '1,uniform 0 1e12,$P_b$,$\mathrm{d}$']
        if strgrtyp == 'wout':
            liststrgdata = ['TRAPPIST', 'EULER', 'WASP']
            dictsett['inst_phot'] = ['TRAPPIST EULER WASP']
        if strgrtyp == 'with':
            liststrgdata = ['SPOC', 'TRAPPIST', 'EULER', 'WASP']
            dictsett['inst_phot'] = ['TESS TRAPPIST EULER WASP']
        #tesstarg.util.init_alle(listtici[0], dictpara=dictpara, dictsett=dictsett, \
        #                                liststrgdata=liststrgdata, strgtarg=liststrgstar[0], pathfold=pathdata, strgalleextn=strgrtyp)
    
    # WASP-62b
    ## Violins
    #allesfitter.postprocessing.plot_viol.plot_viol(listpath, pathpost, liststrgrtyp=liststrgrtyp)
    
    ## eccentricity
    for k, strgstar in enumerate(liststrgstar):
        pass
        #os.system('cp %s%s/allesfit/results/mcmc_fit_b.pdf %s%s_mcmc_fit_b.pdf' % (pathknwn, strgstar, pathknwn, strgstar))
   
    # TTV
    #pathttvr = pathdata + 'WASP-62/ttvr/'
    #os.system('mkdir -p %s' % pathttvr)
    #pathtmpt = pathdata + 'WASP-62/ttvr/allesfit_tmpt/'
    #liststrgplan = ['b']
    #epoc = [listepoc[0][0]]
    #peri = [listperi[0][0]]
    #liststrginst = ['TESS']
    #offs = [0]
    #ttvr.util.ttvr(pathttvr, pathtmpt, epoc, peri, offs, liststrginst)

    indxepoc = []
    indxperi = []
    numbsamp = np.empty(3, dtype=int)
    numbpara = np.empty(3, dtype=int)
    chan = []
    ## JWST predictions
    for k, strgrtyp in enumerate(liststrgrtyp):
        print 'listpath[k]'
        print listpath[k]
        config.init(listpath[k])
        path = listpath[k] + 'results/mcmc_save.h5'
        emceobjt = emcee.backends.HDFBackend(path, read_only=True)
        chan.append(emceobjt.get_chain())
        numbpara[k] = chan[k].shape[2]
        chan[k] = chan[k].reshape((-1, numbpara[k]))
        numbsamp[k] = chan[k].shape[0]
        indxepoc.append(np.asscalar(np.where(config.BASEMENT.fitkeys == 'b_epoch')[0]))
        indxperi.append(np.asscalar(np.where(config.BASEMENT.fitkeys == 'b_period')[0]))
    
    listyear = [2021, 2023, 2025]
    numbyear = len(listyear)
    timejwst = np.empty((numbyear, numbsamp[1]))
    numbepoc = np.empty((numbyear, 3))
    print 'listepoc[0][0]'
    print listepoc[0][0]
    print 'listperi[0][0]'
    print listperi[0][0]
    print 'chan[1, indxepoc]'
    summgene(chan[1][:, indxepoc[1]])
    print 'chan[1, indxperi]'
    summgene(chan[1][:, indxperi[1]])
    epocjwst = np.empty(numbyear)
    for k, year in enumerate(listyear):
        epocjwst[k] = astropy.time.Time('%d-01-01 00:00:00' % year, format='iso').jd
        for n in range(3):
            #numbepoc[k, n] = (epocjwst[k] - listepoc[0][0]) // listperi[0][0]
            numbepoc[k, n] = (epocjwst[k] - np.mean(chan[n][:, indxepoc[n]])) / np.mean(chan[n][:, indxperi[n]])
            print 'numbepoc[k, n]'
            print numbepoc[k, n]
            numbepoc[k, n] = round(numbepoc[k, n])
            print 'numbepoc[k, n]'
            print numbepoc[k, n]
            print
        timejwst[k, :] = chan[1][:, indxepoc[1]] + chan[1][:, indxperi[1]] * numbepoc[k, 1]# - epocjwst[k]
        timejwst[k, :] -= np.mean(timejwst[k, :])
        timejwst[k, :] *= 24. * 60.
        print 'epocjwst[k]'
        print epocjwst[k]
        print 'numbepoc[k, :]'
        print numbepoc[k, :]
        print 'timejwst[k, :]'
        summgene(timejwst[k, :])
    figr, axis = plt.subplots(figsize=(12, 6))
    axis.violinplot(timejwst.T, listyear)
    axis.set_xlabel('Year')
    axis.set_ylabel('Transit time residual [min]')
    plt.tight_layout()
    path = pathpost + 'jwstfull.pdf'
    print 'Writing to %s...' % path
    plt.savefig(path)
    plt.close()
    
    ## JWST prediction comparison for WASP-62b
    figr, axis = plt.subplots(figsize=(12, 6))
    timejwst = np.empty((3, numbsamp[1]))
    print 'numbsamp'
    print numbsamp

    timejwst = []
    meantimejwst = 0.
    for k, strgrtyp in enumerate(liststrgrtyp):
        #timejwst[k, :] = chan[k][:, indxepoc[k]] + chan[k][:, indxperi[k]] * numbepoc[1]
        print 'strgrtyp'
        print strgrtyp
        print 'numbepoc[1, :]'
        print numbepoc[1, :]
        print 'chan[k][:, indxepoc[k]]'
        summgene(chan[k][:, indxepoc[k]])
        print 'chan[k][:, indxperi[k]]'
        summgene(chan[k][:, indxperi[k]])

        timejwst.append(chan[k][:, indxepoc[k]] + chan[k][:, indxperi[k]] * numbepoc[1, k])
        timejwst[k] -= np.mean(timejwst[k])
        timejwst[k] *= 24. * 60
        print 'np.std(timejwst[k])'
        print np.std(timejwst[k])
        print
    
    axis.violinplot(timejwst, range(len(liststrgrtyp)), points=2000)
    axis.set_xticks(range(len(liststrgrtyp)))
    axis.set_xticklabels(['w/o TESS', 'w/ TESS', 'only TESS'])
    axis.set_ylabel('Transit time residual in 2023 [min]')
    axis.set_ylim([-300, 300])
    plt.tight_layout()
    path = pathpost + 'jwstw62b.pdf'
    print 'Writing to %s...' % path
    plt.savefig(path)
    plt.close()
    
    ## known planets on the map 
   


# plot target locations using tessmaps

def main():
   
    '''
    Pipeline to analyze previously known planets detected by TESS
    '''
    
    # set the numpy random seed
    np.random.seed(0)
    
    # global object
    gdat = gdatstrt()

    # base data path for the TESS known planet project
    gdat.pathdata = os.environ['TESS_KNWN_DATA_PATH']
    os.system('mkdir -p %s' % gdat.pathdata)
    
    os.system('mkdir -p %s' % (gdat.pathdata + '/postproc/'))
    gdat.strgtimestmp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    
    print 'Previously known characterization pipeline for TESS started at %s...' % gdat.strgtimestmp
        
    # list of known planet names
    liststrgplan = [ \
                    #'hats-3', \
                    'hats-13', \
                   ]
    
    print 'Will analyze the following known exoplanets:'
    for strgplan in liststrgplan:
        print strgplan
    print

    # threshold p value to conclude significant difference between posteriors with and without TESS
    pvalthrs = 1e-6
    
    indxcomp = np.arange(3)
    
    # list of rtyp
    # 0 : wout
    # 1 : with
    # 2 : only
    # list of comp
    # wout with
    # wout only
    # with only
    indxrtypfrst = [0, 0, 1]
    indxrtypseco = [1, 2, 2]
    
    liststrgrtyp = ['wout', 'with', 'only'] 
    numbrtyp = len(liststrgrtyp)
    indxrtyp = np.arange(numbrtyp)
    #demo_kosm(gdat)
    
    lpos = [[] for i in indxrtyp]
    chan = [[] for i in indxrtyp]
    lablpara = [[] for i in indxrtyp]
    numbwalk = [[] for i in indxrtyp]
    numbswep = [[] for i in indxrtyp]
    numbsamp = [[] for i in indxrtyp]
    numbburn = [[] for i in indxrtyp]
    factthin = [[] for i in indxrtyp]
    numbpara = [[] for i in indxrtyp]
    numbplan = len(liststrgplan)
    indxplan = np.arange(numbplan)
    for k in indxplan:
        
        print 'Processing %s...' % liststrgplan[k] 
        
        pathdataplan = gdat.pathdata + liststrgplan[k] + '/'
    
        for i, strgrtyp in enumerate(liststrgrtyp):
            pathdata = pathdataplan + '%stess_mcmc/' % strgrtyp
            print 'Reading from %s...' % pathdata
            config.init(pathdata)
            lablpara[i] = config.BASEMENT.fitlabels
            print 'config.BASEMENT.fitkeys'
            print config.BASEMENT.fitkeys

            pathsave = pathdata + 'results/mcmc_save.h5'
            if not os.path.exists(pathsave):
                # sample from the posterior excluding the TESS data
                print 'Calling allesfitter to fit the data...'
                allesfitter.mcmc_fit(pathdata)
            
            emceobjt = emcee.backends.HDFBackend(pathsave, read_only=True)
            chan[i] = emceobjt.get_chain()
            
            # parse configuration
            numbswep[i] = config.BASEMENT.settings['mcmc_total_steps']
            numbburn[i] = config.BASEMENT.settings['mcmc_burn_steps']
            factthin[i] = config.BASEMENT.settings['mcmc_thin_by']
            numbpara[i] = config.BASEMENT.ndim
            
            numbwalk[i] = chan[i].shape[1]
            numbswep[i] = chan[i].shape[0] * factthin
            lpos[i] = emceobjt.get_log_prob()
            numbsamp[i] = lpos[i].size
            chan[i] = chan[i].reshape((-1, numbpara[i]))
            print 'Found %d samples and %d parameters...' % (chan[i].shape[0], chan[i].shape[1])
            print

        lablparacomp = [[] for u in indxcomp]
        for u in indxcomp:
            
            lablparacomp[u] = list(set(lablpara[indxrtypfrst[u]]).intersection(lablpara[indxrtypseco[u]]))

            # post-processing
            ## calculate the KS test statistic between the posteriors
            numbparacomp = len(lablparacomp[u])
            pval = np.empty(numbparacomp)
            for j in range(numbparacomp):
                kosm, pval[j] = stats.ks_2samp(chan[indxrtypfrst[u]][:, j], chan[indxrtypseco[u]][:, j])

            ## find the list of parameters whose posterior with and without TESS are unlikely to be drawn from the same distribution
            indxparagood = np.where(pval < pvalthrs)[0]
            if indxparagood.size > 0:
            
                figr, axis = plt.subplots(figsize=(12, 5))
                indxparacomp = np.arange(numbparacomp)
                axis.plot(indxparacomp, pval, ls='', marker='o')
                axis.plot(indxparacomp[indxparagood], pval[indxparagood], ls='', marker='o', color='r')
                axis.set_yscale('log')
                axis.set_xticks(indxparacomp)
                axis.set_xticklabels(lablparacomp[u])
                if u == 0:
                    axis.set_title('Posteriors with TESS vs. without TESS')
                if u == 1:
                    axis.set_title('Posteriors without TESS vs. only TESS')
                if u == 2:
                    axis.set_title('Posteriors with TESS vs. only TESS')

                axis.axhline(pvalthrs, ls='--', color='black', alpha=0.3)
                plt.tight_layout()
                path = gdat.pathdata + 'postproc/%s_kosm_com%d_%s.png' % (gdat.strgtimestmp, u, liststrgplan[k])
                print 'Writing to %s...' % path
                figr.savefig(path)
                plt.close()
        
        indxrtypmaxm = 1
        
        lablparatotl = np.unique(np.concatenate(lablpara))

        for l, lablparatemp in enumerate(lablparatotl):
            figr, axis = plt.subplots()
            numbbins = 50

            listindxrtypcomp = []
            for i in indxrtyp:
                if lablparatemp in lablpara[i]:
                    listindxrtypcomp.append(i)
            
            #minm = 1e99
            #maxm = -1e99
            #for i in listindxrtypcomp:
            #    m = np.where(lablpara[i] == lablparatemp)[0][0]
            #    minm = min(minm, np.amin(chan[i][:, m]))
            #    maxm = max(maxm, np.amax(chan[i][:, m]))
            #bins = np.linspace(minm, maxm, numbbins + 1)
            ## list of comp
            ## wout with
            ## wout only
            ## with only
            #for i in listindxrtypcomp:
            #    m = np.where(lablpara[i] == lablparatemp)[0][0]
            #    if i == 0:
            #        labl = 'w/o TESS'
            #    if i == 1:
            #        labl = 'with TESS'
            #    if i == 2:
            #        labl = 'only TESS'
            #    axis.hist(chan[i][:, m], alpha=0.4, lw=10, histtype='stepfilled', bins=bins, label=labl)
            #axis.set_xlabel(lablparatemp) 
            #axis.legend()
            #plt.tight_layout()
            #path = gdat.pathdata + 'postproc/pdfn_pr%02d_%s.png' % (l, liststrgplan[k])
            #print 'Writing to %s...' % path
            #figr.savefig(path)
            #plt.close()
            
            figr, axis = plt.subplots()
            chanlist = []
            for i in listindxrtypcomp:
                m = np.where(lablpara[i] == lablparatemp)[0][0]
                
                if lablparatemp == '$T_{0;b}$':
                    offs = np.mean(chan[i][:, m])
                    chantemp = chan[i][:, m] - offs
                    print 'Subtracting %g from T_b for TESS...' % offs
                else:
                    chantemp = chan[i][:, m]

                chanlist.append(chantemp)
            axis.violinplot(chanlist, showmeans=True)
            
            if listindxrtypcomp == [0, 1]:
                labltick = ['w/o TESS', 'w/ TESS']
            if listindxrtypcomp == [0, 2]:
                labltick = ['w/o TESS', 'only TESS']
            if listindxrtypcomp == [1, 2]:
                labltick = ['w/ TESS', 'only TESS']
            if listindxrtypcomp == [0, 1, 2]:
                labltick = ['w/o TESS', 'w/ TESS', 'only TESS']
            valutick = np.arange(len(labltick)) + 1

            axis.set_xticks(valutick)
            axis.set_xticklabels(labltick)
            axis.set_ylabel(lablparatemp)
            plt.tight_layout()

            path = gdat.pathdata + 'postproc/%s_viol_pr%02d_%s.png' % (gdat.strgtimestmp, l, liststrgplan[k])
            print 'Writing to %s...' % path
            figr.savefig(path)
            plt.close()
        
        #figr, axis = plt.subplots()
        #chanlist = []
        #for i, strglabl in []:
        #    m = np.where(lablpara[i] == lablparatemp)[0][0]
        #    if i == 2 and lablparatemp == '$T_{0;b}$':
        #        offs = np.mean(chan[i][2, m] - chan[0][:, m])
        #        chantemp = chan[i][:, m] - offs 
        #        print 'Subtracting %g from T_b for TESS...' % offs
        #    else:
        #        chantemp = chan[i][:, m]

        #    chanlist.append(chantemp)
        #axis.violinplot(chanlist, showmeans=True)
        #
        #if listindxrtypcomp == [0, 1]:
        #    labltick = ['w/o TESS', 'w/ TESS']
        #if listindxrtypcomp == [0, 2]:
        #    labltick = ['w/o TESS', 'only TESS']
        #if listindxrtypcomp == [1, 2]:
        #    labltick = ['w/ TESS', 'only TESS']
        #if listindxrtypcomp == [0, 1, 2]:
        #    labltick = ['w/o TESS', 'w/ TESS', 'only TESS']
        #valutick = np.arange(len(labltick)) + 1

        #axis.set_xticks(valutick)
        #axis.set_xticklabels(labltick)
        #axis.set_ylabel(lablparatemp)
        #plt.tight_layout()

        #path = gdat.pathdata + 'postproc/%s_viol_pr%02d_%s.png' % (gdat.strgtimestmp, l, liststrgplan[k])
        #print 'Writing to %s...' % path
        #figr.savefig(path)
        #plt.close()
        

def retr_listexop(gdat):
   
    pass

    # temp
    #coords = SkyCoord(['124.532 -68.313', '42.42, -42.42'], unit='deg')
    #df = gts.get_time_on_silicon(coords)
    #print 'Time on silicon'
    #print 'df'
    #print df
    

def demo_kosm(gdat):
    
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

allesfitter.post_proc(pathdata)
#demo_kosm()

listindxsect = range(13)
coords = []
names = []
is_transiting = False
title = ''
savname = 'save.pdf'
for indxsect in listindxsect:
    tessmaps.tessmaps.make_rect_map(indxsect, coords, names=names,
                 annotate_bools=is_transiting, title=title,
                 bkgnd_cmap='Blues', savname=savname)



