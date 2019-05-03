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

    
    
