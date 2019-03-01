import tesstarg.util
import os, sys
import numpy as np
from astroquery.mast import Catalogs

# get RA and DEC of known exoplanets from the Exoplanet Archive
pathdata = os.environ['KNWNTESS_DATA_PATH'] + '/'
path = pathdata + 'exar.csv'
objtfile = open(path)
liststrgstar = ['WASP-62', 'WASP-100', 'WASP-119', 'WASP-126']
listtici = [149603524, 38846515, 388104525, 25155310]
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
if len(sys.argv) > 1:
    indxruns = liststrgstar.index(sys.argv[1])
    liststrgstar = [liststrgstar[indxruns]]
    listtici = [listtici[indxruns]]
    listepoc = [listepoc[indxruns]]
    listperi = [listperi[indxruns]]

print 'Known TESS planet pipeline started.'
print 'liststrgstar'
print liststrgstar
print 'listepoc'
print listepoc


listrasc = []
listdecl = []
#listtici = []
#listepoc = []
#liststrgsrch = []
#pathticisave = pathdata + 'tici.csv'
#if os.path.exists(pathticisave):
#    arry = np.loadtxt(pathticisave, delimiter=',')
#    listrasc = arry[:, 0]
#    listdecl = arry[:, 1]
#    listtici = arry[:, 2]
#else:
#    objtfilesave = open(pathticisave, 'w')
#    for k, line in enumerate(objtfile):
#        linesplt = line.split(',')
#        if k == 0:
#            continue
#       
#        print 'k'
#        print k
#        rasc = float(linesplt[47])
#        decl = float(linesplt[49])
#        listrasc.append(rasc)
#        listdecl.append(decl)
#
#        strgsrch = '%g %g' % (rasc, decl)
#        print 'strgsrch'
#        print strgsrch
#        if not strgsrch in liststrgsrch:
#            catalogData = Catalogs.query_region(strgsrch, radius='0.1m', catalog = "TIC")
#            if len(catalogData) > 0:
#                tici = int(catalogData[0]['ID'])
#                listtici.append(tici)
#                objtfilesave.write('%g, %g, %d\n' % (rasc, decl, tici))
#            liststrgsrch.append(strgsrch)
#
#    objtfilesave.close()

if len(sys.argv) == 1:
    raise Exception('One argument is expected.')

boolstar = False
for k, strgplan in enumerate(liststrgplan):
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
        tesstarg.util.init_alle(tici, dictpara=dictpara, dictsett=dictsett, strgtarg=liststrgstar[k], pathfold=pathdata)

elif sys.argv[1] == 'plot':

    # table
    ## planets
    pathdata = os.environ['KNWNTESS_DATA_PATH'] + '/'
    path = pathdata + 'tabl.tex'
    objtfile = open(path, 'w')
    objtfile.write('& WASP-62b & WASP-100b & WASP-119b & WASP-126\n')
    objtfile.write('Sectors & 1-4, 6-12 & 1-4, 1-12 & 1-4,7,11 & 1-12\n')
    for k, strgplan in enumerate(liststrgplan):
        objtfile.write('')
    objtfile.write('\n')
    objtfile.close()
    
    ## posteriors
    for k, strgstar in enumerate(liststrgplan):
        pass
        #os.system('cp %s%s/allesfit/results/mcmc_fit_b.pdf %s%s_mcmc_fit_b.pdf' % (pathknwn, strgstar, pathknwn, strgstar))
    
    ## prediction
    
    figr, axis = plt.subplots(figsize=(12, 6))
    print strgpara
    timetran = 'mcmc_save.h5'
    axis.hist(dictpara[strgpara].flatten())
    axis.set_yscale('log')
    axis.set_xlabel(listlablpara[k])
    plt.tight_layout()
    path = gdat.pathdata + 'hist_%s.pdf' % strgpara
    print 'Writing to %s...' % path
    plt.savefig(path)
    plt.close()

    
    
    
    
    # plots
    ## known planets on the map 
    
    ## posteriors





