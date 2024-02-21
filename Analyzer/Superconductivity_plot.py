import numpy as np
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
import f90nml



def FillDict(pandasFile, y_dict, firstIter):
    for row in range(len(pandasFile.spin)):
        dictKey = (int(pandasFile.spin[row]) , int(pandasFile.neighbor[row]) , int(pandasFile.sublat[row]) , int(pandasFile.orbital[row]))
        if firstIter:
            y_dict[dictKey] = [np.sqrt(pandasFile.gammaR[row]**2 + pandasFile.gammaIm[row]**2)]
            #print(dictKey)
        else:
            y_dict[dictKey].append(np.sqrt(pandasFile.gammaR[row]**2 + pandasFile.gammaIm[row]**2))
   



def GetRunsGamma(runsPath, matchPattern, xKeywords):
    directories = [dir for dir in os.listdir(runsPath) if re.match(matchPattern, dir)]
    #print(directories)

    y_dict = {}
    x = []

    firstIter = True

    for dir in directories:
        filePathGammaConverged = os.path.join(runsPath, dir, 'OutputData', 'Gamma_SC_final.dat')
        filePathGammaIter = os.path.join(runsPath, dir, 'OutputData', 'Gamma_SC_iter.dat')
        namelistPath = os.path.join(runsPath, dir, 'input.nml')
        
        with open(namelistPath) as nmlFile:
            nml = f90nml.read(nmlFile)
            paramsValuesList = []
            for xKey in xKeywords:
                paramsValuesList.append(nml['physical_params'][xKey])
            x.append(tuple(paramsValuesList))
            #print(nml['physical_params'][xKeyword])
        #Gamma is printed in [meV]
        #If simulation converged final file should exists
        if os.path.exists(filePathGammaConverged):        
            currentGamma = pd.read_fwf(filePathGammaConverged, skiprows = 1, colspecs = [(0,6), (7,11), (12, 16), (17,21), (22,36), (37,51)], \
                                       names = ['spin', 'neighbor', 'sublat', 'orbital', 'gammaR', 'gammaIm'], dtype = np.float64)
            FillDict(currentGamma, y_dict, firstIter)   
            
        #If simulation did NOT converge, iteration file should exists
        elif os.path.exists(filePathGammaIter):        
            currentGamma = pd.read_fwf(filePathGammaIter, skiprows = 1, colspecs = [(0,6), (7,11), (12, 16), (17,21), (22,36), (37,51)], names = ['spin', 'neighbor', 'sublat', 'orbital', 'gammaR', 'gammaIm'], dtype = np.float64)
            FillDict(currentGamma, y_dict, firstIter)        
        
        else:
            print("No data in ", dir)
            continue
        firstIter = False
    return (x,y_dict)



def main():
    # J, Gamma = GetRuns(runsPath = 'RUNS_J', matchPattern = 'RUN_N_steps_12k_alphaB_0.3_J_.*', xKeywords = ('j_sc',))
    # print(J)
    # plt.figure()
    # plt.plot(J,Gamma[(1,1,1,1)], '.')
    # plt.xlabel(r'$J$ (meV)')
    # plt.ylabel(r'$\Gamma$ (meV)')
    # plt.title(r'$E_{Fermi} = 1000$ (meV)')
    # plt.savefig('GammaJ.png')

    # Ef, Gamma = GetRuns(runsPath = 'RUNS_EF', matchPattern = 'RUN_N_steps_12k_alphaB_0.3_Ef_.*', xKeywords = ('e_fermi',))
    # plt.figure()
    # plt.plot(Ef,Gamma[(1,1,1,1)], '.')
    # plt.xlabel(r'$E_{Fermi}$ (meV)')
    # plt.ylabel(r'$\Gamma$ (meV)')
    # plt.title(r'$J_{SC} = 200$ (meV)')
    # plt.savefig('GammaEF.png')
    X, Y = GetRunsGamma(runsPath= '/home/jczarnecki/LAO-STO-results/RUNS_map_EF_J', matchPattern= 'RUN_.*', xKeywords=('e_fermi', 'j_sc'))
    key = (2,3,2,3)


    j_tab = np.linspace(100, 200, 11)
    print(j_tab)
    plt.figure()
    for j in j_tab:
        X_plot = []
        Y_plot = []
        for i in range(len(X)):
            if int(X[i][1]) == j:
                X_plot.append(X[i][0])
                Y_plot.append(Y[key][i])
        X_plot, Y_plot = zip(*sorted(zip(X_plot, Y_plot)))       
        plt.plot(X_plot, Y_plot, '-', label = j)
    
    plt.legend(title = r"$J_{SC}$ (meV)")
    plt.xlabel(r'$E_{Fermi}$ (meV)')
    plt.ylabel(r'$\Gamma$ (meV)')
    plt.savefig('../Plots/GammaEF.png')

    ef_tab = np.arange(-1010, -950, 15)
    print(ef_tab)
    plt.figure()
    for ef in ef_tab:
        X_plot = []
        Y_plot = []
        for i in range(len(X)):
            if int(X[i][0]) == ef:
                X_plot.append(X[i][1])
                Y_plot.append(Y[key][i])
        X_plot, Y_plot = zip(*sorted(zip(X_plot, Y_plot)))       
        plt.plot(X_plot, Y_plot, '-', label = ef)

    plt.legend(title = r"$E_{Fermi}$ (meV)")
    plt.xlabel(r'$J_{SC}$ (meV)')
    plt.ylabel(r'$\Gamma$ (meV)')
    plt.savefig('../Plots/GammaJ.png')
   


if __name__ == "__main__":
    main()