# -*- coding: utf-8 -*-

"""
@author: alexl
"""

import numpy as np
import PionPhotoLib as ppl
from matplotlib import pyplot as plt
from matplotlib import rcParams
rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'


def main():
    omegaLab=160
    thetas = np.arange(0, 181, 10)
    nucs=["pp0","nn0"]
    poleData = ppl.SaidPoles()

    fig,axs=plt.subplots(2,1,figsize=(8,8))
    #Values of "ccs" come from experiment
    for j, nuc in enumerate(nucs):
        if nuc=="pp0":
            mNucl=938.272 

            ccs=np.array([
                0.1007E+00,
                0.1027E+00,
                0.1084E+00,
                0.1175E+00,
                0.1293E+00,
                0.1432E+00,
                0.1581E+00,
                0.1729E+00,
                0.1866E+00,
                0.1982E+00,
                0.2070E+00,
                0.2126E+00,
                0.2151E+00,
                0.2151E+00,
                0.2133E+00,
                0.2106E+00,
                0.2080E+00,
                0.2061E+00,
                0.2054E+00])
        else:
            mNucl=939.565
            assert nuc=="nn0"

            ccs=np.array([
                0.2351E+00,  
                0.2362E+00, 
                0.2391E+00,
                0.2430E+00,
                0.2464E+00,
                0.2477E+00,
                0.2456E+00,
                0.2387E+00,
                0.2265E+00,
                0.2091E+00,
                0.1871E+00,
                0.1618E+00,
                0.1348E+00,
                0.1078E+00,
                0.8286E-01,
                0.6152E-01,
                0.4522E-01,
                0.3501E-01,
                0.3154E-01])

        sqrtS=np.sqrt(2*omegaLab*mNucl+mNucl*mNucl)

        print("nuc=",nuc)
        print("sqrtS=", sqrtS)
        print("")
        ccsCalc=[]
        for i, theta in enumerate(thetas):
            S = sqrtS**2
            x = np.cos(theta * np.pi / 180)
            # crossSec = ppl.calcCrossSectionFromRaw(S, x, nuc, poleData)
            crossSec = ppl.calcCrossSection(S, x, nuc, poleData)
            print("crossSec=",crossSec)
            ccsCalc.append(crossSec)
            ccExp=ccs[i]

            diff    = abs(ccExp - crossSec)
            pct_err = 100.0 * diff / ccExp
            print(
                "theta={:4d}, cross sec={:8.6f} ---  Experimental Val={:8.6f}, Diff={:8.6f} --> {:6.3f}% Error"
                .format(theta,crossSec, ccExp, diff, pct_err)
            )
        ccCalc=np.array(ccsCalc)
        titles={0:"Proton, Neutral Pion",1:"Neutron, Neutral Pion"}
        axs[j].set_title(titles[j])
        axs[j].plot(thetas,ccs,c='b',label="Experiment")
        axs[j].plot(thetas,ccCalc,c='r',label="Calculated")
        # plt.legend()
        axs[j].legend()
        if j==0:
            print("\n"+100*"#")
    fig.supxlabel(r"$\theta$")
    fig.supylabel(r"$\mathrm{d} \sigma/\mathrm{d}\Omega$")
    plt.show()




if __name__ == "__main__":
    main()
