#!/usr/bin/env python

"""
David Rodriguez
Oct 2, 2012
Use the convergent point analysis on a group of stars
"""

from math import exp, cos, sin, tan, atan, atan2, acos, sqrt
from astropy.io.votable import parse
import sys
import numpy as np


def cpanalysis(rai,deci,pmrai,pmdei,pmra_ei,pmde_ei,raj, dej,vel,sigint,groupd,sigtheta):
    # Function to calculate prob, dkin, rv, theta, and lambda for a given target for a given CP
    rfa = 3.1415/180
    sigint = sigint * 1000 / (groupd * 4.74047)
    sigtheta = sigtheta * rfa

    p1 = sin( (raj - rai)*rfa )
    p2 = cos(deci*rfa)*tan(dej*rfa) - sin(deci*rfa)*cos( (raj-rai)*rfa )

    theta = atan(p1/p2)
    theta = atan2(p1,p2)

    mupar = pmrai*sin(theta) + pmdei*cos(theta)
    mutan = -1*pmrai*cos(theta) + pmdei*sin(theta)

    sigtan2 = (pmra_ei*cos(theta))**2 + (pmde_ei*sin(theta))**2
    sigtan2 += (sigtheta*mupar)**2

    lamb = sin(deci*rfa)*sin(dej*rfa) + cos(deci*rfa)*cos(dej*rfa)*cos((raj-rai)*rfa)
    lamb = acos(lamb) / rfa

    t2 = mutan**2 / (sigtan2 + sigint**2)

    pro = exp(-0.5*t2)

    para = mupar*4.74047 / (vel * sin(lamb*rfa))
    dkin = 1000/para
    vrad = vel*cos(lamb*rfa)

    return pro, dkin, vrad

# NYMG CP data
nymg = list()
nymg.append({
'name':'beta Pic',
'raj':90,
'dej':-28,
'vel':20.8,
'groupd':40,
'sigint':1,
'sigtheta':2
})
nymg.append({
'name':'AB Dor',
'raj':92,
'dej':-47,
'vel':31.2,
'groupd':50,
'sigint':2,
'sigtheta':2
})
nymg.append({
'name': 'Tuc-Hor',
'raj':119,
'dej':-27,
'vel':23.2,
'groupd':50,
'sigint':1,
'sigtheta':2
})
nymg.append({
'name': 'TWA',
'raj':95,
'dej':-26,
'vel':21.6,
'groupd':50,
'sigint':1,
'sigtheta':2
})
nymg.append({
'name':'Car-Near',
'raj':98,
'dej':0,
'vel':31.3,
'groupd':30,
'sigint':2.6,
'sigtheta':2
})
nymg.append({
'name':'Columba',
'raj':106,
'dej':-30,
'vel':26.5,
'groupd':80,
'sigint':1,
'sigtheta':2
})

# Run the code if executed as a file
if __name__ == "__main__":
    if len(sys.argv)==1:
        print """
        cpmethod.py runs the convergent point code on a list of stars
        
        Usage: cpmethod.py INFILE [OUTFILE RUNTYPE]
        INFI    E: Name of VOTable containing RA, Dec, pmRA, pmDE, and PM errors
        OUTFILE: Output text file name (Default: cp_table.txt)
        RUNTYPE: Type of run to execute; 0- return all, 1- return those with Prob>80% in any group (Default: 0)
        """
        sys.exit()

    file = sys.argv[1]

    if len(sys.argv)>2:
        outfile = sys.argv[2]
    else:
        outfile = 'cp_table.txt'

    if len(sys.argv)>3:
        runtype = int(sys.argv[3])
    else:
        runtype = 0

    votable = parse(file)
    t = votable.get_first_table()

    """
    ra = t.array['ra']
    dec = t.array['dec']
    pmra = t.array['pmRA_ppmxl']
    pmde = t.array['pmDE_ppmxl']
    pmra_e = t.array['e_pmRA_ppmxl']
    pmde_e = t.array['e_pmDE_ppmxl']
    name = t.array['designation']
    """

    """
    ra = t.array['ra']
    dec = t.array['dec']
    pmra = t.array['pmra']
    pmde = t.array['pmde']
    pmra_e = t.array['pmra_err']
    pmde_e = t.array['pmde_err']
    name = t.array['Name']
    """


    ra = t.array['ra']
    dec = t.array['dec']
    pmra = t.array['pmra']
    pmde = t.array['pmde']
    pmra_e = t.array['e_pmra']
    pmde_e = t.array['e_pmde']
    name = t.array['designation']


    num = len(ra)

    pra1 = [] #probability array
    rva1 = [] #rv array
    dka1 = [] #distance array

    pra2 = [] #probability array
    rva2 = [] #rv array
    dka2 = [] #distance array

    pra3 = [] #probability array
    rva3 = [] #rv array
    dka3 = [] #distance array

    pra4 = [] #probability array
    rva4 = [] #rv array
    dka4 = [] #distance array

    pra5 = [] #probability array
    rva5 = [] #rv array
    dka5 = [] #distance array

    pra6 = [] #probability array
    rva6 = [] #rv array
    dka6 = [] #distance array

    # Calculate t2 for the test convergent point
    for i in range(len(ra)):

       print ra[i], dec[i]

       # beta Pic
       raj = 90
       dej = -28
       vel = 20.8
       groupd = 40
       sigint = 1
       sigtheta = 2

       pro, dkin, rv = cpanalysis(ra[i],dec[i],pmra[i],pmde[i], pmra_e[i],pmde_e[i],raj, dej,vel,sigint,groupd,sigtheta)

       if dkin<0: pro = pro*-1.

       pra1.append(pro)
       rva1.append(rv)
       dka1.append(dkin)

       # AB Dor
       raj = 92
       dej = -47
       vel = 31.2
       groupd = 50
       sigint = 2
       sigtheta = 2

       pro, dkin, rv = cpanalysis(ra[i],dec[i],pmra[i],pmde[i], pmra_e[i],pmde_e[i],raj, dej,vel,sigint,groupd,sigtheta)

       if dkin<0: pro = pro*-1.

       pra2.append(pro)
       rva2.append(rv)
       dka2.append(dkin)

       # Tuc-Hor
       raj = 119
       dej = -27
       vel = 23.2
       groupd = 50
       sigint = 1
       sigtheta = 2

       pro, dkin, rv = cpanalysis(ra[i],dec[i],pmra[i],pmde[i], pmra_e[i],pmde_e[i],raj, dej,vel,sigint,groupd,sigtheta)

       if dkin<0: pro = pro*-1.

       pra3.append(pro)
       rva3.append(rv)
       dka3.append(dkin)

       # TWA
       raj = 95
       dej = -26
       vel = 21.6
       groupd = 50
       sigint = 1
       sigtheta = 2

       pro, dkin, rv = cpanalysis(ra[i],dec[i],pmra[i],pmde[i], pmra_e[i],pmde_e[i],raj, dej,vel,sigint,groupd,sigtheta)

       if dkin<0: pro = pro*-1.

       pra4.append(pro)
       rva4.append(rv)
       dka4.append(dkin)

       # Car-Near
       raj = 98
       dej = 0
       vel = 31.3
       groupd = 30
       sigint = 2.6
       sigtheta = 2

       pro, dkin, rv = cpanalysis(ra[i],dec[i],pmra[i],pmde[i], pmra_e[i],pmde_e[i],raj, dej,vel,sigint,groupd,sigtheta)

       if dkin<0: pro = pro*-1.

       pra5.append(pro)
       rva5.append(rv)
       dka5.append(dkin)

       # Columba
       raj = 106 #110
       dej = -30 #-29
       vel = 26.5
       groupd = 80
       sigint = 1
       sigtheta = 2

       pro, dkin, rv = cpanalysis(ra[i],dec[i],pmra[i],pmde[i], pmra_e[i],pmde_e[i],raj, dej,vel,sigint,groupd,sigtheta)

       if dkin<0: pro = pro*-1.

       pra6.append(pro)
       rva6.append(rv)
       dka6.append(dkin)

       #print ra[i], dec[i], theta/rfa, mupar, mutan, sigtan2, t2, pro


    ra   = np.array(ra)
    dec  = np.array(dec)
    pmra = np.array(pmra)
    pmde = np.array(pmde)
    name = np.array(name)
    pra1 = np.array(pra1)*100.
    dka1 = np.array(dka1)
    rva1 = np.array(rva1)
    pra2 = np.array(pra2)*100.
    dka2 = np.array(dka2)
    rva2 = np.array(rva2)
    pra3 = np.array(pra3)*100.
    dka3 = np.array(dka3)
    rva3 = np.array(rva3)
    pra4 = np.array(pra4)*100.
    dka4 = np.array(dka4)
    rva4 = np.array(rva4)
    pra5 = np.array(pra5)*100.
    dka5 = np.array(dka5)
    rva5 = np.array(rva5)
    pra6 = np.array(pra6)*100.
    dka6 = np.array(dka6)
    rva6 = np.array(rva6)

    # Clean up some things (distance, probabilities)
    if runtype==1:
        index = (pra1 >= 80) | (pra2 >= 80) | (pra3 >= 80) | (pra4 >= 80) | (pra5 >= 80) | (pra6 >=80)
        ra = ra[index]
        dec = dec[index]
        pmra = pmra[index]
        pmde = pmde[index]
        name = name[index]
        pra1 = pra1[index]
        dka1 = dka1[index]
        rva1 = rva1[index]
        pra2 = pra2[index]
        dka2 = dka2[index]
        rva2 = rva2[index]
        pra3 = pra3[index]
        dka3 = dka3[index]
        rva3 = rva3[index]
        pra4 = pra4[index]
        dka4 = dka4[index]
        rva4 = rva4[index]
        pra5 = pra5[index]
        dka5 = dka5[index]
        rva5 = rva5[index]
        pra6 = pra6[index]
        dka6 = dka6[index]
        rva6 = rva6[index]
        alldata = zip(pra1,pra2,pra3,pra4,pra5,pra6, ra,dec,pmra,pmde,name, dka1,dka2,dka3,dka4,dka5,dka6, rva1,rva2,rva3,rva4,rva5,rva6 )
        alldata.sort()
        alldata.reverse()
        pra1,pra2,pra3,pra4,pra5,pra6, ra,dec,pmra,pmde,name, dka1,dka2,dka3,dka4,dka5,dka6, rva1,rva2,rva3,rva4,rva5,rva6 = zip(*alldata)


    # Make some output too
    f = open(outfile, 'w')
    tout = '#Index\tName\tRA\tDec\tpmRA\tpmDE\t'
    tout += 'Prob_bP\tVrad_bP\tDkin_bP\t'
    tout += 'Prob_ABD\tVrad_ABD\tDkin_ABD\t'
    tout += 'Prob_TH\tVrad_TH\tDkin_TH\t'
    tout += 'Prob_TWA\tVrad_TWA\tDkin_TWA\t'
    tout += 'Prob_CN\tVrad_CN\tDkin_CN\t'
    tout += 'Prob_Col\tVrad_Col\tDkin_Col\n'
    f.write(tout)

    for i in range(len(ra)):

        mtot = sqrt( pmra[i]**2 + pmde[i]**2 )

        tout = '%02d\t"%s"\t%3.10f\t%3.10f\t%4.1f\t%4.1f\t' % ((i+1), name[i], ra[i], dec[i], pmra[i], pmde[i])
        tout += '%03.1f\t%03.1f\t%03.1f\t' % (pra1[i], rva1[i], dka1[i])
        tout += '%03.1f\t%03.1f\t%03.1f\t' % (pra2[i], rva2[i], dka2[i])
        tout += '%03.1f\t%03.1f\t%03.1f\t' % (pra3[i], rva3[i], dka3[i])
        tout += '%03.1f\t%03.1f\t%03.1f\t' % (pra4[i], rva4[i], dka4[i])
        tout += '%03.1f\t%03.1f\t%03.1f\t' % (pra5[i], rva5[i], dka5[i])
        tout += '%03.1f\t%03.1f\t%03.1f\n' % (pra6[i], rva6[i], dka6[i])

        print tout

        f.write(tout)

    f.close()

