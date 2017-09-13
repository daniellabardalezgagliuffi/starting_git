#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 17:17:34 2017

@author: daniella
"""

import sys
sys.path
sys.path.append('/Users/daniella/Python/GitHub/splat/')

import pandas as pd
import numpy as np
import splat
from splat import photometry
from splat import empirical
import time
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astropy.table import Table

## piece of code for regenerating library if needed
'''
class Thingy(object):

def __init__(self):
self._complexValue = None

@property
def complexValue(self):
if self._complexValue is None:
# do expensive thing
self._complexValue = <result>
return self._complexValue
'''

### MAKE SINGLE TEMPLATES DATAFRAME - FOR THIS CASE ONLY NEED L2-L7. FUTURE SPEXBINARYFIT WILL INCLUDE FULL LIBRARY. THIS PART IS TIME-CONSUMING.

start = time.time()
prim = splat.getSpectrum(spexspt=['L2','L7'], VLM=True)
end = time.time()
print(end-start)
## 1111.6824870109558 seconds
## 378.6076831817627 new computer

### ONLY ONE SECONDARY NEEDED FOR THIS EXAMPLE
sec = splat.getSpectrum(shortname='1110+0116')
sec[0].fluxCalibrate('UKIDSS J', 16.16)
sec[0].fluxCalibrate('UKIDSS H', 16.20)
sec[0].fluxCalibrate('UKIDSS K', 16.05)

### COMPARING ONLY ONE SPECTRUM AGAINST THE TEMPLATE LIBRARY. SAVING SOME PROPERTIES WITHIN THE SPECTRUM CLASS
w1355 = splat.Spectrum('/Users/daniella/Research/Papers/BASSUC62/J1355-8258_clean_corrfourstar_spexified_noNaNs.fits')
w1355.opt_type = '--'
w1355.nir_type = 'L9pec'
w1355.spex_type = splat.classifyByStandard(w1355)[0]
w1355.library = 'young'
w1355.name = 'WISE J135501.90-825838.9'
w1355.shortname = 'J1355-8258'
w1355.fluxCalibrate('FOURSTAR J', 16.473)
w1355.fluxCalibrate('FOURSTAR H', 15.447)
w1355.fluxCalibrate('FOURSTAR K', 15.007)
c = SkyCoord('13 55 01.90 -82 58 38.9',  unit=(u.hourangle, u.deg))
w1355.ra = c.ra.value
w1355.dec = c.dec.value

#FOURSTAR mags
#J = 16.473 ± 0.05
#H = 15.447 ± 0.035
#Ks = 15.007 ± 0.047
#w1355.jmag = 16.14
#w1355.hmag = 15.31
#w1355.kmag = 14.72


### MAKE SINGLE TEMPLATE DATAFRAME
specdf = pd.DataFrame(index=np.arange(len(prim)+2),columns=['spectra'])
specdf['spectra'] = prim+sec+[w1355]
specdf['component'] = 'prim'
specdf['component'].ix[len(specdf)-2] = 'sec'
specdf['component'].ix[len(specdf)-1] = 'spec'

### CHECK THAT SPECTRA SPAN FULL JHK WAVELENGTH RANGE
out = []
for i in range(len(specdf)):
    if specdf['spectra'].ix[i].wave.value[0] > photometry.filterProperties('2MASS J')['lambda_min'].value:
        out = np.append(out,i)
out = out[0].astype(int)

specdf = specdf.drop(out,0).reset_index().drop('index',1)

specdf['ra'] = specdf['spectra'].map(lambda x: x.ra)
specdf['dec'] = specdf['spectra'].map(lambda x: x.dec)
specdf['optspt'] = specdf['spectra'].map(lambda x: x.opt_type)
specdf['nirspt'] = specdf['spectra'].map(lambda x: x.nir_type)
specdf['spexspt'] = specdf['spectra'].map(lambda x: x.spex_type)
specdf['library'] = specdf['spectra'].map(lambda x: x.library)

specdf = specdf.replace('--',np.nan)

specdf['name'] = specdf['spectra'].map(lambda x: x.name)
specdf['shortname'] = specdf['spectra'].map(lambda x: x.shortname)
specdf['jmag'], specdf['hmag'], specdf['kmag'] = np.nan, np.nan, np.nan
specdf['jmag_e'], specdf['hmag_e'], specdf['kmag_e'] = np.nan, np.nan, np.nan
specdf['ymg_membership'] = np.nan
abdor = pd.Series(['J0019+4614','J0047+6803','J0326-2102','J0355+1133','J1425-3650','J2206-4217','J2244+2043'])
abdormem = np.where(specdf.shortname.isin(abdor) == True)[0]
specdf['ymg_membership'].ix[abdormem] = 'AB Dor' 


### DOWNLOAD 2MASS MAGNITUDES
for i in np.arange(len(specdf)):
    des = i
    c = SkyCoord(specdf['ra'].ix[i],specdf['dec'].ix[i],unit=(u.deg,u.deg),frame='icrs')
    result = Vizier.query_region(c,radius=10*u.arcsec, catalog='2MASS')    
    if 'II/246/out' in result.keys():
        tmass = result['II/246/out']
    else:
        tmass = np.nan
    if isinstance(tmass,Table):
        print('\nSource {} Designation = {} {} match(es) in 2MASS'.format(i+1,des,len(tmass)))
        #many sources found
        n_tmass = len(tmass)
        if len(tmass) > 1:      # take the closest position
            sep = [c.separation(SkyCoord(tmass['RAJ2000'][lp],tmass['DEJ2000'][lp],unit=(u.deg,u.deg))).arcsecond for lp in np.arange(len(tmass))]
            tmass['sep'] = sep
            tmass.sort('sep')
            while len(tmass)>1:
                tmass.remove_row(1)
            else:
                tmass['sep'] = [c.separation(SkyCoord(tmass['RAJ2000'][0],tmass['DEJ2000'][0],unit=(u.deg,u.deg))).arcsecond]
        print(tmass)
        specdf['jmag'][i] = tmass['Jmag'][0]
        specdf['jmag_e'][i] = tmass['e_Jmag'][0]
        specdf['hmag'][i] = tmass['Hmag'][0]
        specdf['hmag_e'][i] = tmass['e_Hmag'][0]
        specdf['kmag'][i] = tmass['Kmag'][0]
        specdf['kmag_e'][i] = tmass['e_Kmag'][0]

### CALIBRATE ALL SPECTRAL FLUX TO THEIR H BAND APPARENT MAGNITUDES
for i in range(len(specdf)):
    if pd.notnull(specdf['hmag'].ix[i]):
        tmp = specdf['spectra'][i].fluxCalibrate('2MASS H', specdf['hmag'][i])

specdf['MJ'] = specdf['spexspt'].map(lambda x: empirical.typeToMag(x, '2MASS J', ref='dupuy')[0])
specdf['MH'] = specdf['spexspt'].map(lambda x: empirical.typeToMag(x, '2MASS H', ref='dupuy')[0])
specdf['MKs'] = specdf['spexspt'].map(lambda x: empirical.typeToMag(x, '2MASS KS', ref='dupuy')[0])
specdf['MJe'] = specdf['spexspt'].map(lambda x: empirical.typeToMag(x, '2MASS J', ref='dupuy')[1])
specdf['MHe'] = specdf['spexspt'].map(lambda x: empirical.typeToMag(x, '2MASS H', ref='dupuy')[1])
specdf['MKse'] = specdf['spexspt'].map(lambda x: empirical.typeToMag(x, '2MASS KS', ref='dupuy')[1])

#specdf['jmag'] = specdf['spectra'].map(lambda x: splat.filterMag(x, '2MASS J')[0])
#specdf['hmag'] = specdf['spectra'].map(lambda x: splat.filterMag(x, '2MASS H')[0])
#specdf['kmag'] = specdf['spectra'].map(lambda x: splat.filterMag(x, '2MASS KS')[0])
#specdf['jmagerr'] = specdf['spectra'].map(lambda x: splat.filterMag(x, '2MASS J')[1])
#specdf['hmagerr'] = specdf['spectra'].map(lambda x: splat.filterMag(x, '2MASS H')[1])
#specdf['kmagerr'] = specdf['spectra'].map(lambda x: splat.filterMag(x, '2MASS Ks')[1])

nospexspt = np.where(specdf['spexspt'].isnull())[0]
specdf['spexspt'].ix[nospexspt] = specdf['spectra'].ix[nospexspt].map(lambda x: splat.classifyByStandard(x)[0])
specdf['spexsptn'] = specdf['spexspt'].map(lambda x: splat.typeToNum(x))
specdf['distmodH'] = specdf['hmag'] - specdf['MH']

### MASK ELEMENTS WITH NAN MAGNITUDES - THOSE CANNOT BE SCALED


### SAVE SINGLE TEMPLATES DATAFRAME AS PICKLE FILE - ONLY NEED TO READ IT IN NEXT TIME
specdf.to_pickle('/Users/daniella/Research/M7L5Sample/BinaryPopSimulations/wise1355_template_sample_J1110.pickle')

### MAKE BINARY TEMPLATES DATAFRAME
bindf = pd.DataFrame(index=np.arange(len(specdf)),columns=['binspec','spt1','combined_spt','library'])
ssp = specdf['spectra'].ix[545]
ssp.scale(10**(0.4*specdf['distmodH'].ix[545]))

for i in range(len(specdf)):
    psp = specdf['spectra'].ix[i]
    psp.scale(10**(0.4*specdf['distmodH'].ix[i]))
#    bindf['binspec'].ix[i] = specdf['spectra'].ix[i]+specdf['spectra'].ix[545]
    bindf['binspec'].ix[i] = psp+ssp

### REPLACE LAST TWO BINARY TEMPLATE SPECTRA WITH SDSS 1110 AND WISE 1355
bindf['binspec'].ix[len(bindf)-2:] = specdf['spectra'].ix[len(specdf)-2:]

bindf['spt1'] = specdf['spexspt']
bindf['combined_spt'] = bindf['binspec'].map(lambda x: splat.classifyByStandard(x)[0])
bindf['library'] = specdf['library']
bindf.to_pickle('/Users/daniella/Research/M7L5Sample/BinaryPopSimulations/wise1355_binary_sample_J1110.pickle')

### COMPARE SPECTRUM TO SINGLE AND BINARY TEMPLATES. QUANTIFY GOODNESS-OF-FIT WITH CHI SQUARE

def spexbinaryfit(sp):
    
    import splat
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    singdf = pd.read_pickle('/Users/daniella/Research/M7L5Sample/BinaryPopSimulations/wise1355_template_sample_J1110.pickle')
    bindf = pd.read_pickle('/Users/daniella/Research/M7L5Sample/BinaryPopSimulations/wise1355_binary_sample_J1110.pickle')
    
    fitregions = [[0.95, 1.35], [1.45, 1.8], [2.0, 2.35]]
    
    sinchisq = np.zeros(len(bindf))*np.nan
    binchisq = np.zeros(len(bindf))*np.nan
    sinscale = np.zeros(len(bindf))*np.nan
    binscale = np.zeros(len(bindf))*np.nan
    
    for i in range(len(bindf)-1):
        if pd.notnull(singdf['hmag'].ix[i]):
            sinchisq[i], sinscale[i] = splat.compareSpectra(sp,singdf['spectra'].ix[i],statistic='chisqr',fit_ranges=fitregions)
            binchisq[i], binscale[i] = splat.compareSpectra(sp,bindf['binspec'].ix[i],statistic='chisqr',fit_ranges=fitregions)
        else:
            sinchisq[i], sinscale[i] = np.nan,np.nan
            binchisq[i], binscale[i] = np.nan,np.nan

    bindf['sinchisq'] = sinchisq
    bindf['binchisq'] = binchisq
    bindf['sinscale'] = sinscale
    bindf['binscale'] = binscale

    ind = bindf['binchisq'].idxmin()

    
    
    return

spexbinaryfit(w1355)

### MAKE PLOTTING ROUTINES FOR SINGLE AND BINARY TEMPLATES

singdf = pd.read_pickle('/Users/daniella/Research/M7L5Sample/BinaryPopSimulations/wise1355_template_sample_J1110.pickle')
bindf = pd.read_pickle('/Users/daniella/Research/M7L5Sample/BinaryPopSimulations/wise1355_binary_sample_J1110.pickle')

ind = bindf['binchisq'].idxmin()
ind2 = np.where(singdf['component'] == 'sec')[0][0]
ind_abdor = bindf['binchisq'][singdf['ymg_membership'] == 'AB Dor'].idxmin()

ysp1 = singdf['spectra'].ix[ind_abdor]
sp1 = singdf['spectra'].ix[ind]
sp2 = singdf['spectra'].ix[ind2]

ysp1plx = 86.45 #± 0.83 mas Faherty et al. 2016
sp2dist = 19.19 #± 0.44 pc Dupuy & Liu 2012

ysp1.scale(10**(0.4*(5*np.log10(1000/ysp1plx)-5)))

sp1.scale(10**(0.4*singdf['distmodH'].ix[ind]))
#sp2.scale(10**(0.4*singdf['distmodH'].ix[ind2]))
sp2.scale(10**(0.4*(5*np.log10(19.19)-5)))
w1355.scale(10**(0.4*singdf['distmodH'].ix[len(singdf)-1]))

bf = sp1+sp2
ybf = ysp1+sp2

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(w1355.wave.value,w1355.flux.value,'k')
ax.plot(bf.wave.value,bf.flux.value,'g')
ax.plot(sp1.wave.value,sp1.flux.value,'r')
ax.plot(sp2.wave.value,sp2.flux.value,'b')
ax.plot(w1355.wave.value,w1355.noise.value,color='gray')
ax.fill_between(w1355.wave.value,0,w1355.noise.value,color='gray')
ax.plot(bf.wave.value,bf.flux.value,'g')
ax.plot(sp1.wave.value,sp1.flux.value,'r')
ax.plot(sp2.wave.value,sp2.flux.value,'b')
plt.annotate('$\chi^2_R =$ '+str((bindf['binchisq'].ix[ind]/bindf['binspec'].ix[ind].dof).round(2)),xy=(0.6,2.6*10**-11))
plt.ylim([0,3.5*10**-11])
plt.legend([w1355.name+' ('+w1355.spex_type+')','Binary template',sp1.name+' ('+sp1.spex_type+')',sp2.name+' ('+sp2.spex_type+')'],loc=2)
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Flux (erg/$\mu$m s cm$^2$)')

axins = zoomed_inset_axes(ax, 2.5, loc=1)
axins.plot(w1355.wave.value,w1355.flux.value,'k')
axins.plot(bf.wave.value,bf.flux.value,'g')
axins.plot(sp1.wave.value,sp1.flux.value,'r')
axins.plot(sp2.wave.value,sp2.flux.value,'b')
axins.plot(w1355.wave.value,w1355.noise.value,color='gray')
axins.fill_between(w1355.wave.value,0,w1355.noise.value,color='gray')
axins.plot(bf.wave.value,bf.flux.value,'g')
axins.plot(sp1.wave.value,sp1.flux.value,'r')
axins.plot(sp2.wave.value,sp2.flux.value,'b')
axins.set_xlim(1.5, 1.75)
axins.set_ylim(1.2*10**-11, 1.65*10**-11)

plt.savefig('testW1355binaryfit.jpg',dpi=1000)


### BINARYFIT WITH AB DOR PRIMARY

fig, ax = plt.subplots(figsize=(8,6))
ax.plot(w1355.wave.value,w1355.flux.value,'k')
ax.plot(ybf.wave.value,ybf.flux.value,'g')
ax.plot(ysp1.wave.value,ysp1.flux.value,'r')
ax.plot(sp2.wave.value,sp2.flux.value,'b')
ax.plot(w1355.wave.value,w1355.noise.value,color='gray')
ax.fill_between(w1355.wave.value,0,w1355.noise.value,color='gray')
ax.plot(ybf.wave.value,ybf.flux.value,'g')
ax.plot(ysp1.wave.value,ysp1.flux.value,'r')
ax.plot(sp2.wave.value,sp2.flux.value,'b')
plt.annotate('$\chi^2_R =$ '+str((bindf['binchisq'].ix[ind_abdor]/bindf['binspec'].ix[ind_abdor].dof).round(2)),xy=(0.6,2.6*10**-11))
plt.ylim([0,3.5*10**-11])
plt.legend([w1355.name+' ('+w1355.spex_type+')','Binary template',ysp1.name+' ('+ysp1.spex_type+')',sp2.name+' ('+sp2.spex_type+')'],loc=2)
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Flux (erg/$\mu$m s cm$^2$)')

axins = zoomed_inset_axes(ax, 2.5, loc=1)
axins.plot(w1355.wave.value,w1355.flux.value,'k')
axins.plot(ybf.wave.value,ybf.flux.value,'g')
axins.plot(ysp1.wave.value,ysp1.flux.value,'r')
axins.plot(sp2.wave.value,sp2.flux.value,'b')
axins.plot(w1355.wave.value,w1355.noise.value,color='gray')
axins.fill_between(w1355.wave.value,0,w1355.noise.value,color='gray')
axins.plot(ybf.wave.value,ybf.flux.value,'g')
axins.plot(ysp1.wave.value,ysp1.flux.value,'r')
axins.plot(sp2.wave.value,sp2.flux.value,'b')
axins.set_xlim(1.5, 1.75)
axins.set_ylim(1.2*10**-11, 1.65*10**-11)

plt.savefig('testW1355_youngbinaryfit.jpg',dpi=1000)
