"""Spectral radiance to photopic-weighted linear sRGB, as the sky atlas uses.

Numba-free so that code outside the solver can convert spectra the same way;
render_data.py imports everything here.
"""
import numpy as np

# Analytic approximation to the CIE 1931 2-degree colour matching functions.
# Wyman, Sloan & Shirley (2013), JCGT 2(2), cited in SOURCES.md.
def cmf(w):
    def g(mu,left,right):
        t=(w-mu)*np.where(w<mu,left,right);return np.exp(-.5*t*t)
    return np.array([1.056*g(599.8,.0264,.0323)+.362*g(442.,.0624,.0374)-.065*g(501.1,.0490,.0382),
      .821*g(568.8,.0213,.0247)+.286*g(530.9,.0613,.0322),
      1.217*g(437.,.0845,.0278)+.681*g(459.,.0385,.0725)]).T
M_XYZ_RGB=np.array([[3.2406,-1.5372,-.4986],[-.9689,1.8758,.0415],[.0557,-.2040,1.0570]])
def colour_weights(lam):
    dw=np.zeros_like(lam);dw[1:-1]=(lam[2:]-lam[:-2])*.5;dw[0]=(lam[1]-lam[0])*.5;dw[-1]=(lam[-1]-lam[-2])*.5
    return 683.*dw[:,None]*cmf(lam)@M_XYZ_RGB.T

def to_rgb(spec,lam):return spec@colour_weights(lam)
