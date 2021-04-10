import os
from matplotlib import pyplot as plt
import numpy as np
import time as tm
import warnings

folder_name = 'Data Input'
path = f'/Users/ricoxi/Desktop/Coding/Fresnel Factor/{folder_name}'
os.chdir(path)

def to_wn(value):
    wn = 10 ** 7 / value
    return wn

def to_wl(value):
    wl = 10 ** 7 / value
    return wl

# Appending water optical constants...

wn_water = []
# this is n
n_water = []
# this is k
k_water = []

with open('d2o_optical.txt', 'r') as f:
    for line in f.readlines()[1:]:
        field = line.split(' ')
        wn_water.append(float(field[0]))
        wn_water.append(float(field[5]))
        n_water.append(float(field[2]))
        n_water.append(float(field[7]))
        k_water.append(float(field[1]))
        k_water.append(float(field[6]))

# this is a list tuple containing wavenumbers 
# and their corresponding n,k values
water = [(wn_water[i],n_water[i],k_water[i]) for i in np.arange(len(wn_water))]
water.sort(key = lambda x: x[0])
# print('Water optical constants are appended in the list "water".')

# Appending CaF2 optical constants...

wl_um = []
n_caf2 = []
with open('caf2_optical.txt') as f:
    for line in f.readlines():
        field = line.split('\t')
        wl_um.append(float(field[0]))
        n_caf2.append(float(field[1]))

wl_nm = [i * 1000 for i in wl_um]
wn_caf2 = [to_wn(i) for i in wl_nm]
# This is a list tuple containing caf2's optical constants
caf2 = [(wn_caf2[i], n_caf2[i]) for i in np.arange(len(wn_caf2))]
caf2.sort(key = lambda x: x[0])
# print('caf2 optical constants are appended in the list "caf2".')

# Appending nickel optical constants...

wl_um = []
n_nickel = []
k_nickel = []
with open('nickel_optical.txt') as f:
    for line in f.readlines():
        field = line.split('\t')
        wl_um.append(float(field[0]))
        n_nickel.append(float(field[1]))
        k_nickel.append(float(field[2]))

wl_nm = [i * 1000 for i in wl_um]
wn_nickel = [to_wn(i) for i in wl_nm]
nickel = [(wn_nickel[i],n_nickel[i],k_nickel[i]) for i in np.arange(len(wn_nickel))]
nickel.sort(key = lambda x: x[0])
# print('Nickel optical constants are appended in the list "nickel".')

'''
selecting data from 2000 to 4000 cm-1
'''

def wn_select(data, lo, hi):
    pts = []
    for i in np.arange(len(data)):
        if lo <= data[i][0] <= hi:
            if len(data[i]) >= 3:
                pt = (data[i][0],data[i][1],data[i][2])
            else:
                pt = (data[i][0],data[i][1],0)
            pts.append(pt)
    return pts

lo = 2400.
hi = 4000.
water_new = wn_select(water, lo, hi)
caf2_new = wn_select(caf2, lo, hi)
nickel_new = wn_select(nickel, lo, hi)

# curve_fitting to find common x_axis
'''
IMPORTANT: 

This is just a simple polynomial fit for well-behaved materials
within the selected wavenumber range. One may need a better fit
if the features are more complex in the desired frequency range.

'''
def fit_curve(data, zone, degree):
    wn = [i[0] for i in data]
    n  = [i[1] for i in data]
    k  = [i[2] for i in data]

    wn_n = np.polyfit(wn, n, degree)
    wn_k = np.polyfit(wn, k, degree)
    f_n  = np.poly1d(wn_n)
    f_k  = np.poly1d(wn_k)

    wn_new = zone
    n_new  = f_n(wn_new)
    k_new  = f_k(wn_new)

#     plt.subplot(1, 2, 1)
#     plt.plot(wn, n, 'o', wn_new, n_new)
#     plt.title('wn vs n')
#     plt.subplot(1, 2, 2)
#     plt.plot(wn, k, 'o', wn_new, k_new)
#     plt.title('wn vs k')

    new = [(wn_new[i], n_new[i], k_new[i]) for i in range(len(wn_new))]

    return new

ir_zone = [i[0] for i in water_new]
nickel_fitted_ir = fit_curve(nickel_new, ir_zone, 10)
caf2_fitted_ir = fit_curve(caf2_new, ir_zone, 10)

sfg_lo = 10000
sfg_hi = 20000
caf2_sfg = wn_select(caf2, sfg_lo, sfg_hi)
vis_wl = 795
sfg_zone = [i[0] + to_wn(vis_wl) for i in water_new]
caf2_fitted_sfg = fit_curve(caf2_sfg, sfg_zone, 10)

sfg_lo = 10000
sfg_hi = 20000
nickel_sfg = wn_select(nickel, sfg_lo, sfg_hi)
vis_wl = 800
sfg_zone = [i[0] + to_wn(vis_wl) for i in water_new]
nickel_fitted_sfg = fit_curve(nickel_sfg, sfg_zone, 10)

caf2_vis = (12500.0, 1.4305, 0)
nickel_vis = (12500.0, 2.2180, 4.8925)
water_vis = (12500.0, 1.3290, 0)

water_sfg = caf2_fitted_sfg
water_sfg = [(i[0], 1.33, 0) for i in water_sfg]

m1 = (caf2_fitted_ir, caf2_vis, caf2_fitted_sfg)
m2 = (nickel_fitted_ir, nickel_vis, nickel_fitted_sfg)
m3 = (water_new, water_vis, water_sfg)

class opt_cst_lib():
    def __init__(self, m1, m2, m3):
        self.m1_ir = m1[0]
        self.m2_ir = m2[0]
        self.m3_ir = m3[0]
        self.m1_vis = m1[1]
        self.m2_vis = m2[1]
        self.m3_vis = m3[1]
        self.m1_sfg = m1[2]
        self.m2_sfg = m2[2]
        self.m3_sfg = m3[2]
        
mylib = opt_cst_lib(m1, m2, m3)
print('D2O, nickel, CaF2 data appended.')