import numpy as np
# essential functions
def to_wn(value):
    wn = 10 ** 7 / value
    return wn

def to_wl(value):
    wl = 10 ** 7 / value
    return wl

def to_rad(value):
    rad = value * np.pi /180
    return rad

def to_degree(value):
    degree = value * 180 / np.pi
    return degree

# full expression of sfg angle
def sfg_angle(wn_ir, wn_vis, wn_sfg, a_ir, a_vis, n_ir=1, n_vis=1, n_sfg=1,flag='COP'):
    '''
    flag = 'COP', input IR and VIS beams co-propagate
    flag = 'CTOP', input IR and VIS beams counter-propagate
    
    Wavenumber arguments can be substituted with wavevectors.
    '''
    if flag =='COP':
        angle = np.arcsin((n_ir * wn_ir * np.sin(a_ir) + n_vis * wn_vis * np.sin(a_vis)) / 
                          (n_sfg * wn_sfg))
    elif flag =='CTOP':
        angle = np.arcsin((n_vis * wn_vis * np.sin(a_vis) - n_ir * wn_ir * np.sin(a_ir)) / 
                          (n_sfg * wn_sfg))
    return angle

# for equilateral prism sample geometry
def theta_1(n_1, gamma):
    θ = np.pi / 3 - np.arcsin(1 * np.sin(np.pi / 3 - gamma)/ n_1 )
    return θ

# for flat window sample geometry
def from_air(theta_air, n_substrate, n_air=1):
    θ_1 = np.pi / 2 - np.arcsin((n_air / n_substrate) * np.sin(theta_air))
    return θ_1

def refract(θ_i, ind_i, ind_j):
    cos_theta_j = np.sqrt(1 - (ind_i ** 2) / (ind_j ** 2) * (np.sin(θ_i) ** 2))
    theta_j = np.arccos(cos_theta_j)
    return theta_j 

'''
Functions of Fresnel Parameters
'''
def r_ij_p(n_i, n_j, theta_i, theta_j):
    r = (n_j * np.cos(theta_i) - n_i * np.cos(theta_j)) / (n_j * np.cos(theta_i) + n_i * np.cos(theta_j))
    return r

def r_ij_s(n_i, n_j, theta_i, theta_j):
    r = (n_i * np.cos(theta_i) - n_j * np.cos(theta_j)) / (n_i * np.cos(theta_i) + n_j * np.cos(theta_j))
    return r

def t_ij_p(n_i, n_j, theta_i, theta_j):
    t = (2 * n_i * np.cos(theta_i)) / (n_j * np.cos(theta_i) + n_i * np.cos(theta_j))
    return t

def t_ij_s(n_i, n_j, theta_i, theta_j):
    t = (2 * n_i * np.cos(theta_i)) / (n_i * np.cos(theta_i) + n_j * np.cos(theta_j))
    return t

# =========================================== #
def beta(wavelength, n_2, thickness, theta_2):
    β = 2 * np.pi / wavelength * n_2 * thickness * np.cos(theta_2)
    return β

def delta_ir(ir_wl, vis_wl, n_1_ir, n_2_ir, thickness, theta_1_ir, theta_2_ir, theta_2_sfg):
    Δ = (2 * np.pi * n_2_ir * thickness) / (ir_wl * np.cos(theta_2_ir)) - (2 * np.pi * n_1_ir * thickness) / vis_wl * (np.tan(theta_2_ir) + np.tan(theta_2_sfg)) * np.sin(theta_1_ir)
    return Δ

def delta_vis(vis_wl, n_1_vis, n_2_vis, thickness, theta_1_vis, theta_2_vis, theta_2_sfg):
    Δ = (2 * np.pi * n_2_vis * thickness) / (vis_wl * np.cos(theta_2_vis)) - (2 * np.pi * n_1_vis * thickness) / vis_wl * (np.tan(theta_2_vis) + np.tan(theta_2_sfg)) * np.sin(theta_1_vis)
    return Δ

def delta_sfg(sfg_wl, n_2_sfg, thickness, theta_2_sfg):
    Δ = (2 * np.pi * n_2_sfg * thickness) / (sfg_wl * np.cos(theta_2_sfg)) 
    return Δ

# =========================================== #
def L_I_xx(t_12_p, r_12_p, r_23_p, beta, theta_2, theta_1):
    L = (t_12_p) / (1 + r_12_p * r_23_p * np.e **(2j * beta)) * (1 - r_23_p * np.e ** (2j * beta)) * np.cos(theta_2) / np.cos(theta_1)
    return L

def L_I_yy(t_12_s, r_12_s, r_23_s, beta):
    L = (t_12_s) / (1 + r_12_s * r_23_s * np.e **(2j * beta)) * (1 + r_23_s * np.e **(2j * beta))
    return L

def L_I_zz(t_12_p, r_12_p, r_23_p, beta, n_1, n_2, n_int_I):
    L = (t_12_p) / (1 + r_12_p * r_23_p * np.e **(2j * beta)) * (1 + r_23_p * np.e **(2j * beta)) * (n_1 * n_2) / (n_int_I ** 2)
    return L

def L_II_xx(delta, t_12_p, r_12_p, r_23_p, beta, theta_2, theta_1):
    L = np.e ** (1j * delta) * (t_12_p) / (1 + r_12_p * r_23_p * np.e **(2j * beta)) * (1 - r_23_p) * np.cos(theta_2) / np.cos(theta_1)
    return L

def L_II_yy(delta, t_12_s, r_12_s, r_23_s, beta):
    L = np.e ** (1j * delta) * (t_12_s) / (1 + r_12_s * r_23_s * np.e **(2j * beta)) * (1 + r_23_s)
    return L

def L_II_zz(delta, t_12_p, r_12_p, r_23_p, beta, n_1, n_2, n_int_II):
    L = np.e ** (1j * delta) * (t_12_p) / (1 + r_12_p * r_23_p * np.e **(2j * beta)) * (1 + r_23_p) * (n_1 * n_2) / (n_int_II ** 2)
    return L
# =========================================== #
def check(r, t, n_1, theta_1, n_2, theta_2):
    R = r * r.conjugate()
    T = t * t.conjugate() * (n_2 * np.cos(theta_2)) / (n_1 * np.cos(theta_1))
    return T + R

'''
SFG intensity Scalars
'''
def ssp(yy_sfg, yy_vis, zz_ir, theta_ir):
    a = np.sin(theta_ir) * np.sin(theta_ir).conjugate()
    scalar = yy_sfg.real * yy_vis.real * zz_ir.real * a.real
    return scalar

def sps(yy_sfg, zz_vis, yy_ir, theta_vis):
    a = np.sin(theta_vis) * np.sin(theta_vis).conjugate()
    scalar = yy_sfg.real * zz_vis.real * yy_ir.real * a.real
    return scalar

def pss(zz_sfg, yy_vis, yy_ir, theta_sfg):
    a = np.sin(theta_sfg) * np.sin(theta_sfg).conjugate()
    scalar = zz_sfg.real * yy_vis.real * yy_ir.real * a.real
    return scalar

def ppp_xxz(xx_sfg,xx_vis,zz_ir,theta_sfg,theta_vis,theta_ir):
    a = np.cos(theta_sfg) * np.cos(theta_sfg).conjugate()
    b = np.cos(theta_vis) * np.cos(theta_vis).conjugate()
    c = np.sin(theta_ir) * np.sin(theta_ir).conjugate()
    scalar = xx_sfg.real*xx_vis.real*zz_ir.real*a.real*b.real*c.real
    return scalar

def ppp_xzx(xx_sfg,zz_vis,xx_ir,theta_sfg,theta_vis,theta_ir):
    a = np.cos(theta_sfg) * np.cos(theta_sfg).conjugate()
    b = np.sin(theta_vis) * np.sin(theta_vis).conjugate()
    c = np.cos(theta_ir) * np.cos(theta_ir).conjugate()
    scalar = xx_sfg.real*zz_vis.real*xx_ir.real*a.real*b.real*c.real
    return scalar

def ppp_zxx(zz_sfg,xx_vis,xx_ir,theta_sfg,theta_vis,theta_ir):
    a = np.sin(theta_sfg) * np.sin(theta_sfg).conjugate()
    b = np.cos(theta_vis) * np.cos(theta_vis).conjugate()
    c = np.cos(theta_ir) * np.cos(theta_ir).conjugate()
    scalar = zz_sfg.real*xx_vis.real*xx_ir.real*a.real*b.real*c.real
    return scalar

def ppp_zzz(zz_sfg,zz_vis,zz_ir,theta_sfg,theta_vis,theta_ir):
    a = np.sin(theta_sfg) * np.sin(theta_sfg).conjugate()
    b = np.sin(theta_vis) * np.sin(theta_vis).conjugate()
    c = np.sin(theta_ir) * np.sin(theta_ir).conjugate()
    scalar = zz_sfg.real*zz_vis.real*zz_ir.real*a.real*b.real*c.real
    return scalar