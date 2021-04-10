import numpy as np

def to_lambda(omega):
    '''Converts omega in cm^-1 to lambda in nm'''
    return (1.0 / omega) * 1e7

def to_wvnm(lambd):
    '''Converts lambda in nm to omega in cm^-1'''
    return 1.0 / (lambd * 1e-7)

def parse_gold_data(data_file):
    omega = []
    n = []
    k = []
    with open(data_file,'r') as f:
        for line in f.readlines():
            L = line.split()
            omega.append(to_wvnm(float(L[0])*1000))
            n.append(float(L[1]))
            k.append(float(L[2]))
    f.close()
    return np.asarray(omega), np.asarray(n), np.asarray(k)

def parse_h2o_data(data_file):
    omega1 = []
    n1 = []
    k1 = []
    omega2 = []
    n2 = []
    k2 = []
    omega = []
    n = []
    k = []
    with open(data_file,'r') as f:
        for line in f.readlines():
            L = line.split()
            omega1.append(float(L[0]))
            n1.append(float(L[2]))
            k1.append(float(L[1]))
            if len(L) > 5: 
                omega2.append(float(L[5]))
                n2.append(float(L[7]))
                k2.append(float(L[6]))
    f.close()
    omega1 = np.asarray(omega1)
    n1 = np.asarray(n1)
    k1 = np.asarray(k1)
    omega2 = np.asarray(omega2)
    n2 = np.asarray(n2)
    k2 = np.asarray(k2)

    omega = np.squeeze(np.concatenate((np.flipud(omega2[:,np.newaxis]),np.flipud(omega1[:,np.newaxis])),axis=0))
    n = np.squeeze(np.concatenate((np.flipud(n2[:,np.newaxis]),np.flipud(n1[:,np.newaxis])),axis=0))
    k = np.squeeze(np.concatenate((np.flipud(k2[:,np.newaxis]),np.flipud(k1[:,np.newaxis])),axis=0))

    return omega, n, k


def n_sapphire(omega):
    '''Returns the complex refractive index of sapphire for a given omega in cm^-1'''
    lambd = to_lambda(omega)
    lambd_um = lambd / 1000
    n = 1.0 + 1.023798*lambd_um**2/(lambd_um**2 - 0.061448212**2) \
            + 1.058264*lambd_um**2/(lambd_um**2 - 0.1106997**2) \
            + 5.280792*lambd_um**2/(lambd_um**2 - 17.92656**2)
    k = 0.0 * omega
    return n + 1j * k

def n_gold(omega):
    from scipy.interpolate import interp1d
    omega_dat, n_dat, k_dat = parse_gold_data('au_refractive_index_data.txt')
    n_fcn = interp1d(omega_dat,n_dat,kind='cubic', fill_value = 'extrapolate')
    k_fcn = interp1d(omega_dat,k_dat,kind='cubic', fill_value = 'extrapolate')
    nval = n_fcn(omega)
    kval = k_fcn(omega)
    return nval + 1j * kval

def n_water(omega):
    from scipy.interpolate import interp1d
    omega_dat, n_dat, k_dat = parse_h2o_data('h2o_refractive_index_data.txt')
    n_fcn = interp1d(omega_dat,n_dat,kind='cubic', fill_value = 'extrapolate')
    k_fcn = interp1d(omega_dat,k_dat,kind='cubic', fill_value = 'extrapolate')
    nval = n_fcn(omega)
    kval = k_fcn(omega)
    return nval + 1j * kval

def main():
    from matplotlib import pyplot as plt

    omega = np.linspace(2500,4000,10000)
    au = np.zeros((len(omega),2))
    h2o = np.zeros((len(omega),2))
    for i in range(len(omega)):
        nval,kval = n_gold(omega[i])
        au[i,0] = nval
        au[i,1] = kval
        nval,kval = n_water(omega[i])
        h2o[i,0] = nval
        h2o[i,1] = kval
    plt.plot(omega,au[:,0],omega,au[:,1],omega,h2o[:,0],omega,h2o[:,1])
    plt.show()

    return


if __name__ == '__main__':
    main()




