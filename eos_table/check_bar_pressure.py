exec(open('/home/leonardo/Templates/import_py_var.py').read())

filename = "./baryons/DD2_bar.h5"

nb = read_h5(filename, "nb")
t  = read_h5(filename, "t")
yq = read_h5(filename, "yq")
Q1 = read_h5(filename, "Q1")

nb_test = 2.62063249e-01
t_test  = 3.29291436e+00
ye_test = 1.60587635e-02
ym_test = 2.79000776e-03
yq_test = ye_test + ym_test

P = Q1 * nb[:,np.newaxis,np.newaxis]

def P_interp(Pshift):
    logP = np.log(P+Pshift)

    id_nb = np.argmin(abs(nb - nb_test))
    if (nb[id_nb] > nb_test):   id_nb = id_nb - 1

    id_t = np.argmin(abs(t - t_test))
    if (t[id_t] > t_test):   id_t = id_t - 1

    id_yq = np.argmin(abs(yq - yq_test))
    if (yq[id_yq] > yq_test):   id_yq = id_yq - 1

    w1_n = (np.log(nb_test) - np.log(nb[id_nb])) / (np.log(nb[id_nb+1]) - np.log(nb[id_nb]))
    w0_n = 1. - w1_n

    w1_t = (np.log(t_test) - np.log(t[id_t]))    / (np.log(t[id_t+1]) - np.log(t[id_t]))
    w0_t = 1. - w1_t

    w1_y = (yq_test - yq[id_yq])                 / (yq[id_yq+1] - yq[id_yq])
    w0_y = 1. - w1_y

    logP_test = w0_n * (w0_y * (w0_t * logP[id_nb,id_yq,id_t] +      \
                                w1_t * logP[id_nb,id_yq,id_t+1]) +   \
                        w1_y * (w0_t * logP[id_nb,id_yq+1,id_t] +    \
                                w1_t * logP[id_nb,id_yq+1,id_t+1]))  \
              + w1_n * (w0_y * (w0_t * logP[id_nb+1,id_yq,id_t] +    \
                                w1_t * logP[id_nb+1,id_yq,id_t+1]) + \
                        w1_y * (w0_t * logP[id_nb+1,id_yq+1,id_t] +  \
                                w1_t * logP[id_nb+1,id_yq+1,id_t+1]))

    P_test = (np.exp(logP_test) - Pshift) * MeV * 1.0E+39
    return(P_test)

Pshift = abs(np.amin(P)) + 1.0E-10
print("Pshift = %.2lf --> P = %.8e erg/cm3" %(Pshift, P_interp(Pshift)))

Pshift = 2 * abs(np.amin(P))
print("Pshift = %.2lf --> P = %.8e erg/cm3" %(Pshift, P_interp(Pshift)))

Pshift = 0.
print("Pshift = %.2lf --> P = %.8e erg/cm3" %(Pshift, P_interp(Pshift)))
