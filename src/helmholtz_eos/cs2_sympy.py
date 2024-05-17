import numpy as np
import sympy as sp
import sys
from fractions import Fraction

## Define Sympy symbols
eta = sp.symbols("eta", real=True)
theta = sp.symbols("theta", real=True, positive=True)

F = sp.Function("F")
dF_deta = sp.Function("dFdeta")

## Recursive relations for GFF derivatives (Gong+2001) 
## The independent variables are F12, F32, F52, dFdeta_12, dFdeta_32

# wrt theta
def dFdth(k, eta, theta):
    if k == 0.5:
        return sp.expand( (dFdeta(k + 1, eta, theta) - (k + 1) * F(k, eta, theta)) / theta )
    elif k == 1.5:
        return sp.expand( 0.5 * (F(k, eta, theta) - 4. * dFdth(k - 1, eta, theta)) / theta )
    elif k == 2.5:
        return sp.expand( 0.5 * (F(k, eta, theta) - 4. * dFdth(k - 1, eta, theta)) / theta )
    else:
        exit("k must be one between 0.5, 1.5, 2.5, but it is '{k:s}'.")

# wrt eta
def dFdeta(k, eta, theta):
    if k == 0.5:
        return dF_deta(k, eta, theta)
    elif k == 1.5:
        return dF_deta(k, eta, theta)
    elif k == 2.5:
        return sp.expand( 2. * ((k - 1) * F(k - 2, eta, theta) + (0.5 * (k - 1) + 0.75) * theta * F(k - 1, eta, theta) - dF_deta(k - 1, eta, theta)) / theta )
    else:
        exit("k must be one between 0.5, 1.5, 2.5, but it is '{k:s}'.")

## Redefine names of GFFs and GFF derivatives
FDI_names = {}
FDI_deta_names = {}

for i in range(1,7,2):
    for arg in (eta, - eta - 2. / theta):
        k = Fraction(i,2)

        # GFFs
        name = f"f_{str(Fraction(i,2)).replace('/',''):s}"
        if arg == - eta - 2. / theta:
            name = "a_" + name
        FDI_names[F(float(k), arg, theta)] = sp.symbols(name, real=True, positive=True)
        
        # GFF derivatives
        name = f"f_{str(Fraction(i,2)).replace('/',''):s}_deta"
        if arg == - eta - 2. / theta:
            name = "a_" + name
        FDI_deta_names[dF_deta(float(k), arg, theta)] = sp.symbols(name, real=True)
        

## Positron degeneracy parameter
def eta_pos(eta, theta):
    return - eta - 2. / theta

########
# dPdT #
########

def dPdT_partial(eta, theta):
    a_eta = eta_pos(eta, theta)
    
    expr_32 = F(1.5, eta, theta) + F(1.5, a_eta, theta)
    expr_52 = F(2.5, eta, theta) + F(2.5, a_eta, theta)
    expr_32_dth = dFdth(1.5, eta, theta) + dFdth(1.5, a_eta, theta)
    expr_52_dth = dFdth(2.5, eta, theta) + dFdth(2.5, a_eta, theta)
    
    return theta**1.5 * (5. * expr_32 + theta * (3.5 * expr_52 + 2. * expr_32_dth) + theta * theta * expr_52_dth + 2. * dFdeta(1.5, a_eta, theta) / theta + 2. * dFdeta(2.5, a_eta, theta))
    
def dPdeta_partial(eta, theta):
    a_eta = eta_pos(eta, theta)
    
    expr_32_deta = dFdeta(1.5, eta, theta) - dFdeta(1.5, a_eta, theta)
    expr_52_deta = dFdeta(2.5, eta, theta) - dFdeta(2.5, a_eta, theta)
    
    return theta**2.5 * (2. * expr_32_deta + theta * expr_52_deta)

def detadT(eta, theta):
    a_eta = eta_pos(eta, theta)
    
    expr_12 = F(0.5, eta, theta) - F(0.5, a_eta, theta)
    expr_32 = F(1.5, eta, theta) - F(1.5, a_eta, theta)
    expr_12_dth = dFdth(0.5, eta, theta) - dFdth(0.5, a_eta, theta)
    expr_32_dth = dFdth(1.5, eta, theta) - dFdth(1.5, a_eta, theta)

    num = theta**0.5 * (1.5 * expr_12 + theta * (2.5 * expr_32 + expr_12_dth) + theta * theta * expr_32_dth - 2. * dFdeta(0.5, a_eta, theta) / theta - 2. * dFdeta(1.5, a_eta, theta))

    expr_12_deta = dFdeta(0.5, eta, theta) + dFdeta(0.5, a_eta, theta)
    expr_32_deta = dFdeta(1.5, eta, theta) + dFdeta(1.5, a_eta, theta)

    den = theta**1.5 * (expr_12_deta + theta * expr_32_deta) #/ mL[id_L];
   
    num = sp.horner(sp.nsimplify(sp.expand(num * theta**0.5).subs(FDI_names).subs(FDI_deta_names)), wrt=theta)
    den = sp.horner(sp.nsimplify(sp.expand(den / theta**1.5).subs(FDI_names).subs(FDI_deta_names)), wrt=theta)
  
    return - num / den / theta**2

def dPdT(eta, theta):
    return dPdT_partial(eta, theta) + dPdeta_partial(eta, theta) * detadT(eta, theta)

########
# dSdT #
########

def dSdT_partial(eta, theta):
    a_eta = eta_pos(eta, theta)
    
    expr_12 = - 1.5 * eta * F(0.5, eta, theta) - (4.5 * a_eta - 3. / theta) * F(0.5, a_eta, theta) 
    expr_32 = 2.5 * (1. - eta * theta) * F(1.5, eta, theta) + (43./6. - 2.5 * a_eta * theta) * F(1.5, a_eta, theta) 
    expr_52 = 10. / 3. * theta * (F(2.5, eta, theta) + F(2.5, a_eta, theta))
    expr_12_dth = - eta * theta * dFdth(0.5, eta, theta) + (10. / 3. - 3. * a_eta * theta) * dFdth(0.5, a_eta, theta)
    expr_32_dth = theta * (5. / 3. - eta * theta) * dFdth(1.5, eta, theta) + theta * (13. / 3. - a_eta * theta) * dFdth(1.5, a_eta, theta)
    expr_52_dth = dFdth(2.5, eta, theta) + dFdth(2.5, a_eta, theta)

    return theta**0.5 * (expr_12 + expr_32 + expr_52 + expr_12_dth + expr_32_dth + 4. / 3. * theta * theta * expr_52_dth - 2. * a_eta * dFdeta(0.5, a_eta, theta) / theta)

def dSdeta_partial(eta, theta):
    a_eta = eta_pos(eta, theta)
    
    expr_12 = F(0.5, eta, theta) - F(0.5, a_eta, theta)
    expr_32 = F(1.5, eta, theta) - F(1.5, a_eta, theta)
    expr_12_deta = - eta * dFdeta(0.5, eta, theta) + a_eta * dFdeta(0.5, a_eta, theta)
    expr_32_deta = (5. / 3. - eta * theta) * dFdeta(1.5, eta, theta) - (5. / 3. - a_eta * theta) * dFdeta(1.5, a_eta, theta)
    expr_52_deta = dFdeta(2.5, eta, theta) - dFdeta(2.5, a_eta, theta)

    return theta**1.5 * (- expr_12 - theta * expr_32 + expr_12_deta + expr_32_deta + 4. / 3. * theta * expr_52_deta)

def dSdT(eta, theta):
    return dsdT_partial(eta, theta) + dsdeta_partial(eta, theta) * detadT(eta, theta)

########
# dPdn #
########

def dPdn(eta, theta):
    a_eta = eta_pos(eta, theta)

    expr_32_deta = dFdeta(1.5, eta, theta) - dFdeta(1.5, a_eta, theta)
    expr_52_deta = dFdeta(2.5, eta, theta) - dFdeta(2.5, a_eta, theta)

    num = theta * (2. * expr_32_deta + theta * expr_52_deta)

    expr_12_deta = dFdeta(0.5, eta, theta) + dFdeta(0.5, a_eta, theta)
    expr_32_deta = dFdeta(1.5, eta, theta) + dFdeta(1.5, a_eta, theta)

    den = expr_12_deta + theta * expr_32_deta

    num = sp.horner(sp.nsimplify(sp.expand(num / theta).subs(FDI_names).subs(FDI_deta_names)), wrt=theta)
    den = sp.horner(sp.nsimplify(sp.expand(den).subs(FDI_names).subs(FDI_deta_names)), wrt=theta)

    return theta * num / den


########
# dSdn #
########

def dSdn(eta, theta):
    a_eta = eta_pos(eta, theta)

    expr_12 = - F(0.5, eta, theta) + F(0.5, a_eta, theta)
    expr_32 =   F(1.5, eta, theta) - F(1.5, a_eta, theta)
    expr_12_deta = - eta * dFdeta(0.5, eta, theta) + a_eta * dFdeta(0.5, a_eta, theta)
    expr_32_deta_1 = dFdeta(1.5, eta, theta) - dFdeta(1.5, a_eta, theta)
    expr_32_deta_2 = eta * dFdeta(1.5, eta, theta) - a_eta * dFdeta(1.5, a_eta, theta)
    expr_52_deta = dFdeta(2.5, eta, theta) - dFdeta(2.5, a_eta, theta)

    num = expr_12 + expr_12_deta + 5. / 3. * expr_32_deta_1 - theta * (expr_32 + expr_32_deta_2 - 4. / 3. * expr_52_deta)

    expr_12_deta = dFdeta(0.5, eta, theta) + dFdeta(0.5, a_eta, theta)
    expr_32_deta = dFdeta(1.5, eta, theta) + dFdeta(1.5, a_eta, theta)

    den = expr_12_deta + theta * expr_32_deta

    num = sp.horner(sp.nsimplify(sp.expand(num*theta).subs(FDI_names).subs(FDI_deta_names)), wrt=theta)
    den = sp.horner(sp.nsimplify(sp.expand(den).subs(FDI_names).subs(FDI_deta_names)), wrt=theta)

    return num / den / theta


## Get simplified expressions using Sympy
dpdt_partial = dPdT_partial(eta, theta).subs(FDI_names).subs(FDI_deta_names)
dsdt_partial = dSdT_partial(eta, theta).subs(FDI_names).subs(FDI_deta_names)
dpdeta_partial = dPdeta_partial(eta, theta).subs(FDI_names).subs(FDI_deta_names)
dsdeta_partial = dSdeta_partial(eta, theta).subs(FDI_names).subs(FDI_deta_names)

dpdt_partial = sp.horner(sp.nsimplify(sp.expand(dpdt_partial * theta**(-0.5))), wrt=theta)
dsdt_partial = sp.horner(sp.nsimplify(sp.expand(dsdt_partial * theta**(1.5))), wrt=theta)
dpdeta_partial = sp.horner(sp.nsimplify(sp.expand(dpdeta_partial * theta**(-2.5))), wrt=theta)
dsdeta_partial = sp.horner(sp.nsimplify(sp.expand(dsdeta_partial * theta**(-0.5))), wrt=theta)

detadt = detadT(eta, theta) * theta**2.

dpdn = dPdn(eta, theta) / theta
dsdn = dSdn(eta, theta) * theta

## Print derivatives
print("detadT * theta**2:")
print(detadt)
print("")

print("dPdT_partial * theta**(-0.5):")
print(dpdt_partial)
print("")

print("dPdeta_partial * theta**(-2.5):")
print(dpdeta_partial)
print("")

print("dSdT_partial * theta**(1.5):")
print(dsdt_partial)
print("")

print("dSdeta_partial * theta**(-0.5):")
print(dsdeta_partial)
print("")

print("dPdn / theta:")
print(dpdn)
print("")

print("dSdn * theta:")
print(dsdn)
print("")

## Print derivatives using reduced expressions
replacements, reduced_exprs = sp.cse([dpdt_partial, dsdt_partial, dpdeta_partial, dsdeta_partial, detadt, dpdn, dsdn])

print("Reduced expressions:")
for r in replacements:
    print("const double", r[0], "=", r[1], ";")
print("")

print("detadT * theta**2              :", reduced_exprs[4], "\n")
print("dPdT_partial * theta**(-0.5)   :", reduced_exprs[0], "\n")
print("dPdeta_partial * theta**(-2.5) :", reduced_exprs[2], "\n")
print("dSdT_partial * theta**(1.5)    :", reduced_exprs[1], "\n")
print("dSdeta_partial * theta**(-0.5) :", reduced_exprs[3], "\n")
print("dPdn / theta                   :", reduced_exprs[5], "\n")
print("dSdn * theta                   :", reduced_exprs[6], "\n")
