from core import Tmatrix, GF, Stable_Bmatrix, Sweep_backward, Sweep_forward
import numpy as np
import numpy.linalg as LA
import scipy.linalg as sla
import sys
# Check if we are in Jupyter
if 'ipykernel' in sys.modules:
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm

def run(N, Ntau, dtau, Nstable, U, t, mu, Nwarm, Nsweep):
    L = N ** 2
    alpha = np.arccosh(np.exp(0.5*dtau*U))
    T = Tmatrix(t,mu, Nx= N, Ny = N).real
    ExpT = sla.expm(-dtau * T)
    Vdiag = np.zeros(L)
    Vdiag[:] = alpha


    phi = np.random.choice([1, -1], size=(Ntau, L))
    #The part of spin up and spin down are not coupled but direct product
    Uup,Dup,Vup = Stable_Bmatrix(phi,ExpT,Vdiag)
    Udn,Ddn,Vdn = Stable_Bmatrix(phi,ExpT,-Vdiag)
    GFup = GF(Vup[0],Dup[0],Uup[0],Uup[-1],Dup[-1],Vup[-1])
    GFdn = GF(Vdn[0],Ddn[0],Udn[0],Udn[-1],Ddn[-1],Vdn[-1])
    sign  = np.zeros(2 * Nsweep)
    Number = np.zeros(2 * Nsweep)

    # Warm up
    print("Start warm up", flush = True)
    for _ in tqdm(range(Nwarm), desc="Warming up"):
        GFup, GFdn, Uup, Dup, Vup, Udn, Ddn, Vdn, phi = Sweep_backward(GFup,GFdn,Uup,Dup,Vup,Udn,Ddn,Vdn,phi,ExpT,Vdiag,Nstable,alpha)
        GFup, GFdn, Uup, Dup, Vup, Udn, Ddn, Vdn, phi = Sweep_forward(GFup,GFdn,Uup,Dup,Vup,Udn,Ddn,Vdn,phi,ExpT,Vdiag,Nstable,alpha)

    # Sweep
    print("Start sweep")
    for i in tqdm(range(Nsweep), desc="Sweep"):
        GFup, GFdn, Uup, Dup, Vup, Udn, Ddn, Vdn, phi = Sweep_backward(GFup,GFdn,Uup,Dup,Vup,Udn,Ddn,Vdn,phi,ExpT,Vdiag,Nstable,alpha)
        
        # Measurement 
        sign[i*2] = np.sign(LA.det(GFup) * LA.det(GFdn))
        Number[i*2] = np.diag(GFup).sum() + np.diag(GFdn).sum()

        GFup, GFdn, Uup, Dup, Vup, Udn, Ddn, Vdn, phi = Sweep_forward(GFup,GFdn,Uup,Dup,Vup,Udn,Ddn,Vdn,phi,ExpT,Vdiag,Nstable,alpha)
        
        # Measurement 
        sign[i*2 + 1] = np.sign(LA.det(GFup) * LA.det(GFdn))
        Number[i*2 + 1] = np.diag(GFup).sum() + np.diag(GFdn).sum()

    return sign, Number