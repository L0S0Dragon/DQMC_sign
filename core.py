import numpy as np
import numpy.linalg as LA
from pythtb import tb_model
import sys
# Check if we are in Jupyter
if 'ipykernel' in sys.modules:
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm

def Tmatrix(t,mu, Nx, Ny):

    # set model parameters
    lat=[[1.0,0.0],[0.0,1.0]]
    orb=[[0.0,0.0]]

    # create TB model
    my_model=tb_model(2,2,lat,orb)
    
    # set on-site energies
    my_model.set_onsite([-mu])
    # set hopping terms (one for each nearest neighbor)
    my_model.set_hop(t, 0, 0, [1, 0])
    my_model.set_hop(t, 0, 0, [0, 1])
    # create a 4x4 lattice
    my_ham = my_model.cut_piece(Nx,0,glue_edgs=False)
    my_ham = my_ham.cut_piece(Ny,1,glue_edgs=False)
    # print Hamiltonian
    return my_ham._gen_ham()


def ExpM(phi,ExpT,Vdiag):
    # T = exp(-T*dtau)
    return np.diag(np.exp(phi * Vdiag)).dot(ExpT)


def SVD_stablizer_back(A,U,D,V):
    # B = U * D * V
    # B * A = U * D * V * B = U * U1 * D1 * V1 = U2 * D1 * V1
    U1, D, V = LA.svd(np.diag(D).dot(V).dot(A))
    return U.dot(U1), D, V
    
def SVD_stablizer_forward(A,U,D,V):
    # B = U * D * V
    # A * B = A * U * D * V = U1 * D1 * V1 * V = U1 * D1 * V2
    U, D, V1 = LA.svd(A.dot(U.dot(np.diag(D))))
    return U, D, V1.dot(V)


def GF(UL,DL,VL,UR,DR,VR):
    # B(t,0) = UL * DL * VL
    # B(beta,t) = VR * DR * UR Take care of it !!!
    U,D,V = LA.svd(LA.inv(UL.dot(UR))+np.diag(DR).dot(VR.dot(VL)).dot(np.diag(DL)))
    return LA.inv(V.dot(UL)).dot(np.diag(1/D)).dot(LA.inv(UR.dot(U)))


def Stable_Bmatrix(phi,ExpT,Vdiag):
    # use SVD algorithm to calculate the stable B matrix
    # ExpT: the matrix of exp(-dtau * T)
    # Vdiag: Ntau x L matrix with each column being - alpha 
    Ntau, L = phi.shape
    U = np.zeros((Ntau+1,L,L))
    D = np.zeros((Ntau+1,L))
    V = np.zeros((Ntau+1,L,L))

    U[0],D[0],V[0] = LA.svd(np.eye(L))
    for i in range(Ntau):
        U[i+1], D[i+1], V[i+1] = SVD_stablizer_forward(ExpM(phi[i], ExpT, Vdiag),U[i],D[i],V[i])
    return U,D,V

def Update_HS(GFup, GFdn, phi, alpha):
    Ident = np.eye(len(phi))
    for i in range(len(phi)):
        dup = np.exp(-2 * alpha * phi[i]) -1 
        ddn = np.exp(2 * alpha * phi[i]) -1
        Rup = 1 + dup * (1 - GFup[i,i])
        Rdn = 1 + ddn * (1 - GFdn[i,i])
        
        if  np.abs(Rup * Rdn) > np.random.rand():
            
            GFup -= dup * np.outer(GFup[:,i], Ident[i]-GFup[i,:])/Rup
            GFdn -= ddn * np.outer(GFdn[:,i], Ident[i]-GFdn[i,:])/Rdn
            phi[i] = -phi[i]
    
    return GFup, GFdn, phi


def Sweep_forward(GFup,GFdn,Uup,Dup,Vup,Udn,Ddn,Vdn,phi,ExpT,Vdiag,Nstable,alpha):
    # Input GFup, GFdn should be G(0,0) = G(\beta,\beta)
    # The first member of U, D, V should be the result of I = B(0,0) or B(\beta, \beta)

    # reverse the order of U, D, V, and shift B(0,0) to the first one, to make the sweep backward
    
    Ntau = len(phi) 
    for i in range(1,Ntau+1):
        if np.mod(i,Nstable) == 0:
            # recompute equal time GF to avoid accumulation of numerical error
            GFup_recomputed = GF(Vup[i],Dup[i],Uup[i],Uup[i-1],Dup[i-1],Vup[i-1])
            GFdn_recomputed = GF(Vdn[i],Ddn[i],Udn[i],Udn[i-1],Ddn[i-1],Vdn[i-1])

            # compare with advanced one to get the accumulated error
            GFup_error = np.max(np.abs(GFup_recomputed - GFup))
            GFdn_error = np.max(np.abs(GFdn_recomputed - GFdn))
            #print("Devitation between GF by propagating and by recomputed: ", (GFup_error + GFdn_error)/2, flush = True)
            tqdm.write(f"Devitation between GF by propagating and by recomputed: {(GFup_error + GFdn_error)/2}")
            GFup = GFup_recomputed
            GFdn = GFdn_recomputed

        # Advance equal green function first note that phi has not been updated
        Bup = ExpM(phi[i-1], ExpT,  Vdiag)
        Bdn = ExpM(phi[i-1], ExpT, -Vdiag)
        GFup = Bup.dot(GFup).dot(LA.inv(Bup))
        GFdn = Bdn.dot(GFdn).dot(LA.inv(Bdn))
        
        # Update H-S field and Green function
        GFup, GFdn, phi[i-1] = Update_HS(GFup, GFdn, phi[i-1], alpha)

        # Update SVD decompostion for B matrix
        Bup = ExpM(phi[i-1], ExpT,  Vdiag)
        Bdn = ExpM(phi[i-1], ExpT, -Vdiag)

        # Repalce SVD decomposition of B matrix
        Uup[i], Dup[i], Vup[i] = SVD_stablizer_forward(Bup, Uup[i-1], Dup[i-1], Vup[i-1])
        Udn[i], Ddn[i], Vdn[i] = SVD_stablizer_forward(Bdn, Udn[i-1], Ddn[i-1], Vdn[i-1])

    return GFup, GFdn, Uup, Dup, Vup, Udn, Ddn, Vdn, phi


def Sweep_backward(GFup,GFdn,Uup,Dup,Vup,Udn,Ddn,Vdn,phi,ExpT,Vdiag,Nstable,alpha):
    # Input GFup, GFdn should be G(0,0) = G(\beta,\beta)
    # The first member of U, D, V should be the result of I = B(0,0) or B(\beta, \beta)
    Ntau = len(phi) 
    ind_revered = np.roll(np.arange(Ntau+1)[::-1], shift=1)
    for n in range(1,Ntau+1):
        i = ind_revered[n]
        i_1 = ind_revered[n-1]
        if np.mod(n,Nstable) == 0:
            # recompute equal time GF to avoid accumulation of numerical error
            GFup_recomputed = GF(Vup[i_1],Dup[i_1],Uup[i_1],Uup[i],Dup[i],Vup[i])
            GFdn_recomputed = GF(Vdn[i_1],Ddn[i_1],Udn[i_1],Udn[i],Ddn[i],Vdn[i])
            # compare with advanced one to get the accumulated error
            GFup_error = np.max(np.abs(GFup_recomputed - GFup))
            GFdn_error = np.max(np.abs(GFdn_recomputed - GFdn))
            # print("Devitation between GF by propagating and by recomputed: ", (GFup_error + GFdn_error)/2, flush = True)
            tqdm.write(f"Devitation between GF by propagating and by recomputed: {(GFup_error + GFdn_error)/2}")
            GFup = GFup_recomputed
            GFdn = GFdn_recomputed
        
        # Update H-S field and Green function
        GFup, GFdn, phi[i-1] = Update_HS(GFup, GFdn, phi[i-1], alpha)

        # Update SVD decompostion for B matrix
        Bup = ExpM(phi[i-1], ExpT,  Vdiag)
        Bdn = ExpM(phi[i-1], ExpT, -Vdiag)

        # Advance equal time GF using new H-S field
        GFup = LA.inv(Bup).dot(GFup).dot(Bup)
        GFdn = LA.inv(Bdn).dot(GFdn).dot(Bdn)

        # Repalce SVD decomposition of B matrix
        Uup[i], Dup[i], Vup[i] = SVD_stablizer_back(Bup, Uup[i_1], Dup[i_1], Vup[i_1])
        Udn[i], Ddn[i], Vdn[i] = SVD_stablizer_back(Bdn, Udn[i_1], Ddn[i_1], Vdn[i_1])

    return GFup, GFdn, Uup, Dup, Vup, Udn, Ddn, Vdn, phi
