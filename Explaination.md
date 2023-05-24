# Determinant Quantum Monte Carlo

## Algorithm

1. Initialize H-S field for all space-time point by choosing randomly from $[1,1]$.
2. Calculate propagators $B_{\phi}(\tau,0)$ for all time slice between $[\beta, 0]$, use UDV stablization technique and store all U,D,V matrices.
3. Sweep backward then forward
4. Perform measurements after each sweep
5. Iterate 3-4 until measurement operator converge

### Details of sweep

Define propagating direction forward $[0,1,...,N_\tau]$, backward $[N_\tau,N_\tau-1,...,1]$.

For each $B_{\phi}(\tau,0)$, $U$, $D$, $V$ is stored, 

$U = [U_0,U_1,...,U_{N_\tau}] $, $D = [D_0, D_1,..., D_{N_\tau}] $, $V = [V_0, V_1,..., V_{N_\tau}] $.

1. For each time-slice, 

    1. flip each field at space according to acceptance ratio
    2. Update H-S field, Green function if accept proposal
    3. Replace the elements of $U$, $D$, $V$ by the new one after Green function is computed 

2. Advance 