# Determinant Quantum Monte Carlo

## Algorithm

1. Initialize H-S field for all space-time point by choosing randomly from $[1,1]$.
2. Calculate propagators $B_{\phi}(\tau,0)$ for all time slice between $[\beta, 0]$, use UDV stablization technique and store all U,D,V matrices.
3. Sweep backward then forward
4. Perform measurements after each sweep
5. Iterate 3-4 until measurement operator converge

### Details of sweep

Define direction of propagating: forward $[0,1,...,N_\tau]$, backward $[N_\tau,N_\tau-1,...,1]$.

For each $B_{\phi}(\tau,0)$, $U$, $D$, $V$ is stored, 

$U = [U_0,U_1,...,U_{N_\tau}] $, $D = [D_0, D_1,..., D_{N_\tau}] $, $V = [V_0, V_1,..., V_{N_\tau}] $.

1. For each time-slice, 

    1. flip each field at space according to acceptance ratio
    2. Update H-S field, Green function if accept proposal
    3. Replace the elements of $U$, $D$, $V$ by the new one after Green function is computed 


2. Advance equal time Green function by
$$
G(\tau_2,\tau_2) = B(\tau_2,\tau_1) G(\tau_1,\tau_1) B(\tau_2,\tau_1)^{-1}.
$$

3. Every $N_{stable}$ time slice, recompute equal time Green function from propagator using SVD decomposition
$$
G(\tau,\tau) = [1 - B(\tau,0)B(\beta,\tau)]^{-1}.
$$

### Details on update of $U$, $D$, $V$
Starting from $B(\beta,\beta)$, replace the SVD decompostion of $B(\beta,0)$ by that of $I$ after comptating the equal time Green function $G(\beta,\beta)$. 

Then replace the SVD decompostion of $B(\beta-\Delta \tau,0 )$ by that of $B(\beta,\beta -\Delta \tau)$.

For both backward sweep and forward sweep, we just need to do the replacement along the direction defined above.

## Other Mathematical formula
### Propagator

### Deviation of $\exp(-V(\phi))$, $\Delta$

### Acceptance ratio

### Green function update

## Reference

