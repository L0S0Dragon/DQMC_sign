# Determinant Quantum Monte Carlo

## Algorithm

1. Initialize H-S field for all space-time point by choosing randomly from $[-1,1]$.
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
PS. propagator of $B(\tau_2,\tau_1)$ is calculated using the new configuration of H-S field.

3. Every $N_{stable}$ time slice, recompute equal time Green function from propagator using SVD decomposition

$$
G(\tau,\tau) = [I - B(\tau,0)B(\beta,\tau)]^{-1}.
$$

### Details on update of $U$, $D$, $V$
Starting from $B(\beta,\beta)$, replace the SVD decompostion of $B(\beta,0)$ by that of $I$ after comptating the equal time Green function $G(\beta,\beta)$. 

Then replace the SVD decompostion of $B(\beta-\Delta \tau,0 )$ by that of $B(\beta,\beta -\Delta \tau)$.

For both backward sweep and forward sweep, we just need to do the replacement along the direction defined above.

## Measurement
Total number of electrons to determine the doping of electron.
$$
\hat{N} = \sum_i c_i^\dagger c_i.
$$

It is equal to $\text{Tr}(G(\tau,\tau))$

Sign of Monte Carlo weight 
$$
sign[\det[I-B(\beta,0)]]
$$

In pratical, we collect quantity

$$
sign[\det[G(\beta,\beta)^{-1}]] = sign[\det[G(\beta,\beta)]] \\ = sign[\det[G(0,0)]]
$$
## Other Mathematical formula
### Propagator
$$
B_{\phi}(\tau,\tau-d\tau) = e^{-d\tau T} e^{V(\phi(\tau))}
$$

### Deviation of $\exp(-V(\phi))$, $\Delta$
Only flip one site each time

$$\Delta(\tau)_{ij} = \delta_{i,n}\delta_{j,n} \Delta(\tau)_{nn}$$

Furthermore, we have

$$
\Delta(\tau)_{ii} = e^{-\alpha\sigma (\phi(i)' - \phi(i)) } - 1  \\
=  e^{2\alpha\sigma \phi(i) } - 1
$$

with $\phi(i)' = -\phi(i)$

### Acceptance ratio
$$R = \det[I + \Delta(i,\tau) (I - G(\tau,\tau))] \\
 = 1+ \Delta(\tau)_{ii}(1 - G(\tau,\tau)_{ii})$$

### Green function update
$$
G_{jk} \rightarrow G_{jk} - \frac{1}{R} G_{ji} \Delta(\tau)_{ii} (I - G)_{ik}
$$

## Compute Green function by SVD
$$
G(\tau,\tau) = [I-B_{\phi}(\tau,0) B_{\phi}(\beta,\tau)]^{-1}
 = (VU_L)^{-1} D^{-1} (U_RU)^{-1}.
$$

$U$,$D$,$V$ are the decompostion of 
$$
(U_LU_R)^{-1} + D_RV_RV_LD_L.
$$

$$
B_{\phi}(\tau,0) = U_R D_R V_R
$$
and

$$
B_{\phi}(\beta,\tau) = V_L D_L U_L
$$


## Spinful case 
Parts for spin up and spin down are decoupled. $B_\phi$, $G_{\phi}$ are both direct product of spin up and spin down.

The acceptance ratio should be product of contributions of spin up and spin down as 
$$
R = R^{\uparrow} R^{\downarrow}
$$


## Reference

Chapter 7 in [Quantum Monte Carlo Methods](https://doi.org/10.1017/CBO9780511902581)

[A brief introduction of DQMC study in itinerant quantum critical point](http://ziyangmeng.iphy.ac.cn/files/teaching/RECS201809.pdf)


