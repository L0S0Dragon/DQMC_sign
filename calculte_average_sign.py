from run_DQMC import run
import numpy as np
from itertools import product
# Check if we are in Jupyter
if 'ipykernel' in sys.modules:
    from tqdm.notebook import tqdm
else:
    from tqdm import tqdm

N = 4
dtau = 0.1
Nstable = 10
t = -1.0
Nwarm = 500
Nsweep = 20000

chemical_potential = np.linspace(-2,0,10)
imaginary_time_slice = np.linspace(20,50,5,dtype = int)
onsite_interaction = np.linspace(2,8,7)
params = list(product(chemical_potential, imaginary_time_slice, onsite_interaction))

# N params, 2 quantities, Nsweep data
result = np.zeros((len(params),2,Nsweep))
for i,p in tqdm(enumerate(params)):
    mu , Ntau, U = p
    # 0 = sign, 1 = electron number
    result[i,0], result[i,1] = run(N, Ntau, dtau, Nstable, U, t, mu, Nwarm, Nsweep)

np.save("params", params)
np.save("result", result)