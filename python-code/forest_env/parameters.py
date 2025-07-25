# forest_env/parameters.py

import numpy as np

# Core Simulation Settings
NUM_PATCHES = 10
SIM_WEEKS = 300
DT = 1 / 52  # One week

# Patch Characteristics
K = np.full(NUM_PATCHES, 5000)
R = np.full(NUM_PATCHES, 0.06 * DT)
EPS = np.full(NUM_PATCHES, 0.3 * DT)
BETA = 0.5 * DT
D = np.array([0.02 if i < NUM_PATCHES // 2 else 0.005 for i in range(NUM_PATCHES)])

# Strategy & Behavior Parameters
CL = 6.75   # Local cost
CT = 5.0    # Transport cost
N_STRATEGISTS = 1000
SIG = 0.1 * DT
N = 300     # Number of simulation steps

# Utility bias and infestation sensitivity
n = 0.1
f = 0.1
