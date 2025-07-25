# forest_env/model.py

import numpy as np
from forest_env import parameters as p

class ForestPestModel:
    def __init__(self):
        self.num_patches = p.NUM_PATCHES
        self.time_steps = p.N
        self.dt = p.DT

        self.K = p.K.copy()
        self.S = np.zeros((self.time_steps, self.num_patches))
        self.I = np.zeros_like(self.S)
        self.R = np.zeros_like(self.S)
        self.L = np.zeros_like(self.S)
        self.T = np.zeros_like(self.S)
        self.NL = np.zeros_like(self.S)
        self.NT = np.zeros_like(self.S)

        self._initialize()

    def _initialize(self):
        for i in range(self.num_patches):
            self.S[0, i] = self.K[i] - 15 if i == 0 else self.K[i]
            self.I[0, i] = self.K[i] - self.S[0, i]
            self.NL[0, i] = 100
            self.NT[0, i] = p.N_STRATEGISTS - self.NL[0, i]
            self.L[0, i] = self.NL[0, i] / p.N_STRATEGISTS
            self.T[0, i] = self.NT[0, i] / p.N_STRATEGISTS

    def step(self, t):
        for i in range(self.num_patches):
            # Forest dynamics: infection and regeneration
            P_SI = p.BETA * self.I[t, i] / self.K[i]
            P_S = p.R[i] * (1 - (self.S[t, i] + self.I[t, i]) / self.K[i])
            P_IR = p.EPS[i]

            # Utilities
            U_L = -p.CL + p.n * (self.L[t, i] - 0.5)
            U_T = -p.CT + p.n * (0.5 - self.L[t, i]) - p.f * self.I[t, i]

            # Strategy switching probabilities
            if U_L > U_T:
                P_LT = 0
                P_TL = p.SIG * (U_L - U_T)
            else:
                P_TL = 0
                P_LT = p.SIG * (U_T - U_L)

            # Cross-patch infestation: transport strategists move pests
            P_SIT = 0
            for j in range(self.num_patches):
                if i != j:
                    P_SIT += p.BETA * p.D[j] * self.T[t, j] * self.I[t, j] / self.K[j]

            # Stochastic events
            num_SI = np.random.binomial(self.S[t, i], P_SI)
            num_S = np.random.binomial(self.S[t, i], P_S)
            num_IR = np.random.binomial(self.I[t, i], P_IR)
            num_LT = np.random.binomial(self.NL[t, i], P_LT)
            num_TL = np.random.binomial(self.NT[t, i], P_TL)
            num_SIT = np.random.binomial(self.S[t, i], P_SIT)

            # Update forest state
            self.S[t + 1, i] = self.S[t, i] - num_SI + num_S - num_SIT
            self.I[t + 1, i] = max(self.I[t, i] + num_SI - num_IR + num_SIT, 0)
            self.R[t + 1, i] = self.R[t, i] + num_IR

            # Update strategist numbers
            self.NL[t + 1, i] = self.NL[t, i] - num_LT + num_TL
            self.NT[t + 1, i] = p.N_STRATEGISTS - self.NL[t + 1, i]

            self.L[t + 1, i] = self.NL[t + 1, i] / p.N_STRATEGISTS
            self.T[t + 1, i] = self.NT[t + 1, i] / p.N_STRATEGISTS

    def run(self):
        for t in range(self.time_steps - 1):
            self.step(t)

    def get_final_state(self):
        return {
            "S": self.S[-1],
            "I": self.I[-1],
            "R": self.R[-1],
            "L": self.L[-1],
            "T": self.T[-1],
        }

