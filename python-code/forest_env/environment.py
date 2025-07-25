import numpy as np
from forest_env.parameters import *
from forest_env.model import update_state

class ForestEnv:
    def __init__(self):
        self.state = None
        self.time = 0
        self.max_time = 100  # max simulation steps

    def reset(self):
        # Initialize state (e.g., starting pest population)
        self.state = np.array([10.0])  # example initial pest population
        self.time = 0
        return self.state

    def step(self, action):
        """
        Take an action and update the environment state.
        Args:
            action (float): intervention level by agent
        Returns:
            next_state, reward, done, info
        """
        next_state = update_state(self.state, action)
        self.state = next_state
        self.time += 1

        # Define reward: e.g., negative pest population or cost of intervention
        reward = -self.state[0] - 0.1 * action  # example

        done = self.time >= self.max_time
        info = {}

        return self.state, reward, done, info

    def render(self):
        print(f"Time: {self.time}, Pest population: {self.state[0]:.2f}")
