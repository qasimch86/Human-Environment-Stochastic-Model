import numpy as np

class EvolutionaryRLAgent:
    def __init__(self, strategy_space=[0, 1], mutation_rate=0.1):
        self.strategy = np.random.choice(strategy_space)
        self.mutation_rate = mutation_rate
        self.reward_history = []

    def select_action(self, state):
        return self.strategy

    def update(self, state, action, reward, next_state):
        self.reward_history.append(reward)
        if len(self.reward_history) >= 10:
            avg_reward = np.mean(self.reward_history[-10:])
            if np.random.rand() < self.mutation_rate:
                self.strategy = 1 - self.strategy  # flip strategy