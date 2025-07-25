import numpy as np
import random

class QLearningAgent:
    def __init__(self, state_size, action_size, alpha=0.1, gamma=0.95, epsilon=0.1):
        self.q_table = np.zeros((state_size, action_size))
        self.alpha = alpha  # learning rate
        self.gamma = gamma  # discount factor
        self.epsilon = epsilon  # exploration rate
        self.state_size = state_size
        self.action_size = action_size

    def select_action(self, state):
        state_idx = int(state)
        if random.uniform(0, 1) < self.epsilon:
            return random.randint(0, self.action_size - 1)
        return np.argmax(self.q_table[state_idx])

    def update(self, state, action, reward, next_state):
        s = int(state)
        ns = int(next_state)
        best_next_action = np.max(self.q_table[ns])
        td_target = reward + self.gamma * best_next_action
        td_error = td_target - self.q_table[s, action]
        self.q_table[s, action] += self.alpha * td_error