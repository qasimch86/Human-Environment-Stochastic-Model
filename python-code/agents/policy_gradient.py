import numpy as np

class PolicyGradientAgent:
    def __init__(self, state_size, action_size, learning_rate=0.01):
        self.state_size = state_size
        self.action_size = action_size
        self.policy = np.random.rand(state_size, action_size)
        self.lr = learning_rate

    def softmax(self, logits):
        exp_logits = np.exp(logits - np.max(logits))
        return exp_logits / np.sum(exp_logits)

    def select_action(self, state):
        probs = self.softmax(self.policy[int(state)])
        return np.random.choice(len(probs), p=probs)

    def update(self, state, action, reward, next_state):
        s = int(state)
        probs = self.softmax(self.policy[s])
        grad = -probs
        grad[action] += 1
        self.policy[s] += self.lr * reward * grad