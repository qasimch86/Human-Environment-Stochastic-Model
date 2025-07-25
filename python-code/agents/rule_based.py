class RuleBasedAgent:
    def __init__(self, threshold=50):
        self.threshold = threshold

    def select_action(self, state):
        pest_population = state[0]  # assuming state[0] = infestation
        return 1 if pest_population > self.threshold else 0

    def update(self, state, action, reward, next_state):
        pass  # No learning for rule-based agent