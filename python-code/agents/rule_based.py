class RuleBasedAgent:
    def __init__(self):
        pass

    def select_action(self, state):
        # Example policy: intervene if pest population > threshold
        pest_population = state[0]
        if pest_population > 50:
            return 1.0  # intervene strongly
        else:
            return 0.0  # no intervention
