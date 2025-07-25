# main.py

import numpy as np
from forest_env.model import ForestPestModel
from forest_env.utils import plot_time_series
from forest_env import parameters as p

def load_agent(agent_type="rule"):
    if agent_type == "rule":
        from agents.rule_based import RuleBasedAgent
        return RuleBasedAgent()
    elif agent_type == "q":
        from agents.q_learning import QLearningAgent
        return QLearningAgent(state_size=100, action_size=2)
    elif agent_type == "pg":
        from agents.policy_gradient import PolicyGradientAgent
        return PolicyGradientAgent(state_size=100, action_size=2)
    elif agent_type == "evo":
        from agents.evolutionary_rl import EvolutionaryRLAgent
        return EvolutionaryRLAgent()
    else:
        raise ValueError(f"Unknown agent type: {agent_type}")


def main():
    agent_type = "rule"  # Change to "rule", "q", "pg", or "evo"
    agent = load_agent(agent_type)

    model = ForestPestModel()
    time = np.arange(p.N) * p.DT * 52  # convert to weeks

    for t in range(p.N - 1):
        for i in range(p.NUM_PATCHES):
            state = [model.I[t, i]]  # Simple state: infestation level
            action = agent.select_action(state)

            # Override strategist decisions based on agent action
            if action == 1:
                model.L[t, i] = 1.0
                model.T[t, i] = 0.0
            else:
                model.L[t, i] = 0.0
                model.T[t, i] = 1.0

        # Run one time step
        model.step(t)

        # Agent learns (if applicable)
        if hasattr(agent, "update"):
            for i in range(p.NUM_PATCHES):
                reward = -model.I[t + 1, i]  # Reward = low infestation
                next_state = [model.I[t + 1, i]]
                agent.update(state, action, reward, next_state)

    # Plot final results
    plot_time_series(time, model.S, "Susceptible Trees Over Time", "Number of Trees")
    plot_time_series(time, model.I, "Infested Trees Over Time", "Number of Trees")
    plot_time_series(time, model.L, "Proportion of Local Strategists", "Proportion")
    plot_time_series(time, model.T, "Proportion of Transport Strategists", "Proportion")


if __name__ == "__main__":
    main()
