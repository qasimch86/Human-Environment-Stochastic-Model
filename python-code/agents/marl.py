'''Each Patch Can be Its Own Agent (Multi-Agent RL)
If you want to treat each patch as a separate learning agent, then you are entering Multi-Agent RL (MARL). This is already scaffolded in agents/marl.py.

In that case:

Concept	Implementation
Environment	Still ForestPestModel
Agent	One per patch (agents[i])
State	Local to each patch (e.g., I[t, i])
Action	Patch-specific (e.g., strategy switch)
Reward	Patch-specific (e.g., local infestation drop)
Coordination?	Optional — agents could cooperate or compete'''

from agents.q_learning import QLearningAgent

class MultiAgentRL:
    def __init__(self, num_agents):
        self.agents = [QLearningAgent(state_size=10, action_size=2) for _ in range(num_agents)]

    def select_action(self, agent_id, state):
        return self.agents[agent_id].select_action(state)

    def update(self, agent_id, state, action, reward, next_state):
        self.agents[agent_id].update(state, action, reward, next_state)

