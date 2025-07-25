# main.py

import numpy as np
from forest_env.model import ForestPestModel
from forest_env.utils import plot_time_series
from forest_env import parameters as p

def main():
    # Initialize and run model
    model = ForestPestModel()
    model.run()

    # Time vector
    time = np.arange(p.N) * p.DT * 52  # convert to weeks

    # Plot results like Figure 2 in your paper
    plot_time_series(time, model.S, "Susceptible Trees Over Time", "Number of Trees")
    plot_time_series(time, model.I, "Infested Trees Over Time", "Number of Trees")
    plot_time_series(time, model.L, "Proportion of Local Strategists", "Proportion")
    plot_time_series(time, model.T, "Proportion of Transport Strategists", "Proportion")

if __name__ == "__main__":
    main()
