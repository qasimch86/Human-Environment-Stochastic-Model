# forest_env/utils.py

import matplotlib.pyplot as plt
import numpy as np

def plot_time_series(time, data, title, ylabel, patch_labels=None):
    plt.figure(figsize=(10, 5))
    for i in range(data.shape[1]):
        label = patch_labels[i] if patch_labels else f"Patch {i+1}"
        plt.plot(time, data[:, i], label=label)
    plt.xlabel("Time (weeks)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_bar(x, height, title, xlabel, ylabel):
    plt.figure(figsize=(8, 4))
    plt.bar(x, height)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
