import math as m
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Read data from text files
def read_data(file_path):
    with open(file_path, 'r') as file:
        data = np.loadtxt(file)
    return data

# File paths for your text files
file_paths = ['data1.txt', 'data2.txt', 'data3.txt', 'data4.txt']

# Read data from files
data = [read_data(file_path) for file_path in file_paths]

# Create subplots
fig, axs = plt.subplots(2, 2, figsize=(10, 8))
axs = axs.flatten()

# Plot data from text files
for i in range(len(file_paths)):
    x = np.arange(1, len(data[i]) + 1)  # Assuming data is 1-dimensional
    axs[i].plot(x, data[i], marker='o')
    axs[i].set_title(f'Plot {i + 1}')
    axs[i].set_xlabel('X-Axis Label')
    axs[i].set_ylabel('Y-Axis Label')

# Adjust layout and display plots
plt.tight_layout()
plt.show()
