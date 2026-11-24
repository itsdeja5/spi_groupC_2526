#!/usr/bin/env python3
"""Histogram of the peak-to-peak ripple of all the decap configurations (reads all_vpp.txt)."""
import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('all_vpp.txt', skip_header=1, dtype=str)
vpp = data[:, 0].astype(float) * 1e3                    # [mV]
cfg = data[:, 1:]

best = vpp.argmin()
uni10 = np.where((cfg == '10nF').all(axis=1))[0][0]     # five 10 nF capacitors
print(f'{len(vpp)} configurations: min {vpp.min():.1f} mV, median {np.median(vpp):.1f} mV, max {vpp.max():.1f} mV')
print('best :', ' '.join(cfg[best]), f'{vpp[best]:.2f} mV')

fig, ax = plt.subplots(figsize=(7, 3.6))
ax.hist(vpp, bins=80, color='#4477aa')
ax.axvline(vpp[best],  color='#cc3311', ls='--', label=f'best: {vpp[best]:.1f} mV')
ax.axvline(vpp[uni10], color='#228833', ls='--', label=f'5 x 10 nF: {vpp[uni10]:.1f} mV')
ax.set_xlabel('Peak-to-peak ripple at the die [mV]')
ax.set_ylabel('Number of configurations')
ax.grid(alpha=0.3)
ax.legend()
fig.tight_layout()
fig.savefig('ripple_hist.png', dpi=200)