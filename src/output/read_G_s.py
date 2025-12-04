"""
Used to plot MDMC generated normalised g_s space-time pair correlation functions.

This function reads and plots g^s from a file, assumed to be defined as
described in handwritten notes page 29 equation (29), i.e. in theory this
function should converge to 1 for large r.
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def read_G_s(filename):
    """Read and plot MDMC generated normalised g_s space-time pair correlation functions."""
    
    tree = ET.parse(filename)
    root = tree.getroot()
    
    elements = root.findall('G-s')
    
    t_temp = np.zeros(len(elements))
    r_temp = np.zeros(len(elements))
    g_temp = np.zeros(len(elements))
    
    for i, elem in enumerate(elements):
        r_temp[i] = float(elem.get('r'))
        t_temp[i] = float(elem.get('t'))
        g_temp[i] = float(elem.get('G'))
    
    top_element = root.find('G_s-space-time-pair-correlation-function')
    bin_length = float(top_element.get('bin-length'))
    
    n_bin = int(np.ceil(np.max(r_temp) / bin_length))
    
    if len(r_temp) % n_bin != 0:
        raise ValueError('n_bin calculated incorrectly')
    
    n_time = len(r_temp) // n_bin
    
    r = r_temp[:n_bin]
    t = np.zeros(n_time)
    G_s = np.zeros((n_bin, n_time))
    
    time_index = 0
    for j in range(n_time):
        bin_index = 0
        t[time_index] = t_temp[n_bin * j]
        for i in range(n_bin):
            G_s[bin_index, time_index] = g_temp[n_bin * j + i]
            bin_index += 1
        time_index += 1
    
    # Create plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), subplot_kw={"projection": "3d"})
    
    R, T = np.meshgrid(r, t)
    
    # First subplot: log scale
    ax1.plot_surface(R, T, np.log(G_s.T))
    ax1.set_xlabel('r [AA]')
    ax1.set_ylabel('t [10^-13 s]')
    ax1.set_zlabel('log(\\tilde{g}^s) (r,t)')
    ax1.set_title(top_element.get('title'))
    
    # Second subplot: filtered
    t_cut = len(t) // 15
    r_cut = len(r) // 15
    G_s_new = G_s.copy()
    cutoff = np.max(G_s_new) / 3000
    rr, cc = np.where(G_s_new > cutoff)
    G_s_new[rr, cc] = 0
    
    ax2.plot_surface(R, T, G_s_new.T)
    ax2.set_xlabel('r [AA]')
    ax2.set_ylabel('t [10^-13 s]')
    ax2.set_zlabel('\\tilde{g}^s (r,t)')
    
    plt.tight_layout()
    plt.show()
    
    return {'r': r, 't': t, 'G_s': G_s}


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        read_G_s(sys.argv[1])
