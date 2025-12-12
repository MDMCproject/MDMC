"""
Used to plot MDMC generated g_d space-time pair correlation functions.

As of this writing MDMC outputs the "normalised" version of this function
which has the property that (g^norm)_d(r,t) -> 1 when t->infinity and
r->infinity. Note (g^norm)_d(r,t) = N * g_d(r,t) / (N-1).

Note the space-time correlation function is related to the space-time
pair correlation function by G_d = rho * g_d, where rho is N/V.
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def read_G_d(filename):
    """Read and plot MDMC generated G_d space-time pair correlation functions."""
    
    tree = ET.parse(filename)
    root = tree.getroot()
    
    elements = root.findall('G-d')
    
    t_temp = np.zeros(len(elements))
    r_temp = np.zeros(len(elements))
    g_temp = np.zeros(len(elements))
    
    for i, elem in enumerate(elements):
        r_temp[i] = float(elem.get('r'))
        t_temp[i] = float(elem.get('t'))
        g_temp[i] = float(elem.get('G'))
    
    bin_length = float(root.get('bin-length'))    
    
    n_bin = int(np.ceil(np.max(r_temp) / bin_length))
    
    if len(r_temp) % n_bin != 0:
        raise ValueError('n_bin calculated incorrectly')
    
    n_time = len(r_temp) // n_bin
    
    r = r_temp[:n_bin]
    t = np.zeros(n_time)
    G_d = np.zeros((n_bin, n_time))
    
    time_index = 0
    for j in range(n_time):
        bin_index = 0
        t[time_index] = t_temp[n_bin * j]
        for i in range(n_bin):
            G_d[bin_index, time_index] = g_temp[n_bin * j + i]
            bin_index += 1
        time_index += 1
    
    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    R, T = np.meshgrid(r, t)
    ax.plot_surface(R, T, G_d.T)
    ax.set_xlabel('r [AA]')
    ax.set_ylabel('t [10^-13 s]')
    ax.set_zlabel('"normalised" g_d (r,t)')
    ax.set_title(root.get('title'))
    
    plt.show()
    
    return {'r': r, 't': t, 'G_d': G_d}


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        read_G_d(sys.argv[1])
    else:
        print("Usage: python read_G_d.py <mdmc results xml file>")
        # for debugging purposes
        print("Try reading from hardcoded file")
        read_G_d('/root/mdmc/build/bin/output/first_g_d.xml')
