"""
Used to plot MDMC generated h_d functions
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def read_h_d(filename):
    """Read and plot MDMC generated h_d functions."""
    
    tree = ET.parse(filename)
    root = tree.getroot()
    
    elements = root.findall('h-d')
    
    t_temp = np.zeros(len(elements))
    r_temp = np.zeros(len(elements))
    h_temp = np.zeros(len(elements))
    
    for i, elem in enumerate(elements):
        r_temp[i] = float(elem.get('r'))
        t_temp[i] = float(elem.get('t'))
        h_temp[i] = float(elem.get('h'))
    
    top_element = root.find('normalised-h-d-histogram')
    bin_length = float(top_element.get('bin-length'))
    
    n_bin = int(np.ceil(np.max(r_temp) / bin_length))
    
    if len(r_temp) % n_bin != 0:
        raise ValueError('n_bin calculated incorrectly')
    
    n_time = len(r_temp) // n_bin
    
    r = r_temp[:n_bin]
    t = np.zeros(n_time)
    h_d = np.zeros((n_bin, n_time))
    
    time_index = 0
    for j in range(n_time):
        bin_index = 0
        t[time_index] = t_temp[n_bin * j]
        for i in range(n_bin):
            h_d[bin_index, time_index] = h_temp[n_bin * j + i]
            bin_index += 1
        time_index += 1
    
    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    R, T = np.meshgrid(r, t)
    ax.plot_surface(R, T, h_d.T)
    ax.set_xlabel('r [AA]')
    ax.set_ylabel('t [10^-13 s]')
    ax.set_zlabel('h_d (r,t)')
    ax.set_title(top_element.get('title'))
    
    plt.show()
    
    h_d_info = {
        'r': r,
        't': t,
        'val': h_d,
        'bin_length': bin_length,
        'n_bin': n_bin
    }
    
    return h_d_info


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        read_h_d(sys.argv[1])
