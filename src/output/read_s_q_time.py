"""
Used to plot MDMC generated s(q,t) outputs
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def read_s_q_time(filename):
    """Read and plot MDMC generated s(q,t) outputs."""
    
    tree = ET.parse(filename)
    root = tree.getroot()
    
    elements = root.findall('SQt')
    
    q_temp = np.zeros(len(elements))
    t_temp = np.zeros(len(elements))
    Ss_temp = np.zeros(len(elements))
    Sd_temp = np.zeros(len(elements))
    
    n_q = 0
    last_q = -1000  # q should never be negative
    
    for i, elem in enumerate(elements):
        q_temp[i] = float(elem.get('q'))
        t_temp[i] = float(elem.get('t'))
        Ss_temp[i] = float(elem.get('S-self'))
        Sd_temp[i] = float(elem.get('S-diff'))
        if q_temp[i] > last_q:
            n_q += 1
            last_q = q_temp[i]
    
    top_element = root.find('s-q-time')
    
    n_time = len(q_temp) // n_q
    
    q = q_temp[:n_q]
    t = np.zeros(n_time)
    S_d = np.zeros((n_q, n_time))
    S_s = np.zeros((n_q, n_time))
    
    time_index = 0
    for j in range(n_time):
        q_index = 0
        t[time_index] = t_temp[n_q * j]
        for i in range(n_q):
            S_d[q_index, time_index] = Sd_temp[n_q * j + i]
            S_s[q_index, time_index] = Ss_temp[n_q * j + i]
            q_index += 1
        time_index += 1
    
    # Create plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), subplot_kw={"projection": "3d"})
    
    Q, T = np.meshgrid(q, t)
    
    # First subplot
    ax1.plot_surface(Q, T, S_d.T)
    ax1.set_xlabel('q [AA^-1]')
    ax1.set_ylabel('t [10^-13 s]')
    ax1.set_zlabel('S_d(q,t)')
    if top_element is not None and top_element.get('title') is not None:
        ax1.set_title(top_element.get('title'))
    
    # Second subplot
    ax2.plot_surface(Q, T, S_s.T)
    ax2.set_xlabel('q [AA^-1]')
    ax2.set_ylabel('t [10^-13 s]')
    ax2.set_zlabel('S_s(q,t)')
    
    plt.tight_layout()
    plt.show()
    
    s_q_time_info = {
        'S_s': S_s,
        'S_d': S_d,
        'q': q,
        't': t
    }
    
    return s_q_time_info


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        read_s_q_time(sys.argv[1])
