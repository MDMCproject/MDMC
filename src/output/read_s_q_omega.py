"""
Used to plot s(q,omega) XML files
The required XML format is:
<s-q-omega>
    <SQomega q="0.42000" omega="0.00000" S-self="25.50111" S-diff="-24.30186" S="1.23" error="0.1" />
    ... etc
</s-q-omega>

When the file contains S-self and S-diff then S(q,omega) = S-self + S-diff. 
Otherwise S(q,omega) = S
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def read_s_q_omega(filename):
    """Read and plot S(q,omega) XML files."""
    
    tree = ET.parse(filename)
    root = tree.getroot()
    
    elements = root.findall('SQomega')
    
    q_temp = np.zeros(len(elements))
    omega_temp = np.zeros(len(elements))
    Ss_temp = np.zeros(len(elements))
    Sd_temp = np.zeros(len(elements))
    Stot_temp = np.zeros(len(elements))
    
    n_q = 0
    last_q = -1000  # q should never be negative
    
    for i, elem in enumerate(elements):
        q_temp[i] = float(elem.get('q'))
        omega_temp[i] = float(elem.get('omega'))
        
        if elem.get('S') is not None:
            if elem.get('S') != 'no data':
                Stot_temp[i] = float(elem.get('S'))
            else:
                Stot_temp[i] = np.nan
        else:
            Ss_temp[i] = float(elem.get('S-self'))
            Sd_temp[i] = float(elem.get('S-diff'))
            Stot_temp[i] = Ss_temp[i] + Sd_temp[i]
        
        if q_temp[i] > last_q:
            n_q += 1
            last_q = q_temp[i]
    
    top_element = root.find('s-q-omega')
    
    n_time = len(q_temp) // n_q
    
    q = q_temp[:n_q]
    omega = np.zeros(n_time)
    S_tot = np.zeros((n_q, n_time))
    
    time_index = 0
    for j in range(n_time):
        q_index = 0
        omega[time_index] = omega_temp[n_q * j]
        for i in range(n_q):
            S_tot[q_index, time_index] = Stot_temp[n_q * j + i]
            q_index += 1
        time_index += 1
    
    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    Q, Omega = np.meshgrid(q, omega)
    ax.plot_surface(Q, Omega, S_tot.T)
    ax.set_xlabel('q [AA^-1]')
    ax.set_ylabel('omega 1/[10^-13 s]')
    ax.set_zlabel('S')
    if top_element is not None and top_element.get('title') is not None:
        ax.set_title(top_element.get('title'))
    
    plt.show()
    
    s_q_omega_info = {
        'S_tot': S_tot,
        'q': q,
        'omega': omega
    }
    
    return s_q_omega_info


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        read_s_q_omega(sys.argv[1])
