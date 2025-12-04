"""
Used to plot the output of the MDMC subroutine 
print_einstein_diffuse_exp in time_corr_hist_container.f90
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt


def read_einstein_diffuse_exp(filename, format_style='b-'):
    """Read and plot Einstein diffuse experimental data."""
    
    tree = ET.parse(filename)
    root = tree.getroot()
    
    elements = root.findall('D-of-t')
    
    t = np.zeros(len(elements))
    D = np.zeros(len(elements))
    
    for i, elem in enumerate(elements):
        t[i] = float(elem.get('t'))
        D[i] = float(elem.get('D'))
    
    top_element = root.find('einstein-diffuse-exp')
    n_atom = int(top_element.get('n-atom'))
    density = float(top_element.get('density'))
    title = top_element.get('title')
    
    # Create subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    # First subplot
    ax1.plot(t, D, format_style)
    ax1.set_xlabel('t[10^-13 s]')
    ax1.set_ylabel('Einstein diffuse func')
    ax1.set_title(title)
    
    # Second subplot
    # To convert from diffusion constant as defined in code and on pages
    # 19b and 19bb in handwritten notes then the mean value of squared 
    # differences between atoms at t=0 and the atoms at some later t is
    meanSquareDist = D * t * 6
    
    ax2.plot(t, meanSquareDist, format_style)
    ax2.set_xlabel('t[10^-13 s]')
    ax2.set_ylabel('<(r(t)-r(0))^2> [(0.1nm)^2]')
    box_length = (n_atom / density) ** (1/3)
    ax2.set_title(f'Natom={n_atom}  box-length={box_length:.4f}')
    
    plt.tight_layout()
    plt.show()
    
    return {'t': t, 'D': D, 'meanSquareDist': meanSquareDist}


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        format_style = sys.argv[2] if len(sys.argv) > 2 else 'b-'
        read_einstein_diffuse_exp(sys.argv[1], format_style)
