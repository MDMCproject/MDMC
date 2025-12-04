"""
Used to plot MDMC generated rdf functions
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt


def read_rdf(filename, format_style='b-'):
    """Read and plot MDMC generated rdf functions."""
    
    tree = ET.parse(filename)
    root = tree.getroot()
    
    elements = root.findall('g-of-r')
    
    r = np.zeros(len(elements))
    g = np.zeros(len(elements))
    
    for i, elem in enumerate(elements):
        r[i] = float(elem.get('r'))
        g[i] = float(elem.get('g'))
    
    plt.plot(r, g, format_style)
    plt.xlabel('r[0.1nm]')
    plt.ylabel('g(r)')
    plt.show()
    
    return {'r': r, 'g': g}


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        format_style = sys.argv[2] if len(sys.argv) > 2 else 'b-'
        read_rdf(sys.argv[1], format_style)
