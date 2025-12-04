"""
Function to read in target dataset and plot many rdf files
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_many_rdf(how_many):
    """
    Plot multiple RDF files in 3D comparison.
    
    Args:
        how_many (int): Number of rdf files to plot (rdf1.xml, rdf2.xml, etc.)
    """
    
    # Read in target dataset
    tree_target = ET.parse('xml_copy_of_robert_Ar.xml')
    root_target = tree_target.getroot()
    
    elements_target = root_target.findall('g-of-r')
    
    g_target = np.zeros(len(elements_target))
    r_target = np.zeros(len(elements_target))
    
    for i, elem in enumerate(elements_target):
        r_target[i] = float(elem.get('r'))
        g_target[i] = float(elem.get('g'))
    
    # Create 3D plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    up_to = int(10.0 / 0.1)  # r-max element of rdf-fom / bin-length from rdf data file
    
    for ii in range(1, how_many + 1):
        tree = ET.parse(f'rdf{ii}.xml')
        root = tree.getroot()
        
        elements = root.findall('g-of-r')
        
        if ii == 1:
            g = np.zeros(len(elements))
            r = np.zeros(len(elements))
        
        for i, elem in enumerate(elements):
            r[i] = float(elem.get('r'))
            g[i] = float(elem.get('g'))
        
        # Plot comparison
        index = ii  # Use ii as index for 3D visualization
        if len(r) > len(r_target):
            ax.plot(np.full(len(r_target), index), r[:len(r_target)], g[:len(r_target)], 'b-')
            ax.plot(np.full(len(r_target), index), r_target, g_target, 'r-')
        else:
            ax.plot(np.full(len(r), index), r, g, 'b-')
            ax.plot(np.full(len(r), index), r_target[:len(r)], g_target[:len(r)], 'r-')
        
        print(f'File: rdf{ii}.xml')
        print('BE AWARE THAT FOR THIS TO WORK THE OBS DATA AND CAL DATA MUST')
        print('BE STORED ON THE SAME BINNING SETUP')
        
        # Calculate FOM
        FOM = np.sum((g[:up_to] - g_target[:up_to]) ** 2) / up_to
        print(f'FOM = {FOM}\n')
    
    ax.set_ylabel('r[0.1nm]')
    ax.set_zlabel('g(r)')
    ax.set_xlabel('File index')
    
    plt.show()


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        plot_many_rdf(int(sys.argv[1]))
    else:
        print("Usage: python plot_many_rdf.py <number_of_files>")
