"""
Takes a results xml output and plots the accepted/rejected FOM values
and values of the PE parameters
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt


def analyseMDMCrun(filename):
    """
    Read and analyze MDMC run results from XML file.
    
    Returns:
        dict: Contains acceptX, acceptFOM, acceptSigma, acceptEpsilon,
              rejectX, rejectFOM, rejectSigma, rejectEpsilon
    """
    
    tree = ET.parse(filename)
    root = tree.getroot()
    
    accept_elements = root.findall('accept')
    rejected_elements = root.findall('rejected')
    
    nAccept = len(accept_elements)
    nRejected = len(rejected_elements)
    
    acceptXaxis = np.zeros(nAccept)
    acceptYaxis = np.zeros(nAccept)
    acceptYaxisSigma = np.zeros(nAccept)
    acceptYaxisEpsilon = np.zeros(nAccept)
    
    rejectedXaxis = np.zeros(nRejected)
    rejectedYaxis = np.zeros(nRejected)
    rejectedYaxisSigma = np.zeros(nRejected)
    rejectedYaxisEpsilon = np.zeros(nRejected)
    
    # Parse accept elements
    iAcc = 0
    for elem in accept_elements:
        acceptXaxis[iAcc] = iAcc + 1
        acceptYaxis[iAcc] = float(elem.get('val'))
        acceptYaxisSigma[iAcc] = float(elem.get('sigma'))
        acceptYaxisEpsilon[iAcc] = float(elem.get('epsilon'))
        iAcc += 1
    
    # Parse rejected elements
    iRej = 0
    iMC = nAccept + 1
    for elem in rejected_elements:
        rejectedXaxis[iRej] = iMC
        rejectedYaxis[iRej] = float(elem.get('val'))
        rejectedYaxisSigma[iRej] = float(elem.get('sigma'))
        rejectedYaxisEpsilon[iRej] = float(elem.get('epsilon'))
        iRej += 1
        iMC += 1
    
    # Create plots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
    
    ax1.plot(acceptXaxis, np.log10(acceptYaxis), 'bx', label='Accepted')
    ax1.plot(rejectedXaxis, np.log10(rejectedYaxis), 'ro', label='Rejected')
    percentage_accepted = 100 * nAccept / (nAccept + nRejected)
    ax1.set_title(f'percentage accepted = {percentage_accepted:.1f}%')
    ax1.set_ylabel('log10(FOM)')
    ax1.set_xlabel('MC step')
    ax1.legend()
    
    ax2.plot(acceptXaxis, acceptYaxisSigma, 'bx', label='Accepted')
    ax2.plot(rejectedXaxis, rejectedYaxisSigma, 'ro', label='Rejected')
    ax2.set_ylabel(r'$\sigma$')
    ax2.set_xlabel('MC step')
    ax2.legend()
    
    ax3.plot(acceptXaxis, acceptYaxisEpsilon, 'bx', label='Accepted')
    ax3.plot(rejectedXaxis, rejectedYaxisEpsilon, 'ro', label='Rejected')
    ax3.set_ylabel(r'$\epsilon$')
    ax3.set_xlabel('MC step')
    ax3.legend()
    
    plt.tight_layout()
    plt.show()
    
    return {
        'acceptX': acceptXaxis,
        'acceptFOM': acceptYaxis,
        'acceptSigma': acceptYaxisSigma,
        'acceptEpsilon': acceptYaxisEpsilon,
        'rejectX': rejectedXaxis,
        'rejectFOM': rejectedYaxis,
        'rejectSigma': rejectedYaxisSigma,
        'rejectEpsilon': rejectedYaxisEpsilon
    }


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        analyseMDMCrun(sys.argv[1])
