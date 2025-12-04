"""
Takes a results xml output and plots PE epsilon vs sigma etc of the
accepted move to get some insight into correlations etc
"""

import numpy as np
import matplotlib.pyplot as plt
from analyseMDMCrun import analyseMDMCrun


def analyseMDMCcorrelation(filename):
    """
    Analyze correlations in MDMC run results.
    """
    
    result = analyseMDMCrun(filename)
    acceptX = result['acceptX']
    acceptFOM = result['acceptFOM']
    acceptSigma = result['acceptSigma']
    acceptEpsilon = result['acceptEpsilon']
    
    start = 20
    
    fig, axes = plt.subplots(3, 2, figsize=(12, 10))
    
    # Subplot 1: sigma vs epsilon
    axes[0, 0].plot(acceptSigma[start:], acceptEpsilon[start:], 'x')
    axes[0, 0].set_xlabel(r'$\sigma$')
    axes[0, 0].set_ylabel(r'$\epsilon$')
    axes[0, 0].set_title(f"starting from the {start} accepted move")
    
    # Subplot 2: Nmc vs sigma
    axes[0, 1].plot(acceptX[start:], acceptSigma[start:], 'x')
    axes[0, 1].set_xlabel('Nmc')
    axes[0, 1].set_ylabel(r'$\sigma$')
    
    # Subplot 3: epsilon vs FOM
    axes[1, 0].plot(acceptEpsilon[start:], acceptFOM[start:], 'x')
    axes[1, 0].set_xlabel(r'$\epsilon$')
    axes[1, 0].set_ylabel('FOM')
    
    # Subplot 4: sigma vs FOM
    axes[1, 1].plot(acceptSigma[start:], acceptFOM[start:], 'x')
    axes[1, 1].set_xlabel(r'$\sigma$')
    axes[1, 1].set_ylabel('FOM')
    
    # Subplot 5: Nmc vs epsilon
    axes[2, 0].plot(acceptX[start:], acceptEpsilon[start:], 'x')
    axes[2, 0].set_xlabel('Nmc')
    axes[2, 0].set_ylabel(r'$\epsilon$')
    
    # Subplot 6: Nmc vs FOM
    axes[2, 1].plot(acceptX[start:], acceptFOM[start:], 'x')
    axes[2, 1].set_xlabel('Nmc')
    axes[2, 1].set_ylabel('FOM')
    
    plt.tight_layout()
    plt.show()
    
    # Find minimum FOM
    Y = np.min(acceptFOM)
    I = np.argmin(acceptFOM)
    minAtNmc = acceptX[I]
    
    print(f'Minimum FOM: {Y} at Nmc: {minAtNmc}')
    
    # Calculate correlation coefficient
    covariance = np.column_stack([acceptFOM[start:], acceptSigma[start:], acceptEpsilon[start:]])
    corrcoef = np.corrcoef(covariance.T)
    
    print('\nCorrelation coefficient matrix:')
    print(corrcoef)
    
    return {
        'acceptX': acceptX,
        'acceptFOM': acceptFOM,
        'acceptSigma': acceptSigma,
        'acceptEpsilon': acceptEpsilon,
        'minFOM': Y,
        'minAtNmc': minAtNmc,
        'corrcoef': corrcoef
    }


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        analyseMDMCcorrelation(sys.argv[1])
