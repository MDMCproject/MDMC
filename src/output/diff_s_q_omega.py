"""
Take as input two S(q,omega) files and plot these and differences
between these, for example to compare a MDMC generated S(q,omega) file
with a data S(q,omega). The required XML format is:
<s-q-omega>
    <SQomega q="0.42000" omega="0.00000" S-self="25.50111" S-diff="-24.30186" S="1.23" error="0.1" />
    ... etc
</s-q-omega>

It is currently assumed that filename1 contains S-self and S-diff and
S_1(q,omega) = S-self + S-diff. filename2 may either contain S-self and S-diff
or S attribute. If the latter is the case then S_2(q,omega) = S
"""

import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def diff_s_q_omega(filename1, filename2):
    """
    Compare two S(q,omega) files and plot differences.
    """
    
    # Reading info from the 1st file
    tree1 = ET.parse(filename1)
    root1 = tree1.getroot()
    
    elements1 = root1.findall('SQomega')
    
    q_temp = np.zeros(len(elements1))
    omega_temp = np.zeros(len(elements1))
    Ss_temp = np.zeros(len(elements1))
    Sd_temp = np.zeros(len(elements1))
    
    n_q = 0
    last_q = -1000  # q should never be negative
    
    for i, elem in enumerate(elements1):
        q_temp[i] = float(elem.get('q'))
        omega_temp[i] = float(elem.get('omega'))
        Ss_temp[i] = float(elem.get('S-self'))
        Sd_temp[i] = float(elem.get('S-diff'))
        if q_temp[i] > last_q:
            n_q += 1
            last_q = q_temp[i]
    
    n_time = len(q_temp) // n_q
    
    q = q_temp[:n_q]
    omega = np.zeros(n_time)
    S_d = np.zeros((n_q, n_time))
    S_s = np.zeros((n_q, n_time))
    
    time_index = 0
    for j in range(n_time):
        q_index = 0
        omega[time_index] = omega_temp[n_q * j]
        for i in range(n_q):
            S_d[q_index, time_index] = Sd_temp[n_q * j + i]
            S_s[q_index, time_index] = Ss_temp[n_q * j + i]
            q_index += 1
        time_index += 1
    
    S_tot1 = S_s + S_d
    
    # Reading info from the 2nd file
    tree2 = ET.parse(filename2)
    root2 = tree2.getroot()
    
    elements2 = root2.findall('SQomega')
    
    if len(elements1) != len(elements2):
        print('Number of SQomega elements in filename1 and filename2 must be the same')
        return None
    
    q_temp2 = np.zeros(len(elements2))
    omega_temp2 = np.zeros(len(elements2))
    Ss_temp2 = np.zeros(len(elements2))
    Sd_temp2 = np.zeros(len(elements2))
    Stot_temp2 = np.zeros(len(elements2))
    error2 = np.zeros(len(elements2))
    
    # Check if errors in filename2
    errors_in_filename2 = False
    if len(elements2) > 0:
        errors_in_filename2 = elements2[0].get('error') is not None
    
    n_q = 0
    last_q = -1000
    
    for i, elem in enumerate(elements2):
        q_temp2[i] = float(elem.get('q'))
        omega_temp2[i] = float(elem.get('omega'))
        
        if elem.get('S') is not None:
            if elem.get('S') != 'no data':
                Stot_temp2[i] = float(elem.get('S'))
                if elem.get('error') is not None:
                    error2[i] = float(elem.get('error'))
            else:
                Stot_temp2[i] = np.nan
                error2[i] = np.nan
        else:
            Ss_temp2[i] = float(elem.get('S-self'))
            Sd_temp2[i] = float(elem.get('S-diff'))
            Stot_temp2[i] = Ss_temp2[i] + Sd_temp2[i]
            if errors_in_filename2:
                error2[i] = float(elem.get('error'))
        
        if q_temp2[i] > last_q:
            n_q += 1
            last_q = q_temp2[i]
    
    n_time2 = len(q_temp2) // n_q
    
    q2 = q_temp2[:n_q]
    omega2 = np.zeros(n_time2)
    S_tot2 = np.zeros((n_q, n_time2))
    error2Plot = np.zeros((n_q, n_time2))
    
    time_index = 0
    for j in range(n_time2):
        q_index = 0
        omega2[time_index] = omega_temp2[n_q * j]
        for i in range(n_q):
            S_tot2[q_index, time_index] = Stot_temp2[n_q * j + i]
            error2Plot[q_index, time_index] = error2[n_q * j + i]
            q_index += 1
        time_index += 1
    
    # Calculate FOM difference between the two input S(q, omega)
    val = 0
    val_with_error = 0
    
    for i_q in range(len(q)):
        for i_o in range(len(omega)):
            if not np.isnan(S_tot2[i_q, i_o]):
                val += (S_tot2[i_q, i_o] - S_tot1[i_q, i_o]) ** 2
                if errors_in_filename2:
                    val_with_error += (S_tot2[i_q, i_o] - S_tot1[i_q, i_o]) ** 2 / error2Plot[i_q, i_o] ** 2
    
    print(f'FOM_val_without_errors = {val}')
    if errors_in_filename2:
        print(f'FOM_val_with_errors = {val_with_error}')
    
    # Create plots
    fig = plt.figure(figsize=(14, 10))
    
    # Subplot 1
    ax1 = fig.add_subplot(3, 2, 1, projection='3d')
    Q, Omega = np.meshgrid(q, omega)
    ax1.plot_surface(Q, Omega, S_tot1.T)
    ax1.set_xlabel('q [AA^-1]')
    ax1.set_ylabel('omega 1/[10^-13 s]')
    ax1.set_zlabel('S_1(q,omega)')
    
    # Subplot 3
    ax3 = fig.add_subplot(3, 2, 3, projection='3d')
    ax3.plot_surface(Q, Omega, S_tot2.T)
    ax3.set_xlabel('q [AA^-1]')
    ax3.set_ylabel('omega 1/[10^-13 s]')
    ax3.set_zlabel('S_2(q,omega)')
    
    # Subplot 5
    ax5 = fig.add_subplot(3, 2, 5, projection='3d')
    ax5.plot_surface(Q, Omega, (S_tot2 - S_tot1).T)
    ax5.set_xlabel('q [AA^-1]')
    ax5.set_ylabel('omega 1/[10^-13 s]')
    ax5.set_zlabel('S_2 - S_1')
    
    # Subplot 2
    ax2 = fig.add_subplot(3, 2, 2)
    
    Sq1 = np.zeros(len(q))
    for i in range(len(q)):
        for j in range(len(omega) - 1):
            if not np.isnan(S_tot1[i, j]) and not np.isnan(S_tot1[i, j + 1]):
                Sq1[i] += (omega[j + 1] - omega[j]) * (S_tot1[i, j] + S_tot1[i, j + 1]) / 2
    
    ax2.plot(q, Sq1, 'b', label=r'$\int S_1(q) d\omega$')
    
    Sq2 = np.zeros(len(q))
    for i in range(len(q)):
        for j in range(len(omega) - 1):
            if not np.isnan(S_tot2[i, j]) and not np.isnan(S_tot2[i, j + 1]):
                Sq2[i] += (omega[j + 1] - omega[j]) * (S_tot2[i, j] + S_tot2[i, j + 1]) / 2
    
    ax2.plot(q, Sq2, 'r--', label=r'$\int S_2(q) d\omega$')
    ax2.set_xlabel('q [AA^-1]')
    ax2.legend()
    
    # Subplot 4
    ax4 = fig.add_subplot(3, 2, 4, projection='3d')
    if errors_in_filename2:
        ax4.plot_surface(Q, Omega, (np.abs(S_tot2 - S_tot1) / error2Plot).T)
        ax4.set_zlabel('abs(S_2-S_1)/error')
    else:
        ax4.plot_surface(Q, Omega, np.abs(S_tot2 - S_tot1).T)
        ax4.set_zlabel('abs(S_2-S_1)')
    ax4.set_xlabel('q [AA^-1]')
    ax4.set_ylabel('omega 1/[10^-13 s]')
    
    # Subplot 6
    if errors_in_filename2:
        ax6 = fig.add_subplot(3, 2, 6, projection='3d')
        ax6.plot_surface(Q, Omega, error2Plot.T)
        ax6.set_xlabel('q [AA^-1]')
        ax6.set_ylabel('omega 1/[10^-13 s]')
        ax6.set_zlabel('error on data')
    
    plt.tight_layout()
    plt.show()
    
    return {
        'S_tot1': S_tot1,
        'S_tot2': S_tot2,
        'q': q,
        'omega': omega,
        'FOM_without_errors': val,
        'FOM_with_errors': val_with_error if errors_in_filename2 else None
    }


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 2:
        diff_s_q_omega(sys.argv[1], sys.argv[2])
