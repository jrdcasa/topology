import numpy as np
import matplotlib.pyplot as plt
from topology.internal_coordinates import dihedral_angle_purepython, dihedral_angle_purepython_vmd


def periodic_improper(x, k, n, x0):

    # From degrees to radians
    xrad = np.deg2rad(x)
    x0rad = np.deg2rad(x0)

    return k*(1+np.cos(n*xrad-x0rad))


def harmonic_improper(x, k, x0):

    # From degrees to radians
    xrad = np.deg2rad(x)
    x0rad = np.deg2rad(x0)

    return 0.5*k*(xrad-x0rad)**2


if __name__ == "__main__":

    #                   (c4)
    #                    ||
    # HEAD--(c5)--(c3)--(c2)--(c1)--TAIL
    #
    # The line is:
    #       improper: c1 c2 c3 c4 (Type 4) --> Periodic
    #       improper: c2 c1 c3 c4 (Type 2) --> Harmonic
    c1 = np.array([1.749, -5.485, 1.191])
    c2 = np.array([0.981, -4.161, 1.349])
    c3 = np.array([-0.486, -4.175, 1.709])
    c4 = np.array([1.695, -2.846, 1.145])
    c5 = np.array([-0.994,-2.757, 1.8020])

    vmd_d1vmd = 179.951111
    vmd_d1impr = 0.028200
    #           (l)
    #           ||
    #     (k)--(j)--(i)  Type 4 (GROMACS, periodic improper)
    d1vmd1 = dihedral_angle_purepython_vmd(c1, c2, c3, c4, radians=False)
    d1vmd2 = dihedral_angle_purepython_vmd(c3, c2, c1, c4, radians=False)

    #           (l)
    #           ||
    #     (j)--(i)--(k) Type 2 (GROMACS, harmonic improper)
    d1impr1 = dihedral_angle_purepython_vmd(c2, c3, c1, c4, radians=False)
    d1impr2 = dihedral_angle_purepython_vmd(c2, c1, c3, c4, radians=False)

    print("#           (l)")
    print("#           ||")
    print("#     (k)--(j)--(i)  Type 4 (GROMACS, periodic improper)")
    print("PERIODIC IMPROPER (Type 4) => i-j-k-l order: {0:8.2f}, k-j-i-l order: {1:8.2f}".format(d1vmd1, d1vmd2))
    print("#           (l)")
    print("#           ||")
    print("#     (j)--(i)--(k) Type 2 (GROMACS, harmonic improper)")
    print("HARMONIC IMPROPER (Type 2) => i-j-k-l order: {0:8.2f}, k-j-i-l order: {1:8.2f}".format(d1impr1, d1impr2))

    # DIHEDRAL:
    print("#     (i)--(j)--(k)--(l) (GROMACS, proper dihedral)")
    dih1 = dihedral_angle_purepython_vmd(c5, c3, c2, c1, radians=False)
    dih2 = dihedral_angle_purepython_vmd(c1, c2, c3, c5, radians=False)
    print("PROPER DIHEDRAL => i-j-k-l order: {0:8.2f}, k-j-i-l order: {1:8.2f}".format(dih1, dih2))

    x_values = np.linspace(-180, 180, 360)
    y_values1 = periodic_improper(x_values, 43.932, 2, 180.0)
    y_values2 = harmonic_improper(x_values, np.pi*43.932, 0.0)

    plt.plot(x_values, y_values1, label='Periodic improper')
    plt.plot(x_values, y_values2, label='Harmonic improper')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.xlim(-80, 80)
    plt.ylim(0, 90)
    plt.grid(True)
    plt.legend()
    plt.show()