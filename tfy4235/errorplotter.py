import sys
import numpy as np
import h5py as h5

import matplotlib
from matplotlib import pyplot as plt
# Get nicer looking plots than default
plt.style.use('bmh')

import seaborn as sns
sns.set(style='ticks', palette='Set1')


def options_handeler(args):
    """ Sort the run arguments """
    return sys.argv[1:]


def plot2d(files):
    print('\nrender 2d plot')
    plt.figure(figsize = (12,8))
    plt.yscale('log')
    plt.xscale('log')
    # plt.axis([10e-6, 0, 10e-20, 0])


    for filename in files:
        try:
            f = h5.File(filename, 'r')
            X = f['error'][()]
            feil_euler = X[0]
            feil_mid = X[1]
            feil_RK4 = X[2]
            h = X[3]

        except Exception as e:
            raise


        plt.scatter(h, feil_mid, color = 'blue', marker = 'o')
        plt.scatter(h, feil_euler, color = 'red', marker = 'd')
        plt.scatter(h, feil_RK4, color = 'green', marker = 's')


    #plt.axis([-6, 6, -4, 4])
    plt.title('Figur')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig('numerical_error.pdf', filetype = 'pdf')


def main():

    files = options_handeler(sys.argv)

    plot2d(files)



if __name__ == '__main__':
    main()
