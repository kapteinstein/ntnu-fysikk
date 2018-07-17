import sys
import numpy as np
import h5py as h5

import matplotlib
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# Get nicer looking plots than default
plt.style.use('bmh')

import seaborn as sns
sns.set(style='ticks', palette='Set1')


def options_handeler(args):
    """ Sort the run arguments """

    if ('-3d' or '-3D') in args:
        plt3d = True
    else:
        plt3d = False

    if ('-2d' or '-2D') in args:
        plt2d = True
    else:
        plt2d = False

    if ('-f' or '--files') in args:
        index = args.index('-f' or '--files')
        files = args[index + 1::]
    elif ('--files') in args:
        index = args.index('--files')
        files = args[index + 1::]
    else:
        files = None

    return plt2d, plt3d, files

def plot2d(files):
    print('\nrender 2d plot')
    plt.figure(figsize = (12,8))
    linewidth = 2
    alpha = 0.9


    for filename in files:
        try:
            f = h5.File(filename, 'r')
            X = f.get('results')
        except Exception as e:
            raise

        print(filename)
        plt.plot(X[:,0], X[:,1], lw = linewidth, alpha = alpha, label = filename)


    plt.axis([-6, 6, -4, 4])
    plt.title('Figur')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.legend(loc = 'best')
    plt.savefig('bildemothafokka_2d.pdf', filetype = 'pdf')
    plt.show()

def plot3d(files):
    print('\nrender 3d plot')

    linewidth = 2
    alpha = 0.9

    fig = plt.figure(figsize = (12,8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim3d(-5, 5)
    ax.set_ylim3d(-5, 5)
    # ax.set_zlim3d(0, 5)

    for filename in files:
        try:
            f = h5.File(filename, 'r')
            X = f.get('results')
        except Exception as e:
            raise

        print(filename)

        plt.plot(X[:,0], X[:,1], X[:,2], lw = linewidth, alpha = alpha, label = filename)

    plt.title('numerical solution to an ODE')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(loc = 'best')

    plt.savefig('bildemothafokka_3d.pdf', filetype = 'pdf')
    plt.show()

def main():

    plt2d, plt3d, files = options_handeler(sys.argv)

    if plt2d:
        plot2d(files)

    if plt3d:
        plot3d(files)


if __name__ == '__main__':
    main()
