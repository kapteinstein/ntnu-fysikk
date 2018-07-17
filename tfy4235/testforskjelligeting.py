# import numpy as np
# from matplotlib import pyplot as plt
# import netCDF4 as nc
# import time

def datafile(N):
    ''' Use netCDF to store variables '''

    def write(N):
        ''' Write to netCDF file '''
        r = np.random.normal(size = N, loc = 0, scale = 1000)
        # This creates and opens the file random.nc in write mode.
        d = nc.Dataset('random.nc', 'w')
        # Create time-dimension, with name 't' and size len(r)
        t = d.createDimension('t', len(r))
        # The arguments to createVariable are name, datatype
        # and name of dimension(s)
        v = d.createVariable('random', r.dtype, ('t',))
        # Note the use of array indexing to fill variable with data.
        # If you simply say v = r, v will stop being of the type
        # NetCDF variable, and become a numpy array instead.
        v[:] = r
        d.close()

    def read(file):
        ''' Read from netCDF file '''
        d = nc.Dataset(file, 'r')
        r = d.variables['random']
        # r is now a NetCDF variable type, using array indexing to
        # get values. Use r[:] to get the entire array.
        print(r[0:10])
        d.close()

    write(N)
    read('random.nc')

def crossProduct():
    a = np.array([[1,2,3],[2,2,2]])
    b = np.array([2,3,4])

    cross = np.cross(a[1,:],b)
    print(cross)

def Signal(t):
    q = np.sin(t)

class Coordinates(object):
    def __init__(self, x, y):
        self.x = x
        self.y = y
    def __str__(self):
        return('yay')
    def distance(self, other):
        return ((self.x-other.x)**2 + (self.y-other.y)**2)**(0.5)


class Clock(object):
    """docstring for clock"""
    def __init__(self, time):
        self.time = time
    def print_time(self):
        time = 'asdf'
        print(self.time)

class intSet(object):
    """An intSet is a set of integers
    The value is represented by a list of ints, self.vals.
    Each int in the set occurs in self.vals exactly once."""

    def __init__(self):
        """Create an empty set of integers"""
        self.vals = []

    def insert(self, e):
        """Assumes e is an integer and inserts e into self"""
        if not e in self.vals:
            self.vals.append(e)

    def member(self, e):
        """Assumes e is an integer
           Returns True if e is in self, and False otherwise"""
        return e in self.vals

    def remove(self, e):
        """Assumes e is an integer and removes e from self
           Raises ValueError if e is not in self"""
        try:
            self.vals.remove(e)
        except:
            raise ValueError(str(e) + ' not found')

    def __str__(self):
        """Returns a string representation of self"""
        self.vals.sort()
        return '{' + ','.join([str(e) for e in self.vals]) + '}'

    def intersect(self, other):
        """Assumes other is an intSet
           Returns a new intSet containing elements that appear in both sets."""
        # Initialize a new intSet
        commonValueSet = intSet()
        # Go through the values in this set
        for val in self.vals:
            # Check if each value is a member of the other set
            if other.member(val):
                commonValueSet.insert(val)
        return commonValueSet

    def __len__(self):
        return len(self.vals)


def main():



    a = intSet()
    a.insert(1)
    a.insert(2)

    b = intSet()
    b.insert(2)
    b.insert(4)
    print(b)
    print(a)
    a.intersect(b)
    print(a)

    '''
    a = np.random.normal(size = 10000000, loc = 0, scale = 1000)

    tic = time.time()
    signal = np.vectorize(Signal, otypes = [np.float])
    signal(a)
    toc = time.time()
    print(toc-tic)

    tic = time.time()
    for i in a:
        Signal(i)
    toc = time.time()
    print(toc-tic)

    tic = time.time()
    Signal(a)
    toc = time.time()

    print(toc-tic)

    crossProduct()
    print('hei')
    # crossProduct()
    '''

    return 0


if __name__ == '__main__':
    main()
