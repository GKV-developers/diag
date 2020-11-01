#!/usr/bin/python

filename = './little.dat'


def check_endian(filename):
    import numpy as np
    import sys

    print 'Check endian of the file : ', filename

    print 'System byte order : ', sys.byteorder
    if sys.byteorder == 'little':
        bo1 = 'little'
        bo2 = 'big'
    if sys.byteorder == 'big':
        bo1 = 'big'
        bo2 = 'little'

    f = open(filename, mode='rb')

    f.seek(0)
    print 'Test endian : ', bo1
    for i in range(10):
        x = np.fromfile(f, dtype=np.float32, count=1)
        print x

    f.seek(0)
    print 'Test endian : ', bo2
    for i in range(10):
        x = np.fromfile(f, dtype=np.float32, count=1)
        print x.byteswap()

    f.close()


if __name__ == '__main__':
    check_endian(filename)



