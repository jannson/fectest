
import numpy as np
from gf2 import GF2Element, GF2DiscreteLog, GF2Poly, GF2QuotientRing
from gf2 import _bitsOf, _gf2GaussJordan

if __name__ == '__main__':
    A = np.matrix([[1,0,1,0,1],[0,1,1,0,1],[0,1,0,0,1],[1,0,0,0,0],[0,0,1,1,0]], int)
    b = _bitsOf(0b10010)
    x = _gf2GaussJordan(A, b)
    print b
    print x
    print ((A*np.matrix(x))&1).transpose()
