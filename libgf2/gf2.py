'''
Created on Oct 5, 2013

@author: jmsachs

Copyright 2013-2015 Jason M Sachs

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

'''

import numpy as np
try:
    # Speedup if Cython is installed.
    import pyximport; pyximport.install()
    import fastgf2mulmod
except:
    fastgf2mulmod = None

class _GF2(object):
    @staticmethod
    def mod(x,m):
        m2 = m
        i = 0
        while m2 < x:
            m2 <<= 1
            i += 1
        while i >= 0:
            xnew = x ^ m2
            if xnew < x:
                x = xnew
            m2 >>= 1
            i -= 1
        return x
    @staticmethod
    def mul(x,y):
        z = 0
        while x > 0:
            if (x & 1) != 0:
                z ^= y
            y <<= 1
            x >>= 1
        return z
    @staticmethod
    def bitlength(x):
        n = 0
        while x > 0:
            n += 1
            x >>= 1
        return n
    @staticmethod
    def mulmod(x,y,m):
        z = 0
        while x > 0:
            if (x & 1) != 0:
                z ^= y
            y <<= 1
            y2 = y ^ m
            if y2 < y:
                y = y2
            x >>= 1
        return z
    @staticmethod
    def power(x,k):
        z = 0
        while k > 0:
            if (k & 1) != 0:
                z = _GF2.mul(z,x)
            x = _GF2.mul(x,x)
            k >>= 1
        return z
    @staticmethod
    def divmodvect(xvec,dvec):
        nx = _GF2.bitlength(xvec[0])
        nd = _GF2.bitlength(dvec[0])
        i = nx - nd
        q = 0
        test = 1 << (nx-1)
        while i >= 0:
            if (xvec[0] & test) != 0:
                xvec = [x ^ (d << i) for (x,d) in zip(xvec,dvec)]
                q |= (1 << i)
            i -= 1
            test >>= 1
        return (q,xvec)
    @staticmethod
    def exteuc(a,b):
        '''
        based on Blankenship's algorithm
        return (g,x,y) such that g = gcd(a,b) and g = ax+by in GF(2)
        '''
        arow = [a,1,0]
        brow = [b,0,1]
        while True:
            (_,rrow) = _GF2.divmodvect(arow, brow)
            if rrow[0] == 0:
                break
            arow = brow
            brow = rrow
        return tuple(brow)

class GF2Poly(object):
    """A class representing the field of polynomials with coefficients in GF(2)"""
    def __init__(self, poly):
        self.poly = poly
    def __eq__(self, other):
        return self.poly == other.poly
    def __ne__(self, other):
        return self.poly != other.poly
    def __repr__(self):
        return 'GF2Poly(0b{0:b})'.format(self.poly)
    @staticmethod
    def getpoly(p):
        if isinstance(p, GF2Poly):
            return p.poly
        else:
            return p
    def __rmod__(self, x):
        x = GF2Poly.getpoly(x)
        return GF2Poly(_GF2.mod(x, self.poly))
    def __mod__(self, m):
        m = GF2Poly.getpoly(m)
        return GF2Poly(_GF2.mod(self.poly, m))
    @property
    def bitlength(self):
        return _GF2.bitlength(self.poly)
    @property
    def degree(self):
        return self.bitlength - 1

class GF2QuotientRing(GF2Poly):
    """A class representing the ring of polynomials with coefficients in GF(2),
    mod p(x), where p(x) is the characteristic polynomial having coefficients in GF(2).
    """
    @staticmethod
    def cast(p):
        try:
            return GF2QuotientRing(p.poly)
        except:
            return GF2QuotientRing(p)
    def wrap(self,x):
        return GF2Element(x, self)
    def mul(self,x,y):
        x = GF2Poly.getpoly(x)
        y = GF2Poly.getpoly(y)
        return GF2Poly(self.mulraw(x,y))
    def mulraw(self,x,y):
        return _GF2.mulmod(x,y,self.poly)
    def power(self, x, k):
        x = GF2Poly.getpoly(x)
        return GF2Poly(self.powraw(x,k))
    def powraw(self, x, k):
        z = 1
        m = self.poly
        while k > 0:
            if (k & 1) != 0:
                z = _GF2.mulmod(z,x,m)
            x = _GF2.mulmod(x,x,m)
            k >>= 1
        return z 
    def powvectraw(self, x,kvect):
        nk = len(kvect)
        zvect = [1]*nk
        alldone = False
        kvecttmp = [k for k in kvect]
        while not alldone:
            alldone = True
            for i in xrange(nk):
                k = kvecttmp[i]
                if (k & 1) != 0:
                    zvect[i] = _GF2.mulmod(zvect[i],x,self.poly)
                k >>= 1
                if k > 0:
                    alldone = False
                kvecttmp[i] = k
            x = _GF2.mulmod(x,x,self.poly)
        return zvect
    def powvect(self, x, kvect):
        return [self.wrap(z) for z in self.powvectraw(x, kvect)]
    def lshiftraw(self, x, k):
        '''
        return x << k mod m in GF(2)
        '''
        return self.mulraw(x,self.powraw(2,k))
    def lshift(self, x, k):
        return self.wrap(self.lshiftraw(x, k))    
    def rshiftraw(self, x, k):
        '''
        return x >> mod m in GF(2)
        TODO: works for primitive polynomials only
        '''
        r = self.poly >> 1
        return self.mulraw(x,self.powraw(r,k))
    def rshift(self, x, k):
        return self.wrap(self.rshiftraw(x, k))    

    def inv(self, x):
        '''
        return multiplicative inverse mod m, such that xy = 1 mod m in GF(2)
        '''
        x = GF2Poly.getpoly(x)
        (r,y,_) = _GF2.exteuc(x,self.poly)
        if r != 1:
            raise ValueError('%x and %x are not relatively prime but have a common factor of %x' % (x,self.poly,r))
        return y

def _bitlenlt(x,y):
    xy = x | y
    return (x << 1) < xy

def _gf2GaussJordan(A,b):
    '''
    solves for x, where Ax = b, in GF2.
    A is an n x n matrix; b is an n-element vector or an n x m matrix.
    '''
    n = A.shape[0]
    assert n == A.shape[1]
    try:
        s = b.shape
    except:
        s = (len(b),1)
    assert n == s[0]
    m = s[1]
    C = np.zeros((n,n+m),int)
    C[:,:n] = A
    if m == 1:
        C[:,n] = b
    else:
        C[:,n:] = b

    fails = []
    for j in xrange(n):
        # Find pivot
        p = None
        for i in xrange(j,n):
            if C[i,j] == 1:
                p = i
                break
        if p is None:
            fails.append(j)
            continue
        if p != j:
            C[p,:],C[j,:] = np.copy(C[j,:]),np.copy(C[p,:])
        for i in range(n):
            if i == j:
                continue
            if C[i,j] != 0:
                C[i,:] ^= C[j,:]
    if len(fails) > 0:
        raise ValueError('singular matrix: missing indices = %s' % fails)
    x = C[:,n:]
    return x

def _bitlencmp(x,y):
    xy = x | y
    if (x << 1) < xy:
        return -1
    elif (y << 1) <= x:
        return 1
    else:
        return 0    

def _bitsOf(p,n=None):
    def helper(p,n):
        if n is None:
            while p != 0:
                yield p & 1
                p >>= 1
        else:
            for _ in range(n):
                yield p & 1
                p >>= 1
    return tuple(helper(p,n))

def _gatherBits(v):
    x = 0
    p2 = 1
    for b in v:
        if b == 1:
            x |= p2
        p2 <<= 1
    return x
        




if fastgf2mulmod is not None:
    _gf2mulmod = fastgf2mulmod._gf2mulmod

def _gf2divmod(x,d):
    nx = _GF2.bitlength(x)
    nd = _GF2.bitlength(d)
    i = nx - nd
    q = 0
    while i >= 0:
        xnew = x ^ (d << i)
        if xnew < x:
            q |= (1 << i)
            x = xnew
        i -= 1
    return (q,x)

        


            

def _exteuc(a,b):
    '''
    based on Blankenship's algorithm
    return (g,x,y) such that g = gcd(a,b) and g = ax+by
    '''
    arow = [a,1,0]
    brow = [b,0,1]
    while True:
        (q,r) = divmod(arow[0],brow[0])
        if r == 0:
            break
        rrow = [r,arow[1]-q*brow[1],arow[2]-q*brow[2]]
        arow = brow
        brow = rrow
    return tuple(brow)


def _modinv(x,m):
    (r,y,_) = _exteuc(x,m)
    if r != 1:
        raise ValueError('%d and %d are not relatively prime but have a common factor of %d' % (x,m,r))
    return y

def _calculateCofactors(factors):
    '''
    given a vector of factors V,
    calculate vector V' such that the element-by-element
    product of V and V' is a vector of equal elements P 
    where P is the product of all factors V.
    
    In other words, each element of V' is the product of
    all elements of V except for the element in the corresponding position.
    
    for example: if V = [2,3,5,7] then V' = [105, 70, 42, 30] and P = 210. 
    '''
    n = len(factors)
    cofactors = [1]*n
    for (i,factor) in enumerate(factors):
        cofactors = [x * factor if i != j else x for (j,x) in enumerate(cofactors)]    
    return tuple(cofactors) 

def _pullFactor(x, testFactor):
    '''return tuple (n,r,y) such that
    x = r * y
    where y = testFactor^n
    and n is the largest possible integer
    such that r is not divisible by testFactor
    ''' 
    n = 0
    fpower = 1
    while True:
        (q,r) = divmod(x, testFactor)
        if r != 0:
            break
        n += 1
        fpower *= testFactor
        x = q
    return (n,x,fpower)

def _calculatePolynomialFactors(poly):
    '''
    calculate the factors of (2^n) - 1 where n is the degree of
    the desired polynomial
    '''
    n = _GF2.bitlength(poly)-1
    period = (1 << n) - 1
    m = period
    e1 = GF2Element(1, poly)
    if (e1 << m).value != 1:
        raise ValueError('%s not in primitive polynomial' % e1)
    return factorize(m)

def factorize(m):
    factors = []
    def testFactors():
        for f in [3,5,7,11,13,17,19,23,29,31,37,41,43,47]:
            yield f
        f = 53
        while True:
            yield f
            f += 2
    for f in testFactors():
        if m < f*f:
            break
        (_,m,fpower) = _pullFactor(m,f)
        if fpower > 1:
            factors.append(fpower)
    if m > 1:
        factors.append(m)
    return factors

def product(v):
    return reduce(lambda x,y: x*y, v, 1)
                
class GF2DiscreteLog(object):
    '''
    Facilitates computation of discrete logarithms, 
    after Clark and Weng (1994)
    "Maximal and Near-Maximal Shift Register Sequences:
    Efficient Event Counters and Easy Discrete Logarithms" 
    http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.44.118
    '''
    def __init__(self, poly, factors=None, maxtablesize=65536):
        self.qr = GF2QuotientRing.cast(poly)
        if factors is None:
            factors = _calculatePolynomialFactors(poly)
        factors = tuple(factors)
        self.factors = factors
        cofactors = _calculateCofactors(factors)
        self.cofactors = cofactors
        # verify factors
        e2 = self.qr.wrap(1)
        period = cofactors[0]*factors[0]
        self.period = period
        assert (e2 << period).value == 1
        lookup = []
        for (factor, cofactor) in zip(factors,cofactors):
            if factor > maxtablesize:
                raise ValueError('Factor %d exceeds maximum table size %d' % (factor, maxtablesize))
            g = e2 << cofactor
            gx = g
            glog = {1: 0}
            assert g.value != 1
            for i in xrange(1,factor):
                glog[gx.value] = i
                gx = gx * g
            assert gx.value == 1
            v = _modinv(cofactor, factor)
            lookup.append({'factor':factor, 'cofactor':cofactor, 'g':g, 'logtable':glog, 'v':v})
        self.lookup = lookup
    @staticmethod
    def _rem(y,item,period=None):
        cofactor = item['cofactor']
        r = item['logtable'][y]
        if period is None:
            return r*cofactor*item['v']
        else:
            return (((r*cofactor)%period)*item['v'])%period
    def log(self, x):
        if isinstance(x, GF2Element):
            if x.qr != self.qr:
                raise ValueError('Element %s has a different polynomial than %x' % (x, self.qr))
        else:
            x = self.qr.wrap(x)
        #yvect = [(x ** k).value for k in self.cofactors]
        yvect = self.qr.powvectraw(x.value, self.cofactors)     # this is faster, reuses powers of x
        r = [GF2DiscreteLog._rem(y,item,self.period) for (y,item) in zip(yvect,self.lookup)]
        return sum(r)%self.period 

def checkPeriod(poly, n):
    e2 = GF2Element(1,poly)
    if (e2 << n).value != 1:
        return 0
    factors = factorize(n)
    cofactors = _calculateCofactors(factors)
    for cf in cofactors:
        if (e2 << cf).value == 1:
            return cf
    return n
            
class GF2(object):
    '''
    classdocs
    '''
    def __init__(self, x):
        self.value = x
    def __repr__(self):
        return 'GF2(0b{0:b})'.format(self.value)
    def __add__(self, other):
        return GF2(self.value ^ other.value)
    def __sub__(self, other):
        return GF2(self.value ^ other.value)
    def __mul__(self, other):
        return GF2(_GF2.mul(self.value, other.value))
    def __mod__(self, other):
        return GF2(_GF2.mod(self.value, other.value))
    def __eq__(self, other):
        return self.value == other.value
    def __ne__(self, other):
        return self.value != other.value
    def expmod(self, k, m):
        return GF2(_GF2.powmod(self.value, k, m.value))
    
class GF2Element(object):
    def __init__(self, x, p):
        self.value = x
        if isinstance(p, GF2QuotientRing):
            self.qr = p
        else:
            self.qr = GF2QuotientRing(p)
        n = self.qr.degree
        self.fmt = 'GF2Element(0b{0:0%db},0x{1:x})' % n
    @property
    def poly(self):
        return self.qr.poly
    def _wrapraw(self, x):
        return GF2Element(x, self.qr)
    def __add__(self, other):
        return self._wrapraw(self.value ^ other.value)
    def __iadd__(self, other):
        self.value ^= other.value
        return self
    def __sub__(self, other):
        return self._wrapraw(self.value ^ other.value)
    def __isub(self, other):
        self.value ^= other.value
        return self
    def __mul__(self, other):
        return self.qr.wrap(self.qr.mulraw(self.value,other.value))
    def __imul__(self, other):
        self.value = self.qr.mulraw(self.value,other.value)
        return self
    def __pow__(self, k):
        return self.qr.wrap(self.qr.powraw(self.value, k))
    def __ipow__(self, k):
        self.value = self.qr.powraw(self.value, k)
        return self
    def __lshift__(self, k):
        return self.qr.lshift(self.value, k)
    def __ilshift__(self, k):
        self.value = self.qr.lshiftraw(self.value, k)
        return self
    def __rshift__(self, k):
        return self.qr.rshift(self.value, k)
    def __irshift__(self, k):
        self.value = self.qr.rshiftraw(self.value, k)
        return self
    def __repr__(self):
        return self.fmt.format(self.value, self.qr.poly)
    def __eq__(self, other):
        return self.value == other.value and self.qr == other.qr 
    def __ne__(self, other):
        return self.value != other.value or self.qr != other.qr
   
if __name__ == '__main__':
    x1 = GF2(0b101101)
    x2 = GF2(0b110110)
    x3 = GF2(0b101)
    print x1+x2
    print x2
    print x3
    print x2*x3
    e1 = GF2Element(0b110, 137)
    print e1*e1
    print e1 << 6
    print e1 << 127
    print e1 >> 125
    e2 = GF2Element(0b100, 137)
    print e2 ** 2
    print e2 ** 3
    print e2 ** 127
    
    b = 0b11010011
    a = 0b101101
    (g,x,y) = _GF2.exteuc(a,b)
    print g
    print x,y
    print _GF2.mul(a,x) ^ _GF2.mul(b,y)
    
    dlog5a = GF2DiscreteLog(0x23, [3,7])
    dlog5 = GF2DiscreteLog(0x25, [31])
    dlog8 = GF2DiscreteLog(0x11d, [3,5,17])
    dlog14 = GF2DiscreteLog(0x402b, [3,43,127])
    dlog16 = GF2DiscreteLog(0x1002d, [3,5,17,257])
    for dlog in [dlog5,dlog5a,dlog14]:
        e1 = dlog.qr.wrap(1)
        for i in xrange(25):
            x = e1 << i
            logx = dlog.log(x)
            print 'log %s = %d' % (x, logx) 
    for dlog in [dlog14,dlog16]:
        e1 = dlog.qr.wrap(1)
        for i in xrange(0,3000,33):
            x = e1 << i
            logx = dlog.log(x)
            print '%d: log %s = %d' % (i, x, logx)
            assert i == logx 
    e1 = dlog8.qr.wrap(1)
    print e1
    for i in xrange(9):
        print '1 << %d == %s' % (i,e1<<i)
    for i in xrange(100,109):
        print '1 << %d == %s' % (i,e1<<i)
    for i in xrange(200,209):
        print '1 << %d == %s' % (i,e1<<i)
    
    for p in [0x402b, 0x100000000065, 0x10000000000b7]:
        n = _GF2.bitlength(p)-1
        dlog = GF2DiscreteLog(p)
        f = dlog.factors
        print 'degree %d polynomial %x has factors %s' % (n,p,f)
    
    A = np.matrix([[1,0,1,0,1],[0,1,1,0,1],[0,1,0,0,1],[1,0,0,0,0],[0,0,1,1,0]], int)
    b = _bitsOf(0b10010)
    x = _gf2GaussJordan(A, b)
    print b
    print x
    print ((A*np.matrix(x))&1).transpose()

    b2 = np.matrix([[1,0,0,1,0],[0,1,1,1,1],[1,1,0,0,0]]).transpose()
    x2 = _gf2GaussJordan(A, b2)
    print x2
    print ((A*np.matrix(x2))&1)

    I5 = np.matrix(np.eye(5,dtype=int))
    Ainv = _gf2GaussJordan(A,I5)
    print Ainv
    print A*Ainv&1
    
    np.random.seed(123)
    I7 = np.matrix(np.eye(7,dtype=int))
    A7 = np.random.permutation(I7)
    Ainv = np.matrix(_gf2GaussJordan(A7,I7))
    print Ainv
    print A7
    print Ainv*A7&1
    
    poly = 0x20007
    qr = GF2QuotientRing(poly)
    import random
    r = random.Random(22)
    for i in range(1000):
        a = r.getrandbits(16)
        b = r.getrandbits(16)
        y1 = qr.mul(a,b)
        y2 = _GF2.mul(a,b) % qr
        print y1
        print y2
        assert y1 == y2
    print "Yay!"