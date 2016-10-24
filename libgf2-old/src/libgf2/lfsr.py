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

from gf2 import GF2Element, GF2DiscreteLog, GF2Poly, GF2QuotientRing
        
def parity(k):
    i = 1
    while True:
        kshift = k >> i
        if kshift == 0:
            break
        k ^= kshift
        i <<= 1
    return k & 1

class LFSRBase(object):
    def __init__(self, poly):
        self.poly = poly;
        self.n = GF2Poly(poly).degree
    @property
    def degree(self):
        return self.n
    def _stateFromReverseIterator(self, it, atEnd=True):
        '''returns the state of an LFSR with default mask that produces the given output.
        if atEnd is True, the state in question is the state at the end of the output.
        Otherwise it is the state at the beginning of the output.
        '''
        n = self.degree
        s = 0
        # run backwards until the state has no unknown bits
        for bit in it:
            if bit != 0:
                s ^= self.poly
            s >>= 1
        # then run forwards again, if desired
        if atEnd:
            for i in xrange(n-1):
                s <<= 1
                if (s >> n) & 1 == 1:
                    s ^= self.poly
        return s
    
    def stateFromOutput(self, output, atEnd=True):
        '''returns the state of an LFSR with default mask that produces the given output.
        if atEnd is True, the state in question is the state at the end of the output.
        Otherwise it is the state at the beginning of the output.
        '''
        n = self.degree
        return self._stateFromReverseIterator(((output >> i & 1) for i in xrange(n)), atEnd)
    def stateFromOutputBits(self, output, atEnd=True):
        '''returns the state of an LFSR with default mask that produces the given output.
        if atEnd is True, the state in question is the state at the end of the output.
        Otherwise it is the state at the beginning of the output.
        '''
        n = self.degree
        return self._stateFromReverseIterator(output[n-1::-1], atEnd)

class LFSR(LFSRBase):
    '''
    Linear feedback shift register, with an output mask:
    at any given instant, the output is the parity of the state ANDed with the output mask.
    If no mask, then we just take the previous high bit.
    LFSR state s[k] = initstate * x^k 
    '''
    def __init__(self, poly, initstate=1, mask=None):
        LFSRBase.__init__(self, poly)
        self.testmask = 1 << self.n
        self.coefficientsMask = mask
        self.state = initstate
    def __iter__(self):
        return self
    def next(self):
        s = self.state << 1
        (b,nextState) = (0,s) if ((s & self.testmask) == 0) else (1,s^self.poly)
        if self.coefficientsMask is not None:
            b = parity(self.state & self.coefficientsMask)
        self.state = nextState            
        return b
    def __lshift__(self, k):
        newstate = GF2QuotientRing(self.poly).lshift(self.state, k).poly
        return LFSR(self.poly, initstate=newstate, mask=self.coefficientsMask)
    def __ilshift__(self, k):
        self.state = GF2QuotientRing(self.poly).lshift(self.state, k).poly
        return self
    def __rshift__(self, k):
        newstate = GF2QuotientRing(self.poly).rshift(self.state, k).poly
        return LFSR(self.poly, initstate=newstate, mask=self.coefficientsMask)
    def __irshift__(self, k):
        self.state = GF2QuotientRing(self.poly).rshift(self.state, k).poly
        return self
    def __repr__(self):
        return 'LFSR({0:x},{1:x})'.format(self.poly,self.state)

class LFSRAnalyzer(LFSRBase):
    def __init__(self, poly, factors=None):
        LFSRBase.__init__(self, poly)
        self.dlog = GF2DiscreteLog(poly, factors)
    def lookaheadCoefficients(self, k):
        '''returns coefficients a_j such that the sum of a_j*state[j] = output[n+k],
        for an LFSR with default mask'''
        e1 = GF2Element(1,self.poly)
        s = e1 << k
        return self.coefficientsFromState(s.value)
    def stateFromCoefficients(self, c):
        '''returns initial state of an LFSR with default mask, such that the output matches 
        the output of a masked LFSR with given coefficients and initial state of 1'''
         
        n = self.degree
        e1 = GF2Element(1,self.poly)
        out = 0
        for i in xrange(n):
            out <<= 1
            out |= parity((e1<<i).value & c)
        return self.stateFromOutput(out, atEnd=False)
    def timeshiftFromCoefficients(self, c):
        '''inverse of lookaheadCoefficients'''
        s = self.stateFromCoefficients(c)
        return self.dlog.log(s)
    def timeshiftFromState(self, s):
        return self.dlog.log(s)
    def coefficientsFromState(self, s):
        '''Inverse of stateFromCoefficients'''
        n = self.degree
        c = 0
        for i in xrange(n):
            b = (s >> (n-1)) << i
            c |= b
            s <<= 1
            if (s >> n) & 1 == 1:
                s ^= self.poly
        return c        
        

if __name__ == '__main__':
    print parity(0b1101)
    print parity(0b11011000010101101)
    print parity(0b11010000010101101)
    sr = LFSR(0b1101)
    for i in xrange(17):
        b = sr.next()
        print (i,b,sr)
    print 'sr >> 17 = %s' % (sr >> 17)
    sr = LFSR(0x10000000000b7)
    print sr << 17
    print sr << 19
    print sr << 123
    sr <<= 1
    sr <<= 123
    print sr
    print sr >> 124
    print sr << 12
    print sr << 13
    poly = 0x8003
    poly = 0x9091
    lfsrAnalyzer = LFSRAnalyzer(poly)
    n = GF2Poly(poly).degree
    e1 = GF2Element(1,poly)
    for k in [100,110,200,556,9171]:
        o = 0
        for i in xrange(n):
            y = e1 << (i+k)
            o = (o << 1) | (y.value >> n-1)
            print '%d: %s' % (i+k, y)
        print 'output = %s' % bin(o)
        print 'state  = 0b{0:b} (expected 0b{1:b})'.format(lfsrAnalyzer.stateFromOutput(o,True), (e1 << (k+n-1)).value)
        c = lfsrAnalyzer.lookaheadCoefficients(k)
        print bin(c)
        print '{0:b} --> timeshift={1:d} (expected {2:d})'.format(c,lfsrAnalyzer.timeshiftFromCoefficients(c),k)
        S = [e1 << j for j in xrange(30,50)]
        xexact = [(e1 << j).value >> (n-1) for j in xrange(30+k,50+k)]
        xpredict = [parity(s.value & c) for s in S]
        print 'k=%d\n    xexact=%s\n  xpredict=%s' % (k,xexact,xpredict)

    ecomp = e1 << ((1 << n)-n)
    for c in [123, 942, 1000, 2107, 12280, 15092, 21038, 16384, 32767]:
        s = lfsrAnalyzer.stateFromCoefficients(c)
        c2 = lfsrAnalyzer.coefficientsFromState(s)
        print 'c=%d, s=%s, c2=%d' % (c,s,c2)
        
   