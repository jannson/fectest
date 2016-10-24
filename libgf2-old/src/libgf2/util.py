'''
Created on Feb 10, 2014

@author: jmsachs

Copyright 2013-2014 Jason M. Sachs

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

def reversed_bits_generator(bitsx, n=None):
    reversed_bits = 0
    n = float('inf') if n is None else n
    while n > 0:
        n -= 1
        bit = bitsx & 1
        bitsx >>= 1
        reversed_bits <<= 1
        reversed_bits |= bit
        yield reversed_bits

def parity(bitsx):
    '''parity of the bits contained in bitsx'''
    k = 1
    while True:
        y = bitsx >> k
        if y == 0:
            break
        bitsx = bitsx ^ (bitsx >> k)
        k <<= 1
    return bitsx & 1

def berlekamp_massey(bits, N, verbose=False):
    '''
    compute minimal LFSR that produces the bits in question
    bit ordering is little-endian
    '''
    
        
    b = 1
    c = 1
    L = 0
    rbitsgen = reversed_bits_generator(bits)
    m = -1
    
    for n in xrange(N):
        rbits = rbitsgen.next()
        # rbits = the bottom n bits, in reversed order
        mask = (1 << (L+1)) - 1
        d = parity(mask & rbits & c)
        if verbose:
            print 'n={0:d}, m={1:d}, rbits={2:b}, c={3:b}, b={4:b}, L={5:d}'.format(n, m, rbits, c, b, L)
        if d == 1:
            t = c
            nminusm = n - m
            c ^= b << nminusm
            if L <= n/2:
                L = n+1 - L
                m = n
                b = t
    return (c, L)



if __name__ == '__main__':
    pass