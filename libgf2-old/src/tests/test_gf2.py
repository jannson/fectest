'''

Copyright 2015 Jason M Sachs

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
import pytest

from libgf2.gf2 import GF2Element, GF2Poly, GF2QuotientRing

class Test_Poly8(object):
    def setup(self):
        self.qr = GF2QuotientRing(0x187)
    def test_value(self):
        for p in [0x3, 0x87, 0x0, 0x123, 0x12345]:
            assert GF2Poly(p).poly == p
    def test_mod1(self):
        assert 0x100 % self.qr == GF2Poly(0x87)
        assert 0x187 % self.qr == GF2Poly(0)
        assert 0x200 % self.qr == GF2Poly(0x89)
    def test_mod2(self):
        assert self.qr % 0x33 == GF2Poly(0x1f)
        assert self.qr % self.qr == GF2Poly(0)

class Test_Element(object):
    POLY = 0x187
    def setup(self):
        self.e = GF2Element(1, self.POLY)
    def test_lshift(self):
        e = self.e
        assert e.value == 1
        assert e.poly == self.POLY
        assert (e << 1).value == 2
        assert (e << 2).value == 4
        assert (e << 40).value == 0x62
        assert (e << 255).value == 1
        assert (e << 254).value == 0xc3
        assert (e << 253).value == 0xa2
    def test_rshift(self):
        e = self.e
        assert (e >> 1).value == 0xc3
        assert (e >> 2).value == 0xa2
        assert (e >> 215).value == 0x62
        assert (e >> 253).value == 4
        assert (e >> 254).value == 2
        assert (e >> 255).value == 1

