import numpy as np

w = 8
text = "hello world"
elen = len(text)
vlen = 6
emax = elen + vlen

d = np.array([ord(i) for i in text], dtype=np.int)
identify = np.identity(elen)
vander = np.vander(np.array([i+1 for i in range(elen)], dtype=np.int), vlen, increasing=True)
b = np.append(identify, vander.T, axis=0)

#print "b:"
#print b

print "origin d:"
print d

#encode
encode = np.dot(b, d.T)
print "encode:"
print encode.astype(int) % 255
