import pyscamp as mp
import sys

f1 = sys.argv[1] 
f2 = sys.argv[2]
sublen = sys.argv[3]
o1 = sys.argv[4] 
o2 = sys.argv[5]
with open(f1, 'r') as f1p:
    f1c = f1p.readlines()

for i in range(len(f1c)):
    f1c[i] = float(f1c[i])

with open(f2, 'r') as f2p:
    f2c = f2p.readlines()

for i in range(len(f2c)):
    f2c[i] = float(f2c[i])

profile, index = mp.abjoin(f1c, f2c, int(sublen))

with open(o1, 'w') as o1c:
    for i in profile:
        o1c.write(str(i)+'\n')

with open(o2, 'w') as o2c:
    for i in index:
        o2c.write(str(i)+'\n')