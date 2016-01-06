
# Python2 compatibility
from __future__ import print_function, unicode_literals

# ensure we have Python 3 semantics from input, even in Python 2
try:
    input = raw_input
except NameError:
    pass

import sys

arr = [0] * 21000000

for name in sys.stdin:
    
    f=open(name.strip())


    for line in f:


        parts = line.strip().split()
            
        idx = int(parts[0])
            
        arr[idx] += int(parts[1])

    f.close()


for i in range(len(arr)):

    if arr[i] > 0:
        print(str(i) + " " + str(arr[i]))





