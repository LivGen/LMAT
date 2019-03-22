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
        print str(i) + " " + str(arr[i])





