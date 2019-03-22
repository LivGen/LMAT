import sys

tax_dict = {}

tid_dict = {}

# TODO need usage message


f = open ( sys.argv[1])

line = f.readline()
line = f.readline()
line = f.readline()

line = f.readline()

while (len(line) > 0):

    
    parts = line.strip().split()

    tax_dict[parts[0]] = parts[-1]

    line = f.readline()
    line = f.readline()

f.close()



for line in open (sys.argv[2]):

    parts = line.split()

    val = parts[0]

    tid_dict[val] = 1

    while (val != 0):

        parent = tax_dict[val]
        
        if parent in tid_dict.keys(): 
            break
        
        tid_dict[parent] = 1
        
        if parent == 1:
            break

        val = parent


f =  open ( sys.argv[1])

header1 = f.readline()
header2 = f.readline()
header3 = f.readline()

print header1.strip()
print header2.strip()
print header3.strip()

line = f.readline()

while (len(line) > 0):

    
    parts = line.split()

    if parts[0] in tid_dict.keys():
        print line.strip()

        line = f.readline()
        print line.strip()

    else:
        line = f.readline()

f.close()
