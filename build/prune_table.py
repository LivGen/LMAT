import sys

if (len(sys.argv) < 2):

    print "    Usage:  python " + sys.argv[0] + " <input-table-fn> <list-of-ids-table> [1] > <output-table-fn>"
    print " [1] indicates that the file is a taxonomy_rank.txt file"
    exit(0)

dict = {}

for line in open(sys.argv[2]):
    parts = line.split()
    
    dict[parts[0]] = "1"

parse_fmt = False

if (len (sys.argv) > 3 and sys.argv[3] == "1"):
    parse_fmt = True
    print >> sys.stderr, "taxonomy rank fmat"

for line in open(sys.argv[1]):
    
    tid = 0

    if (parse_fmt):
        parts = line.split(',')
        pp = parts[1].split('=')
        tid = pp[1]
        
    else:
        
        parts = line.split()
        tid =parts[0]
        
    if tid in dict.keys():
        print line.strip()
            




