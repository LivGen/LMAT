import sys

if (len(sys.argv) < 2):
	print "Usage:  " + sys.argv[0] + " taxonomy_rank_fn \n"
	print "prints to stdout\n"
	exit(0)

inf = open(sys.argv[1])

items = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "superkingdom"]

rank = [15, 14, 12, 10, 8, 6, 4, 2]

MAX_RANK = 1

match = -1

inter_node = False

for line in inf:



	parts = line.split('\t')

	match = -1

	inter_node = False

	if (len(parts) > 1):
	
		xx = parts[1:]

		xx.reverse()

		for it in xx:
		
			pp = it.split(",")
			
			if (pp[0] in items):

				match = rank[items.index(pp[0])]
				if (inter_node == True):

					match = match + 1
				break
			
			else:
				inter_node = True
					
		if (inter_node and (match == -1)):
			match = MAX_RANK
	
	
	
		pp = parts[0].split(',')

		kt = pp[1].split("=")

		print kt[1] + " " + str(match)



	



