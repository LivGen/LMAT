import sys

parent_map = {}

species_level_lst = []

# load the map to parents

inf = open(sys.argv[1])

out_list = []

inf.readline()
inf.readline()

count = int(inf.readline().strip())

outf = open(sys.argv[3], "w")

#print str(count) +  " entries to load." 

while (inf):

	ids = inf.readline().strip()

	parts = ids.split()

	if (len(parts) < 2):
		break


	parent_map[parts[0]] = parts[-1]
	inf.readline()


print "Loaded Taxonomy - entries: " + str( len(parent_map) )

inf.close()

# now load the species and strains

inf = open(sys.argv[2])

for line in inf:

	parts = line.strip().split('\t')

	if (len(parts) > 3):
	
		
		zz = parts[0].split(',')
		
		zzz = zz[2].split('=')
		
		tid = zzz[1]
		

		pp = parts[-1].split(',')				
		
					
		if pp[0] == 'species':
			
			species_level_lst.append(tid)
			
			outf.write(tid + " " + tid + "\n")
		else:

			for x in reversed(parts[1:-1]):

				pp = x.split(',')				
			
				if pp[0] == 'species':
					out_list.append(tid)
					break


inf.close()

print "loaded lists:"
print str(len(species_level_lst)) + " species and " + str(len(out_list)) + " strains."


species_set = frozenset(species_level_lst)

errout = open("errors", "w")

for xx in out_list:

	running = True

	tid = xx

	while (running):

				
		if parent_map.has_key(tid):
			parent = parent_map[tid]

			if parent in species_set:
				running = False
				outf.write(xx + " " + parent + "\n")
			else:
				if tid == parent:
					running = False
					errout.write("same as parent: " + xx + " at " + tid + "\n")	
				else:
					tid = parent
				
		else:
			running = False
			errout.write("not in map: " + xx + " at "  + tid  + "\n")

outf.close()
