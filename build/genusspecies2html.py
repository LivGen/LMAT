import sys

if len(sys.argv) < 4:
    print "Usage:  python " + sys.argv[0] + " <species-filename> <genus-filename> <taxonomy_rank-txt-file>"
    print "outputs to stdout"

    exit(0)



speciesf=open(sys.argv[1])

genusf=open(sys.argv[2])

rank_f=open(sys.argv[3])

dict = {}

for line in rank_f:

    parts=line.split("\t")

    pp=parts[0].split(",")

    ppp=pp[1].split("=")

    tid = ppp[1]

    dict[tid] = parts[1:]

    
HOMO="#FFCCCC"
VIRSTR="#FF0000"
VIRSPE="#990033"
VIRGEN="#660000"
FUNSTR="#CCCCCC"
FUNSPE="#999966"
FUNGEN="#666633"
PLASMID="#00CC00"
BACSTR="#00000FF"
BACSPE="#0033CC"
BACGEN="#000066"

EUKSTR="#FF00FF"
EUKSPE="#CC00CC"
EUKGEN="#990099"



print "<html>"
print "<body>"

print "<table border=2>"


speciesarr = []
allarr = []

for line in speciesf:


    parts = line.split("\t")

    item = []

    if len(parts) > 5:

        item= parts[4:8]

        item.append("strain")
    else:
        item=parts[0:4]

        item.append("species")

    item.append(dict[item[2]])

    speciesarr.append(item)
  
i = 0 

for line in genusf:
    
    parts = line.split("\t")

    item = parts[0:4]

    item.append("genus")

    item.append(dict[parts[2]])

    while  i < len(speciesarr) and float(speciesarr[i][0]) >= float(parts[0]) :
        allarr.append(speciesarr[i])
        i = i + 1

    allarr.append(item)
                          


for data in allarr:

#    print data

    fgcolor="#FFFFFF"

    if len(data) < 6:
        color = "#FFFFFF"
        fgcolor="#000000"
    elif "plasmid" in data[5][-1]:
        color = PLASMID
        fgcolor="#000000"
    elif "Homo" in data[5][-1]:
        color = HOMO
        fgcolor="#000000"
    elif "Virus" in data[5][0]:
        if data[4] == "genus":
            color = VIRGEN
        elif data[4] == "species":
            color = VIRSPE
        else:
            color=VIRSTR
    elif len (data[5]) > 1 and "Bacteria" in data[5][1]:
        if data[4] == "genus":
            color = BACGEN
        elif data[4] == "species":
            color = BACSPE
        else:
            color=BACSTR
    elif len(data[5]) > 3 and "Fungi" in data[5][3]:

        if data[4] == "genus":
            color = FUNGEN
        elif data[4] == "species":
            color = FUNSPE
        else:
            color=FUNSTR
            fgcolor="#000000"
    else:
        if data[4] == "genus":
            color = EUKGEN
        elif data[4] == "species":
            color = EUKSPE
        else:
            color=EUKSTR
    
    print "  <tr>" 



        
    rank = data[4]

    for n in data[0:3]:

        print '    <td bgcolor="' + color + '">'
        
        print "<b><font color="+ fgcolor +">" + n + "</font></b></td>"



        ident = ""
        
    if ',' in data[3]:
        pair= data[3].split(',')


            
        ident=pair[1]

    else:

        
        ident = data[3]
        
    print '    <td bgcolor="' + color + '">'    
    print "<b><font color="+ fgcolor +">" + rank + "</font></b></td>"
    print '    <td bgcolor="' + color + '">'
    print "<b><font color="+ fgcolor +">" + ident + "</font></b></td>"

    print "  </tr>"

print "</table>"


print "</body>"
print "</html>"
