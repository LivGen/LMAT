#! /bin/sh

#Put seq.mapping.multi as input

if grep ^"32630" $1; then
	echo "32630 in mapping file"
else
	echo "Adding synthetic construct information"
	grep ">" /usr/mic/bio/blastdb/vector/vec.042014.fa| perl -ne 'print "32630\t-1\t-1$_";' >> $1
fi


