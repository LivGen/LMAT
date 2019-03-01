
#! bin/perl 

#*This script creates the catalog.*; it takes the all_* file and parses the data giving name, sequence and taxonomic id
#of sequence. Arguments are : cat.pl all_*        

my $file = $ARGV[0];
my $filename="catalog.$file";
open(my $fh, '<',$file) or die "Could not open file 'filename' $!";
my $count=0;
open(my $sf,'>>',$filename)|| "file not created\n" ;
while( my $row = <$fh>){
	chomp $row;
	$count++;
	if ($row =~ m/\|gi\|(\d+)\|(.*?)\[.*sequence_id (\d+).*tax_node_id (\d+)/){
		my($name,$seq_id, $tax)=($2,$3,$4);
        	print {$sf} "$name\t$seq_id\t$tax\n";
	}elsif($row =~/Glued fragments of sequence (\d+) \((.*?),.*\[sequence_id (\d+).*tax_node_id (\d+)/){
        	my($name,$seq_id, $tax)=($2,$3,$4); 
        	print {$sf} "$name\t$seq_id\t$tax\n";
	}elsif($row=~/>gi\|(\d+)\|(.*?)\|(.*?),.*\[sequence_id (\d+).*tax_node_id (\d+)/){ 
	        my($name,$seq_id, $tax)=($3,$4,$5);
        	print {$sf} "$name\t$seq_id\t$tax\n";
	}elsif($row=~/>gi\|(\d+)\|(.*?)\|(.*?)\| (.*?),(.*?)/){
		my($seqid,$tax,$name)=($1,$3,$4);
		#print "$name\n";
		print {$sf} "$name\t$seqid\t$tax\n";
	}elsif($row=~/>(.*?)\s(.*?),(.*?)/){
		my($tax,$name,$seqid)=($1,$2,$count);
		print {$sf} "$name\t$\$seqid\t$tax\n";
	}elsif($row=~/>(.*?)\s(.*?)$/){
		my($tax,$name,$seqid)=($1,$2,$count);
		print {$sf} "$name\t$\$seqid\t$tax\n";
	}else{ 
		if($line=~ m/>/){  print "$count\n";}
 
	}
}
