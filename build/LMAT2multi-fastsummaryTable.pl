#!/usr/bin/perl

# LMAT2multi-fastsummaryTable.pl  -in flst  -out outfile  -min_reads 1 -min_frac 0.00 -min_score 0.7 

use Getopt::Long; 

my $min_reads=1;  # minimum number of raw reads to count in output
my $min_frac=0.000; # minimum fraction of total reads to count in output
my $file_list;  # file with two columns, listing the full path to each LMAT output file( .fastsummary or .species or .genus or .plasmid files), one path per line.
my $out;  # table summarizing output of files in -in file. taxID's are rows, 
my $avg_score_min=0.5;

GetOptions ( "min_reads:i" =>\$min_reads,
             "min_frac:f" =>\$min_frac,
             "min_score:f" =>\$avg_score_min,
             "in:s" =>\$file_list,
             "out:s" => \$out
    );
        
if($file_list !~ /\w/ || $out  !~ /\w/) {
    die "\nUsage: $0 -in fileList -out outfile \n\t-min_score <minimum average LMAT score per read for the tax id to be tabulated, default 0.5>\n\t-min_frac <minimum fraction of total reads in sample for a tax id to be tabulated, default 0>\n\t-min_reads <minimum number reads in sample for a tax id to be tabulated, default 1>\n\nExample: LMAT2multi-fastsummaryTable.pl  -in fileList  -out TaxIDBySample.table.txt -min_frac 0.0001 -min_score 0.5 -min_reads 5\n\nfileList is a file listing the files to include in the output table, where each line of fileList has two tab-delimited columns: 1) sample name and 2) full path name to the LMAT .fastsummary, .species, .genus or .plasmid file for that sample\n
If outfile has 'megan' in the file name, then it creates a .megan summary file, with the MEGAN \@headers at top and raw read counts as table entries. Otherwise it omits MEGAN \@headers, prints organism name in the first column, taxid in second col, and that organism's fraction of total reads as table entries. Taxid's are in the rows and samples are in the columns. This program requires the LMAT *.log file which gives the total number of reads for that sample to be in the same directory as the LMAT .fastsummary|.species|.genus|.plasmid file.\n";
}


open IN,"$file_list" or die "Cannot open $file_list: $!\n";
my $count=0;
my @File=();
my @Name=();
while ( my $line=<IN>) {
    chomp $line;
    ($name,$file)=split/\t/,$line;
    $Name[$count]=$name;
    $File[$count]=$file;
    $count++;
    print "$name\t$file\n";
}

close IN or warn "Cannot close $file_list: $!\n";

my $total_reads=0;
my %data_by_org=();

my @num_reads_array=();

my %valid;
for (my $i=0 ; $i<@File ; $i++) {
    my $file=$File[$i];
    open IN,"$file" or warn "Cannot open $file: $!\n";
    my ($avg_score,$wrc,$rc,$taxid,$org);
    while (my $line=<IN>) {
	    chomp $line;
       if ($file =~ /fastsummary$/) {
          ($wrc,$rc,$taxid,$org) = split/\t/,$line;
          $avg_score=$wrc/$rc;
       } else {
           ($avg_score,$wrc,$rc,$taxid,$org,@junk) = split/\t/,$line;
       }
       if( $rc >= $min_reads && $avg_score >= $avg_score_min ) {
           $valid{$taxid}=1;
       }
    }
}

for (my $i=0 ; $i<@File ; $i++) {
    my $file=$File[$i];

    $file =~ s/\*//;
    print "file: $file\n";
    my $file_nopath=`basename $file`;
    chomp $file_nopath;
    print "file_nopath: $file_nopath\n";
    my $prefix;
    if ($file_nopath =~ /^(\S+)_output/){
	$prefix=$1; 
    }
    my $thisdir=`dirname $file`;
    chomp $thisdir;
    print "dir: $thisdir\n";
    my $cmd="grep \" reads in\" $thisdir/$prefix*_output.log "  ;
    print "$cmd\n";
    my $z=`$cmd`;
    my $num_reads=0;
    if ($z=~ /(\d+)/) {
	$num_reads=$1;
    }
    print "num reads: $num_reads\n";
    $num_reads_array[$i]=$num_reads;
    $total_reads +=$num_reads;
    open IN,"$file" or warn "Cannot open $file: $!\n";
    my ($avg_score,$wrc,$rc,$taxid,$org);
    while (my $line=<IN>) {
	chomp $line;

	if ($file =~ /fastsummary$/) {
	    ($wrc,$rc,$taxid,$org) = split/\t/,$line;
	    $avg_score=$wrc/$rc;
	} else {
	     ($avg_score,$wrc,$rc,$taxid,$org,@junk) = split/\t/,$line;
	}
	#if ($rc/$num_reads >= $min_frac && $rc >= $min_reads && $avg_score >= $avg_score_min  ) {
	if ($rc/$num_reads >= $min_frac && $avg_score >= $avg_score_min && $valid{$taxid}  ) {
	    # print "$org $rc ",$rc/$num_reads,"\n";
	    $org =~ s/^\s+//;
	    $org =~ s/\s+$//;

	    $observed_orgs{$taxid}=$org;
	    if ($out =~ /megan/i) {
		$data_by_org{$taxid}{$i}=$rc;
	    } else {
		$data_by_org{$taxid}{$i}=$rc/$num_reads;
	    }
	}
    }
    close IN;
} # for (my $i=0 ; $i<@File ; $i++) {
 
open OUT,">$out";
if ($out =~ /megan/i) {
    print OUT "\@Creator\tLMAT\n";
    print OUT "\@CreationDate\t",`date`;
    print OUT "\@ContentType\tSummary4\n";
    print OUT "\@Names\t",join("\t",@Name),"\n";
    #print OUT "\@Uids\t@sra_list\n";
    print OUT "\@Sizes\t",join("\t",@num_reads_array),"\n";
    print OUT "\@TotalReads\t$total_reads\n";
    #print OUT "\@Collapse\n";
    print OUT "\@Algorithm\tTaxonomy\tlmat\n";
    #print OUT "\@NodeStyle\t
} else {
     print OUT "Sequence\ttaxid\t",join("\t",@Name),"\n";
}
foreach my $taxid (sort {$a <=> $b} keys %observed_orgs) {
    if ($out =~ /megan/i) {
	print OUT "TAX\t$taxid";
    } else {
	print OUT "$observed_orgs{$taxid}\t$taxid";
    }
    for (my $i=0 ; $i<@File ; $i++) {
	if (!defined $data_by_org{$taxid}{$i} ) {
	    $data_by_org{$taxid}{$i}=0;
	}
	print OUT "\t$data_by_org{$taxid}{$i}";
    }
    print OUT "\n";
}
close OUT;


