#!/usr/bin/env perl
#
# dustmasker_wrapper.pl 
#
# This script is modified version of "parse_ftp_tracedb_data.pl" to suit for processing PacBio sequencing reads.
#

use strict;
use warnings;
use FAlite;
use Getopt::Long;

###############################
#
#  Command-line options
#
################################

my $dir;              # directory to look for *.fasta files
my $min_bases;        # minimum number off bases that you need in a sequence after clipping to keep it
my $max_n;	          # what is the maximum percentage of N's allowed in the masked sequence
my $verbose;          # turn on extra output - e.g. errors for reads that were rejected
my $help;             # display help

GetOptions ("dir:s"             => \$dir,
			"min_bases:i"       => \$min_bases,
			"max_n:f"           => \$max_n,
			"verbose"           => \$verbose,
			"help"              => \$help);



###############################################################
#
# Set some default values and check commmand line options
#
###############################################################

$min_bases = 1000    if (!$min_bases);
$max_n = 15         if (!$max_n);

my $usage = "
usage: dustmakser_wrapper.pl -dir <directory> <options>
  -dir <directory> : specify a directory containing fasta & anc files for a single species
  -min_bases <int> : minimum number of bases required in a read *after* clipping. Default = 1000
  -max_n <percentage> : maximum allowable percentage of N characters in final sequence. Default = 15%
  -help : this help
\n";

die "$usage" if (!$dir or $help);

# add a trailing slash if none was specified to path (this is a bit of a kludge)
($dir .= "/") if ($dir && $dir !~ m/\/$/);



###############################################################
#
# Main variables for script
#
###############################################################

# need to keep track of:
# 1) total reads parsed
# 2) total bases in the reads
# 3) how many reads are rejected for being too short after clipping low quality bases, or have too few non-ATCG characters after dusting
# 4) how many reads are rejected for containing > $max_n % Ns (calculated both before and after dusting)
my $total_reads = 0;
my $total_bases = 0;
my $total_rejected_high_n = 0;
my $total_rejected_too_short = 0;



################################################################################
#
#                          M A I N   P R O G R A M
#
################################################################################

print STDERR "\n# $0 started at ", `date`, "\n";

# open a file which will contain details of any files that have been processed (so they can be ignored in future)
my $processed_file_name = "parameters.txt";
open(PROCESSED,">$processed_file_name") or die "Can't create $processed_file_name\n";


my @fasta_files = glob("${dir}*.fasta");
my $outdir = "preprocessed";
if(! -e $outdir) {
	system("mkdir -p $outdir");
}

foreach my $fasta_file (@fasta_files){
	# reset file-level counters
	my $file_total_reads = 0;
	my $file_total_bases = 0;
	my $file_rejected_high_n = 0;
	my $file_rejected_too_short = 0;


	# extract species name and file number
	my $fasta_file_name = $fasta_file;
	$fasta_file_name =~ s/^$dir//;
	$fasta_file_name =~ m/^(.+)\.([0-9]+).fasta/;
	my ($species,$file_number) = ($1,$2);

	# create output file name to be used for a couple of output files
	my $output_file = "$outdir/${species}_processed_traces.${file_number}.fa";
	my $dust_temp_file = "${output_file}.tmp";

	###############################################################
	#
	# Run DUST filter over file and create final output file
	#
	###############################################################

	# now need to run the dust program to filter low complexity regions from reads
	print STDERR "Filtering sequences with dustmasker\n";
	system("dustmasker -in $fasta_file -out $dust_temp_file -outfmt fasta") && die "Can't run dust program on $output_file.predust\n";
	


	###############################################################
	#
	# Process DUSTed file and apply final set of selection criteria
	#
	###############################################################


	# open dusted file and filter sequences
	open(DUST,"<$dust_temp_file") or die "Can't open $dust_temp_file\n";	
	my $file = new FAlite (\*DUST);

 	# open main output file to hold processed sequences
	open(OUT,">$output_file") or die "Can't create $output_file\n";	


	print STDERR "Applying sequence filters to dustmasked file\n";
	
    while (my $entry = $file->nextEntry) {
    	my $header = $entry->def;
		my ($ti) = $header =~ m/>(\d+)/;
		my $seq = $entry->seq;
		$seq =~ s/a|c|g|t/N/g;
		$seq = uc($seq);        
		my $length = length($seq);
		my $n = $seq =~ tr/N/N/;
		my $non_n = $length - $n;
		my $percent_n = ($n / $length * 100);
		
		$total_reads++;
		$file_total_reads++;
		$total_bases += $length;
		$file_total_bases += $length;

		# reject sequence if there is not enough space to have a tandem repeat
		if ($non_n < $min_bases){
			print STDERR "ERROR: $ti contains fewer than $min_bases bases that are not Ns in its sequence after dusting\n" if ($verbose);
			$file_rejected_too_short++;			
			$total_rejected_too_short++;			
			next;
		}
		# reject if there are too many N's now overall (some Ns may have been in sequence before dusting)
		elsif ($percent_n > $max_n){
			print STDERR "ERROR: $ti contains more than ${max_n}% Ns in its sequence after dusting\n" if ($verbose);
			$file_rejected_high_n++;			
			$total_rejected_high_n++;			
			next;
		}
		# if we get here, we have a sequence which is OK and we can print to the final output file
		else{
			print OUT "$header\n$seq\n";
		}
    }

	close(DUST);
	unlink("$dust_temp_file") or die "Can't remove $dust_temp_file\n";
	close(OUT);

	###############################################################
	#
	# Print summary stats for this file
	#
	###############################################################
	
	print STDERR "$file_rejected_high_n reads were rejected for containing more than $max_n% Ns after dusting\n" if ($file_rejected_high_n > 0);
	print STDERR "$file_rejected_too_short reads were rejected for being having less than $min_bases non-N characters after dusting\n" if ($file_rejected_too_short > 0);


	# update processed file information, include settings used to process files
	print PROCESSED "$fasta_file min bases=$min_bases max \%n=$max_n\n";
	
	print STDERR "\n";
	
}

close(PROCESSED);



###############################################################
#
# Print summary stats for all processed files
#
###############################################################

print STDERR "\n\n======================================================\n\n";
print STDERR "TOTAL: Processed $total_reads reads containing $total_bases nt\n";
print STDERR "TOTAL: $total_rejected_high_n reads were rejected for containing more than $max_n% Ns after dusting\n";
print STDERR "TOTAL: $total_rejected_too_short reads were rejected for having less than $min_bases non-N characters after dusting\n";
print STDERR "\n======================================================\n\n";

print STDERR "# $0 started at ", `date`, "\n";
exit(0);
	
