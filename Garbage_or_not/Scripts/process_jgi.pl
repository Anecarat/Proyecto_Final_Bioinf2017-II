#!/bin/env perl

#    Copyright (C) [2015] [Sur Herrera Paredes]
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

# Usage:
#	$ module load perl/5.8.9
#	$ module load r/3.0.1
#	$ process_jgi.pl -i <ncontam_nphix_phile> -t <trim_length>

# Load dependencies
use warnings;
use strict;
#use Getopt::Long qw(:config no_ignore_case bundling);
use Getopt::Long qw(:config no_ignore_case);
#use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;

# Declare global variables.
my ($bindir,$infile,$trim_length,$diagnosticsdir,$Fprimer,
	$Rprimer,$keep_N,$mapping_file,$help,$QC_only,
	$min_length,$max_length);

# Get options
set_defaults(\$bindir,\$trim_length,\$diagnosticsdir,
		\$keep_N,\$mapping_file,\$help,\$QC_only,
		\$infile,\$min_length,\$max_length);
my $opts = GetOptions("infile|i=s" => \$infile,
			"trim|t=i" => \$trim_length,
			"bindir|b=s" => \$bindir,
			"diagdir|d=s" => \$diagnosticsdir,
			"Fprimer|F=s" => \$Fprimer,
			"Rprimer|R=s" => \$Rprimer,
			"keep_N|k" => \$keep_N,
			"mapping_file|m=s" => \$mapping_file,
			"QC_only" => \$QC_only,
			"min_length|l=i" => \$min_length,
			"max_length|L=i" => \$max_length,
			"help|h" => \$help);
			
pod2usage(-exival => 1, -verbose => 2) if $help;

# If the option -QC_only was passed, then we skip the first two steps.
unless($QC_only){
	# IF 
	pod2usage(-exival => 1, -verbose => 2) unless $infile;
	
	# Preprocess and merege seuquences
	merge($infile,$trim_length,$bindir);
	
	# Plot diagnostics
	diagnose_run('out.extendedFrags.fastq',
			'ncontam_nphix.forward.fastq',
			'ncontam_nphix.reverse.fastq',
			'out.hist',
			$diagnosticsdir,
			$bindir);
}

# Perform quality control and print out high quality reads
QC_reads($Fprimer,$Rprimer,$keep_N,$mapping_file,
	$min_length,$max_length);


##############SUBROUTINES####################
sub set_defaults{
	# Function that sets default values for program options.
	# It takes a set of references and asigns the default values.
	my ($bindir_ref,$trim_length_ref,$diagnosticsdir_ref,
		$keep_N_ref,$mapping_file_ref,$help_ref,
		$QC_only_ref,$infile_ref,$min_length_ref,$max_length_ref) = @_;
	$$bindir_ref = '/proj/dangl_lab/sur/immune_code/';
	$$trim_length_ref = 165;
	$$diagnosticsdir_ref = './diagnostics/';
	$$keep_N_ref = '';
	$$mapping_file_ref = '';
	$$help_ref = '';
	$$QC_only_ref = '';
	$$infile_ref = '';
	$$min_length_ref = 0;
	$$max_length_ref = 0;
	return;
}

sub merge{
	# TO DO: THE SCRIPT THAT I CALL HERE MUST BE TURNED INTO A FUNCTION
	# FOR THIS SCRIPT TO BE COMPLETEL SELF CONTAINED. TRIMMING WOULD
	# REQUIRE EXTRA DEPENDENCIES.
	my ($ncontam_nphix_file,$trim_length,$bindir) = @_;
	my $bin = "$bindir/merge_miseq_jgi.pl";
	die "Can't find executable script for merging ($bin)" unless -f $bin;
	my $command = "$bin $ncontam_nphix_file $trim_length";
	execute($command);
	return;
}

sub diagnose_run{
	# Function that creates the parameters file, and sets up an output directory for
	# the R script that plots some features of the run.
	my ($merged_reads,$forward_reads,$reverse_reads,$merged_hist,$outdir,$bindir) = @_;
	my $merged_stats = "$outdir/merged_qscores.txt";
	my $forward_stats = "$outdir/forward_qscores.txt";
	my $reverse_stats = "$outdir/reverse_qscores.txt";

	# Create output directory and save summary statistics to it
	mkdir $outdir unless -d $outdir;
	get_stats($merged_reads,$merged_stats,$bindir);
	get_stats($forward_reads,$forward_stats,$bindir);
	get_stats($reverse_reads,$reverse_stats,$bindir);

	# Create parameters file for plotting.
	my $params_file = "$outdir/params.txt";
	open(OUT,'>',$params_file) or die "Can't create $params_file ($!)";
	print OUT "param\tvalue\n";
	print OUT "file.merged\t$merged_stats\n";
	print OUT "file.forward\t$forward_stats\n";
	print OUT "file.reverse\t$reverse_stats\n";
	print OUT "file.hist\t$merged_hist\n";
	print OUT "outdir\t$outdir\n";
	close OUT;

	plot_diagnostics($params_file,$bindir);
	return;
}

sub QC_reads{
	# NOTE: TO FULLY DECIDE ON THIS, I REQUIRE AN ANALYSIS OF DIFERENT PLATES, SHOWING THE DISTRIBUTIONS
	# OF Q-SCORES.
	# There are 6 filters that should be implemented:
	# Fprimer forward primer sequence to be removed. DONE
	# Rprimer reverse primer sequence to be removed. DONE
	# N remove sequences with ambiguous bases (anything but ACGT). DONE
	# mean_score remove sequences with a mean Q-score below a certain threshold
	# short remove sequences shorter than a threshold. DONE
	# long remove sequwnces longer thana a threshold. DONE
	# lowqual remove sequences with a single base below a Q-score threshold.
#	my ($Fprimer,$Rprimer,$keep_N,$remove_mean_score,$remove_short,$remove_long,$remove_lowqual) = @_;
	my ($Fprimer,$Rprimer,$keep_N,$mapping_file,$min_length,$max_length) = @_;
	my ($insert,%Map,%Counter); # this coold be done with OOP, with an object that has all the information as attributes

	print "Starting quality control step\n";	
	
	if($mapping_file){
		print "\tReading mapping file...\n";
		open(MAP,$mapping_file) or die "Can't open $mapping_file ($!)";
		while(<MAP>){
			chomp;
			my @line  = split(/\t/,$_);
			$Map{$line[1]} = $line[0];
			$Counter{$line[1]} = 0;
		}
		print "\tMapping file processed...\n";
	}else{
		print "\tNo mapping file provided. Read ID's will remain the same\n";
	}
	
	# Get regular expression for the primer sequence
	my ($regexp,$pattern) = make_regexp($Fprimer,$Rprimer);
	die "You must provide at least one primer!" unless $pattern;
	print "\tPrimer matching regular expression set:\n\t==$regexp==\n";
	
	# Extract sequences that match the primer sequence.
	my ($i,$match_regexp,$N_seqs,$short_seqs,$long_seqs,$written_seqs,$unknown_barcode);
	my $IN = Bio::SeqIO->new(-format => 'fastq', -file => 'out.extendedFrags.fastq');
	#my $OUT = Bio::SeqIO->new(-format => 'fastq',-file => ">merged.qc.fastq");
	open(OUT,">merged.qc.fasta") or die "Can't create file ($!)";
	print "\tReading input file out.extendedFrags.fastq...\n";
	$i = $match_regexp = $N_seqs = $short_seqs = $long_seqs = $written_seqs = $unknown_barcode = 0;
	
	while(my $Seq = $IN->next_seq){
		# Extract sequence information.
		my $seq = $Seq->seq;
		#my @quals = @{$Seq->qual};
		my $id = $Seq->id;
		print "\tProcessed $i sequences...\n" unless ++$i % 100000;
		
		if ($mapping_file){
			$id = process_read_id($id,\%Map,\%Counter);
			if($id eq ''){
				$unknown_barcode++;
				next;
			}
		}
		
		my $insert_length;
		#$seq = reverse_complement($seq) if $reverse_complement; # I don't think I need this
		if($seq =~ /$regexp/){
			$match_regexp++;
			if($pattern == 1){
				$insert = $3;
				$insert_length = length $3;
				#print "($1) ($2) ($insert_length)\n";
			}elsif($pattern == 2){
				$insert = $1;
				$insert_length = length $1;
				#print "($insert_length) ($2) ($3)\n";
			}elsif($pattern == 3){
				$insert = $3;
				$insert_length = length $3;
				#print "($1) ($2) ($insert_length) ($4) ($5)\n";
			}else{
				die "Unknown pattern flag, please contact Sur";
			}
			#$N++;
			
			# Primer filter is always used so we put the rest of the filters there;
			# Remove N
			unless ($keep_N || $insert =~ /^[ACGT]+$/){
				$N_seqs++;
				next;
			}

			unless($insert_length >= $min_length){
				$short_seqs++;
				next;
			}

			if($max_length && $insert_length >= $max_length){
				$long_seqs++;
				next;
			}
			
			# Print output
			print OUT ">$id\n$insert\n";
			$written_seqs++;
		}
	}
	#print "$N out of $i sequences matched the following regular expression:\n===$regexp===\n";
	print "Quality control finished...\n";
	print "Total sequences processed = $i\n";
	print "Sequences with unknown barcode = $unknown_barcode\n";
	print "Sequences matching regular expression = $match_regexp\n";
	print "Sequences with ambiguous bases = $N_seqs\n";
	print "Sequences shorter than ${min_length}bp = $short_seqs\n";
	print "Sequences longer than ${max_length}bp = $long_seqs\n";
	print "Sequences written to outfile = $written_seqs\n";
}

sub process_read_id{
	# Takes a Miseq read ID, and uses the barcode and the current counter to generate
	# a new ID of the form <sample ID>_<read number>
	my ($read_id,$map_ref,$counter_ref) = @_;
	my ($barcode,$new_id,$read_num,$name);
	if($read_id =~ /#([ACGT]+)$/){
		$barcode = $1;
	}else{
		die "The read id of the followign read does not appear to have a barcode ($read_id)";
	}

	if(exists($map_ref->{$barcode})){
		$read_num = ++$counter_ref->{$barcode};
		$new_id = "$map_ref->{$barcode}_$read_num";
	}else{
		print "Barcode ($barcode) not found on mapping file\n";
		$new_id = '';
	}
	
	return $new_id;
}

sub get_stats{
	# runs script that calculates Qscore statistics from a fastq file.
	my ($fastq_file,$outfile,$bindir) = @_;
	my $bin = "$bindir/fastq_mean_qscore.pl";
	die "Can't find executable script for calculating quality stats ($bin)" unless -f $bin;
	my $command = "$bin $fastq_file $outfile";
	execute($command);
	return;
}

sub plot_diagnostics{
	# Calls the R scprit tha makes the diagnostic plots for a MiSeq run.
	my ($params_file,$bindir) = @_;
	my $bin = "$bindir/plot_run_diagnostics.r";
	die "Can't find script for plotting run diagnostics ($bin)" unless -f $bin;
	my $command = "R CMD BATCH --no-save \"--args $params_file\" $bin";
	execute($command);
	return;
}

sub execute{
	# Executes system call with a passed command.
	my ($command) = @_;
	print "Executing:\n\t>$command\n";
	system($command);
	return;
}

sub reverse_complement{
	# Reverse complements a DNA sequence.
	my ($seq) = @_;
	my $rc_seq = reverse($seq);
	$rc_seq =~ tr/ACGT][/TGCA[]/;
	return $rc_seq;
}

sub make_regexp{
	my ($Fprimer,$Rprimer) = @_;
	my ($Fprimer_seq,$Rprimer_seq,$regexp1,$regexp2,$pattern);
	my $regexp = $regexp1 = $regexp2 = '';
	$pattern = 0;
	
	# Forward primer
	if ($Fprimer){
		$Fprimer_seq = get_primer($Fprimer);
		$regexp1 = "^(.*)($Fprimer_seq)";
		$pattern += 1;
	}else{
		$Fprimer_seq = ''
	}

	# Reverse primer
	if ($Rprimer){
		$Rprimer_seq = get_primer($Rprimer);
		$Rprimer_seq = reverse_complement($Rprimer_seq);
		$regexp2 = "($Rprimer_seq)(.*)\$";
		$pattern += 2;
	}else{
		$Rprimer_seq = '';
	}

	# Make regexp
	$regexp = "$regexp1(.+)$regexp2";
	return $regexp,$pattern;
}

sub get_primer{
	my ($primer_id) = @_;
	my $seq;
	if($primer_id eq '515F'){
		$seq = 'GTGCCAGC[CA]GCCGCGGTAA';
	}elsif($primer_id eq '806R'){
		$seq = 'GGACTAC[ACT][ACG]GGGT[AT]TCTAAT';
	}elsif($primer_id eq '1392R'){
		$seq = 'ACGGGCGGTGTGT[AG]C';
	}elsif($primer_id eq '804F'){
		$seq = 'ATTAGATACCC[AGT][AG]GTAGT';
	}elsif($primer_id eq '926F'){
		$seq = 'AAACT[CT]AAA[GT]GAATTGACGG';
	}elsif($primer_id eq '1114F'){
		$seq = 'GCAACGAGCGCAACCC';
	}elsif($primer_id eq 'ITS_9'){
		$seq = '[ACGT][ACGT][ACGT][ACGT][ACGT]GAACGCAGC[AG]AA[TAG][TAG]G[CT]GA';
	}elsif($primer_id eq 'ITS_4'){
		$seq = 'TCCTCCGCTTATTGATATGC';
	
	}else{
		die "Primer $primer_id not supported";
	}
	return $seq;
}

__END__

=head1 NAME

process_jgi.pl - Process a rDNA MiSeq run from JGI.

=head1 SYNOPSIS

=head1 USAGE

In kure it requires specific versions of Perl and R. Load them with:

	module load perl/5.8.9
	module load r/3.0.1

Starting from a a JGI file without phiX and contaminations, and after identification of barcodes use:

process_jgi.pl -infile I<<ncontam_nphix.fastq>>

=head1 OPTIONS

=over 8

=item -bindir, -b

Directory that contains the scripts for: i) calculating the Q-score statistics, ii) plotting the
run diagnostics, and iii) merges the reads.

Default='/proj/dangl_lab/sur/immune_code/'

=item -diagdir, -d

Name of the directory where the diagnostics plots will be saved. It will be created if it does not
exists.

Default='./diagnostics/'

=item -Fprimer, -F

Forward primer ID. If not passed the script will assume that the forward primer
is not present in the sequences. Currently supported primers: 515F, 804F,
806R, 926F, 1114F, 1392R, ITS_9, ITS_4.

Currenlty at least one primer must be provided, so either the -Fprimer or the -Rprimer
option must be well defined.

=item -help, -h

Prints help.

=item -infile, -i

Input fastq file. This file is required unless the -QC_only option is passed along. Should contain
every read after removing phiX and contamination. The read ID shoudl end with the pattern 
'#I<<Barcode sequence>>/I<<1 or 2>>', where the 1 represents teh forward read and 2 the reverse reads.
Reads should be arranged in a way that every forward read is followed by the reverse read of the same
cluster, the program will die otherwise.

This type of file is provided by JGI normally under the name ncontam_nphix.fastq, but no default
name has been implemented.

=item -keep_N, -k

If passed reads that have ambiguous bases (anything but ACGT), after merging would not be discarded.

The default behavior is to discard those reads and it is not recommended to use this parameter if the
intent is to use the sequences for OTU building, unless there is another quality control step implemented
later.

=item -mapping_file, -m

File name of the mapping file that matches barcodes to sample names. It should be a tab delimited file
where the first column gives the sample names and the second the barcode sequence.

=item -min_length, -l

Minimun length of the amplicon for sequences to be kept. Setting it to zero disables it.

DEFAULT=0

=item -max_length, -L

Maximum length of the amplicon for sequences to be kept. Setting it to zero disables it.

DEFAULT=0

=item -QC_only

If passed this option only performs the quality control, it won't do any preprocessing (trimming and merging).
There must be a file called I<out.extendedFrags.fastq> in the directory where the script is called, otherwise
the program will die.

=item -Rprimer, -R

Reverse primer ID. If not passed the script will assume that the forward primer
is not present in the sequences. Currently supported primers: 515F, 804F,
806R, 926F, 1114F, 1392R, ITS_9, ITS_4.

Currenlty at least one primer must be provided, so either the -Fprimer or the -Rprimer
option must be well defined.

=item -trim, -t

Length to which the sequences must be trimmed before merging. This is probably no neccessary since flash works
well with full length sequences but it is included for consistency with JGI.

Default='165'

=back

=head1 DESCRIPTION

=head1 AUTHOR

Sur from the Dangl Lab

=cut
