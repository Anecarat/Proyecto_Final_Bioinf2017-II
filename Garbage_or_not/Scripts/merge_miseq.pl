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

# Script that takes a ncontam_nphix.fastq file and process it until
# it merges reads.
# USAGE:
#	$ module load perl
#	$ merge_miseq_jgi.pl <nfile> <trim_length>

use warnings;
use strict;

print "====Welcome to merge_miseq_jgi.pl===========\n";

my ($infile,$trim_length) = @ARGV;
my ($command);

# parameters
my $split_bin = '/proj/dangl_lab/bin/fastq_split_pe_reads.pl';
my $trim_bin = '/proj/dangl_lab/bin/fastq_trim.pl';
my $flash_bin = "/proj/dangl_lab/bin/flash";

# Split reads
$command = "$split_bin $infile ncontam_nphix";
print "Executing\n\t>$command\n";
system($command);

# trim forward
$command = "$trim_bin --input ncontam_nphix.forward.fastq --output ncontam_nphix.forward.trim.fastq --first 1 --last $trim_length";
print "Executing\n\t>$command\n";
system($command);

# trim reverse
$command = "$trim_bin --input ncontam_nphix.reverse.fastq --output ncontam_nphix.reverse.trim.fastq --first 1 --last $trim_length";
print "Executing\n\t>$command\n";
system($command);

# Merge reads
$command = "$flash_bin ncontam_nphix.forward.trim.fastq ncontam_nphix.reverse.fastq -m 30 -M $trim_length -x 0.25 -r $trim_length -f 282 -s 20";
print "Executing\n\t>$command\n";
system($command);

print "===================GOOD BYE=================\n";
