#!/usr/bin/perl -w
use strict;
use File::Basename;

my $file = $ARGV[0];
my $check_n_characters = $ARGV[1];

my $lines_per_read = 4;
my ($basename, $parentdir, $extension) = fileparse($file, qr/\.[^.]*$/);
if ($extension =~ /\.fasta$/) { $lines_per_read = 2; }

open INFILE, "<", $file or die "ERROR: Could not open file\n  $file\n$!\n";
my @buffer = ();    my %seen_reads = ();    my @output_reads = ();
while (my $line = <INFILE>) {
    chomp $line;
    push @buffer, $line;
    if (@buffer >= $lines_per_read) {
        # Check if a hash exists for the sequence line
        # If not, write it, and set a value in the corresponding hash.
        my $seq = $buffer[1];
        my $checkseq = $seq;
        # In cases where we know the read quality drops towards the end of the reads, we may want to work out PCR duplicates
        # based on only the first n reads. There is a parameter in the config file for that. (Leave it blank to check everything).
        if ($check_n_characters) { $checkseq = substr($seq, 0, $check_n_characters); }
        
        if (!$seen_reads{$checkseq}) {
            push @output_reads, @buffer;
            $seen_reads{$checkseq} = 1;
        }
        @buffer = ();
    }
}
close INFILE;

foreach my $l (@output_reads) { print "$l\n"; }
