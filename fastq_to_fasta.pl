#!/usr/bin/perl -w
use strict;

my $fastq_filepath = $ARGV[0];

open INFILE, "<", $fastq_filepath or die "ERROR: Could not open fastq file\n  $fastq_filepath\n$!\n";
my @buffer = ();    my @output_reads = ();
my $read_count = 0; my $outfile = ();
while (my $line = <INFILE>) {
    #chomp $line;
    # Keep track of previous lines; that will help us keep pace. (Store them in an array, and pop n' push them down).
    # This is a kind of sliding window approach. When the window can be determined to be in the right place, dump the buffer to the output file.
    push @buffer, $line;
    if (@buffer >= 4) {
        my $header = $buffer[0];
        my $seq = $buffer[1];
        $header =~ s/^\@/\>/;
        #print OUT "$header\n$seq\n";
        push @output_reads, "$header"."$seq";
        $read_count ++;
        @buffer = ();
    }
}
close INFILE;

foreach my $read (@output_reads) { print "$read"; }
    