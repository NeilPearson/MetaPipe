#!/usr/bin/perl -w
use strict;
use Getopt::Long;

use List::Util qw(shuffle);
use File::Basename;
use lib "/usr/users/ga002/pearsonn/Scripts/airseq";
use MetaPipe;
my $funcs = MetaPipe->new();
use List::Compare;
use Text::Wrap;

$| ++;
# This is designed to call all the various programs that make up the MetaPipe pipeline.
# General goal is an easily comprehensible and extensible script that handles all the messy stuff so we don't have to.

# Many of the steps in this pipeline also take a considerable amount of time. That raises two problems:
#  -If repeating steps (particularly while developing and debugging, but plausibly at other times), we should be able to check for existing output before commencing a step.
#  -We should be able to wait for one step to complete before launching the next. This requires a means of understanding bsub jobs and waiting for the right ones to complete.

# Inputs
# Requested capability:
#  -Choose which aligners to run
#  -A good way of specifying input files
#  -A good way of handling subsampling
# May specify some inputs - the ones that control what the thing does, rather than what inputs it receives - via a config file.

# Processes/outputs (they follow on from each other)
# Start with output dir (/out). Presumably, this includes sample name.
# Begin                     Make /out/reads
# Dechunk and unzip         (original read location)
# Run FastQC first time;    make /out/fastqc
# Trim reads;               make /out/reads/trimmed (or whatever)
# Run FastQC again;         results in /out/fastqc
# Kontaminant filter;       make /out/reads/human_filtered
# Run flash;                make /out/reads/flashed
# Do subsampling;           put in /out/reads/subsampled/subsample_size.fastq
# For each of those:
#     Convert to fasta;     put in /out/reads/fasta_prepped/subsample_size.fastq
#     Run various aligners; put results in /out/aligners/aligner_name/subsample_size.fasta.txt
#     Run MEGAN;            put results in /out/megan/aligner_name/subsample_size.whatever

# Get all flags. Check that the right ones have been used.
my $incorrect_flags = 0;
foreach my $cla (@ARGV) {
    chomp $cla;
    if ($cla =~ /^-/) {
        if ($cla =~ /^-[a-zA-Z]/) {
            print "WARN: Possible error in command line argument $cla\n      Single dash detected where double expected\n";
        }
        unless ($cla =~ /^--(config|data|help|subsample|step|num_subsamples|unique_subsamples|num_chunks|exclude_all|output_prefix|overwrite|blastn|blastx|rapsearch|diamond|log_dir)/) {
            $incorrect_flags = 1;
            print "ERROR: Unknown option $cla\n";
        }
    }
}

my ($help, $subsample_start_size, $subsample_step, $number_of_subsamples, $unique_subsamples, $num_chunks, $exclude_all) = ();
my ($datadir, $configfile, $output_prefix, $overwrite, $log_path, $halt_at) = ();
# We have a choice of alignment tools! We may want to run one, or more than one. This should be user-specifiable at the command line.
my ($blastn, $blastx, $rapsearch, $diamond) = ();

GetOptions(
    'config=s'          => \$configfile,
    'data=s'            => \$datadir,       # Should be a directory containing a single sample's reads, with the name "Sample_[sample ID]".
    'help'              => \$help,
    'subsample:s'       => \$subsample_start_size,
    'step:s'            => \$subsample_step,
    'num_subsamples:s'  => \$number_of_subsamples,
    'unique_subsamples' => \$unique_subsamples,
    'num_chunks:s'      => \$num_chunks,
    'exclude_all'       => \$exclude_all,
    'output_prefix=s'   => \$output_prefix,
    'overwrite'         => \$overwrite,
    'blastn'            => \$blastn,
    'blastx'            => \$blastx,
    'rapsearch'         => \$rapsearch,
    'diamond'           => \$diamond,
    'log_dir:s'         => \$log_path,
    'halt_at:s'         => \$halt_at,
);

# The overwrite option is for times a user re-runs this pipeline - for example, using slightly different parameters.
# Normally, pipeline steps will check for existing output at each step before launching into what might be a lengthy bit of processing.
# But sometimes, a user might want to overwrite all that and start from the top. The overwrite parameter instructs the pipeline to do that. 

# Here, I would quite like to do a thorough check of my parameters concerning subsampling.
# However the number of reads is going to change as we trim, filter etc - and not in a particularly predctable way, either. Best I can hope for
# is to generate a kind of rule of thumb (from evidence of many runs) and throw a warning, but I can't do that yet.
# Instead, I have to simply check that the values are correctly formatted, and leave the rest until later.
# Remember that no subsample parameters may be supplied - in which case we just sample everything.
if ($subsample_start_size) {
    print "Checking subsampling inputs\n";
    # If this is set, it must be numeric and an integer.
    unless ($subsample_start_size =~ /^[0-9]+$/) {
        print "ERROR: Subsample initial size (--subsample) must be numeric!\n";
        $incorrect_flags = 1;
    }
    $subsample_start_size = int $subsample_start_size;
    
    # If this is set, the other two subsample inputs may also be set.
    # If one of those is set, both must be; but if neither is set, it's OK.
    # These must also be numeric.
    if ($subsample_step) {
        unless ($number_of_subsamples) {
            print "ERROR: Number of subsamples (--num_subsamples) must also be set if subsample step-up (--step) is set!\n";
            $incorrect_flags = 1;
        }
        unless ($subsample_step =~ /^[0-9]+$/) {
            print "ERROR: Subsample step-up size (--step) must be numeric!\n";
            $incorrect_flags = 1;
        }
        $subsample_step = int $subsample_step;
    }
    
    if ($number_of_subsamples) {
        unless ($subsample_step) {
            print "ERROR: Subsample step-up (--step) must also be set if Number of subsamples (--num_subsamples) is set!\n";
            $incorrect_flags = 1;
        }
        unless ($number_of_subsamples =~ /^[0-9]+$/) {
            print "ERROR: Subsample step-up size (--step) must be numeric!\n";
            $incorrect_flags = 1;
        }
        $number_of_subsamples = int $number_of_subsamples;
    }
}
else {
    # Just to be sure, blank the other parameters if subsample start size isn't supplied.
    print "No subsampling required\n";
    $subsample_step = ();
    $number_of_subsamples = ();
}


# Display help text if requested, or if incorrect flags detected
if (($help) || ($incorrect_flags == 1)) {
    die 
    "METAPIPE PIPELINE
A wrapper script for a set of tools used for the metagenomic analysis of sequencing runs.
The procedure is as follows:
--Measure and correct PCR duplication using NextClip
--Trim adapters from reads (Trimmomatic)
--Assess success of trimming (FastQC)
--Filter out contaminant reads (Kontaminant)
--Overlap paired-end reads (FLASH)
--Align reads to database, using any combination of
    BLASTN
    BLASTX
    RapSearch
    Diamond
--Metagenomic analysis and visualisation (MEGAN)

Input: the pipeline is capable of working with chunked and zipped reads, but these must be present in a single directory.
The directory path is supplied as input.
-----
Command line options
Flag                Description
--config            Path to config file
--data              Path to directory containing input data
--help              This help message
--subsample         Number of reads to subsample. If --step and --num_subsamples are supplied, this
                    specifies the size of the smallest subsample.
--step              Number of reads by which the subsample size increases on successive subsampling              
--num_subsamples    Number of subsamples to take (if greater than 1)
--unique_subsamples Reads used in each subsample should be unique (i.e., used in no other subsample)
--num_chunks        The number of sections a subsample will be re-divided into, in order to speed up alignment processing.
--exclude_all       Do not use the full set of reads in this analysis (requires --num_subsamples value)
--output_prefix     Path to output directory (will be created if does not exist)
--overwrite         Existing data in the output directory will be overwritten with that from new analyses,
                    rather than used as a starting point for continuation of the pipeline
--blastn            Use the BLASTN aligner
--blastx            Use the BLASTX aligner
--RapSearch         Use the RapSearch aligner
--Diamond           Use the Diamond aligner
--log_dir           Path in which logs will be stored
--halt_at           Stop the pipeline before launching the supplied step (options: trimming, filtering, flash, alignment, megan)

\n";
}

# Remove tailing slash from output prefix, if present
$output_prefix =~ s/\/$//; 
# Then check that it's a real directory
print "Check output directory\n";
$funcs->directory_check($output_prefix);

# Use the data directory ($datadir) to figure out the output prefix, if not set already
# Output prefix is a little tricky. I might supply a general one and have the script figure out a specific one for this sample. Or, I might supply a
# specific dir outright. There might be contexts in which I might want to do either. Figure it out.
# Let's try to get the sample ID from the data directory path, and go from there.
#my @sp = split /\//, $datadir;
#my $datadir_sampleid = ();
#if ($sp[-1] =~ /Sample_/) {
#    # If this level matches, add it to the output path
#    $datadir_sampleid = $sp[1];
#    # Add sample ID tag to output prefix - but only if it's not there already!
#    unless ($output_prefix =~ /$datadir_sampleid$/) { $output_prefix = "$output_prefix/$datadir_sampleid"; }
#    $funcs->directory_check($output_prefix);
#}
# NOTE: I want to avoid using this, since there is no real need for it.

# Check the log path; supply if not set
if (!$log_path) { $log_path = $funcs->directory_check("$output_prefix/logs"); }

# Inputs have been checked; I can now add them to an object.
# Config file also gets parsed when this call is made
print "OO-ise input parameters\n";
$funcs->assign_parameters(
                           $configfile,             $datadir,               $subsample_start_size,
                           $subsample_step,         $number_of_subsamples,  $unique_subsamples,
                           $num_chunks,             $exclude_all,           $output_prefix,
                           $overwrite,              $log_path,              $blastn,
                           $blastx,                 $rapsearch,             $diamond,
                           $halt_at);

# Let's work on the reasonable assumption that we're using data from a PAP directory, which follows the same organisational structure.
# A run's data is probably going to be contained within a single project directory, though there may be multiple multiplexed samples in there.
# To keep this simple, shall we assume that $datadir points to a single sample's directory, and therefore directly contains all the data we want?

# That also helps us set a default value for the output directory. But I'll get to that shortly.

# Merge input files - they may be chunked.
# I can copy in my dechunk'n'unzip code from the NextClip wrapper to do that.
$datadir =~ s/\/$//; # Remove tailing slash
print "Looking for input files in $datadir\n";
my $pattern = '*.fastq*';
my $readsfiles = $funcs->find_files($datadir, $pattern);
if ($readsfiles->[0] =~ /No such file or directory/) {
    die "ERROR: Cannot find read files in directory\n  $datadir\n";
}
print "Found ".@$readsfiles." input files:\n";
foreach my $infile (@$readsfiles) { print "$infile"; }

# Determine what needs to be done to get a single, unzipped file for reads 1 and 2 (and then do it).
# NOTE: Both of these hold the complete path from root.
# ANOTHER NOTE: Going to assume for the moment that we're doing paired end only.
print "Dechunk and unzip input files\n";
my $read1file = $funcs->dechunk_and_unzip(1, $readsfiles);
my $read2file = $funcs->dechunk_and_unzip(2, $readsfiles);
# For the moment, we can leave these unzipped/dechunked read files where they are. But later, we create several different sets -
# fasta format, flashed, filtered, subsampled, etc - that ought to go in the output directory somewhere. 

# Run NextClip, in order to measure the level of PCR duplicates, and to remove them if called for.
# Note that $read1file and $read2file continue to hold full path names.
my $nextclip_path = $funcs->directory_check("$output_prefix/nextclip");
print "Run NextClip on reads";
($read1file, $read2file) = $funcs->run_nextclip($read1file, $read2file, $nextclip_path);

# Run FastQC (first time); get relevant info out of it
# Run FastQC (and the other stuff) from in here, rather than via Richard's shell scripts. It's more flexible, not to mention easily maintainable, if I do it here.
# I can copy-paste bits to get the right commands as needed.
# Run it on read 1 and 2 independently?
# Now, here's a problem. I'm going to have to figure out how to read the FastQC data to get what I want. (Write a sub - I use it again later).
# Moreover, I'm going to need some sort of directory structure so I can organise this stuff
$funcs->directory_check("$output_prefix/fastqc");
my $fastqc_firstpass_path = $funcs->directory_check("$output_prefix/fastqc/untrimmed");


# Because I'm bsubbing jobs (to take advantage of our cluster's opportunities for parallelisation) I need to copy my waiting sub in.
# That gets a job ID and waits for it to complete before cintinuing with the script. Without it, we'll charge ahead to stuff wer're not ready for.
print "Run FastQC on unmodified input files\n";
my ($fastqcfileR1, $fastqcfileR2) = $funcs->run_fastqc($read1file, $read2file, $fastqc_firstpass_path);

# Look at those FastQC results; determine appropriate action
# (Amounts to throwing a minor warning, at most, at this stage)
my $trimming_done = ();
$funcs->examine_fastqc_results($fastqcfileR1, $trimming_done);
$funcs->examine_fastqc_results($fastqcfileR2, $trimming_done);

#ÊTrim the reads 
# Use trimmomatic. It's a bit of a complex call by the look of it, but I can handle it.
my $threads = 4;
$funcs->directory_check("$output_prefix/reads");
my $trimming_data_dir = $funcs->directory_check("$output_prefix/reads/trimmed_fastq");
print "Run trimming on input data\n";
my ($read1file_trimmed, $read2file_trimmed) = $funcs->run_trimming($read1file, $read2file, $trimming_data_dir);

# Re-run FastQC, to check adapters have been removed.
# (That is, if trimming has been specified!)
if ($funcs->{config}{run_trimming} =~ /yes/) {
    my $fastqc_secondpass_path = $funcs->directory_check("$output_prefix/fastqc/trimmed");
    print "Run FastQC again on trimmed data\n";
    ($fastqcfileR1, $fastqcfileR2) = $funcs->run_fastqc($read1file_trimmed, $read2file_trimmed, $fastqc_secondpass_path);
    
    # Look at those second FastQC results; determine appropriate action
    $trimming_done = 1;
    $funcs->examine_fastqc_results($fastqcfileR1, $trimming_done);
    $funcs->examine_fastqc_results($fastqcfileR2, $trimming_done);
}

# Run kontaminant - remove human reads
my $filtering_data_dir = $funcs->directory_check("$output_prefix/reads/human_filtered_fastq");
print "Run kontaminant on trimmed reads\n";
my ($read1file_filtered, $read2file_filtered) = $funcs->run_kontaminant($read1file_trimmed, $read2file_trimmed, $filtering_data_dir, $funcs->{config}{human_reference});

# Run FLASH - join overlapping paired-end reads into single reads
my $flashed_data_dir = $funcs->directory_check("$output_prefix/reads/flashed_fastq");
my $readsfile_flashed = "flashed_reads"; # Note: Not quite the complete filename; no extension. FLASH adds one.
print "Run FLASH on filtered reads\n";
$readsfile_flashed = $funcs->run_flash($read1file_filtered, $read2file_filtered, $flashed_data_dir, $readsfile_flashed);
print "Flashed reads file:\n$readsfile_flashed\n";

# Subsampling
# This is where it starts to get a little more complex.
# Inputs can specify nothing (sample everything), a subsample of size x (one subsample of x reads), or a starting subsample size, a step, and a number of subsamples to do.
# In the latter case, we'll step through subsamples of that size and process each in turn. The pipeline needs to check that there are sufficient runs to follow the suggested plan. (That is, as well as the normal analysis of all reads)
# One more thing: reads used in one subsample then can't be used in others. 
# (We want to do that in order to assess the level of species diversity - it'll tend towards an asymptote as sample size increases)
my $aligners = $funcs->list_aligners();

my $subsample_dir = $funcs->directory_check("$output_prefix/reads/subsampled");
print "Make subsamples of the data\n";
my $subsamples = $funcs->make_subsamples($readsfile_flashed, $subsample_dir);
# Chunk those subsets up too
foreach my $subsample_file (@$subsamples) {
    my $chunked_files = $funcs->rechunk($subsample_file);
}

my $megan_dir = $funcs->directory_check("$output_prefix/MEGAN");

my $subset_sizes = $funcs->list_subset_sizes();
foreach my $subset (@$subset_sizes) {
    print "-----\nSubset $subset\n";
    
    # Chunked reads: output_prefix/reads/subsampled/[number of reads]/[chunk number].fastq
    # Need results files out, separated by aligner used to produce them.
    # Also need list of jobs
    # Set a base aligner output directory too
    # Make a general output directory
    $funcs->directory_check("$output_prefix/alignments");
    $funcs->{config}{alignment_base_dir} = $funcs->directory_check("$output_prefix/alignments/$subset");
    
    print "Convert to FASTA\n";
    my (@fasta_files, @jobs) = ();
    foreach my $chunk (1..$funcs->{param}{num_chunks}) {
        print "    chunk $chunk\n";
        my $chunkdir = "$subsample_dir/$subset";
        my $chunk_o_reads = "$chunkdir/$chunk.fastq";
        
        my ($fasta_file, $conversion_jobs) = $funcs->convert_fastq_to_fasta($chunk_o_reads, $chunkdir);
        push @fasta_files, $fasta_file;
        if ($conversion_jobs) {
            push @jobs, @$conversion_jobs;
        }
    }
    # Wait for these all to finish before going further
    if (@jobs) { my $done = $funcs->done_when_its_done(\@jobs); }
    
    my %alignment_results = ();
    my %alignment_jobs = ();
    print "Preparing alignments for...\n";
    foreach my $chunk (1..$funcs->{param}{num_chunks}) {
        print "    chunk $chunk\n";
        my $chunkdir = "$subsample_dir/$subset";
        my $chunk_o_reads = "$chunkdir/$chunk.fasta";
        
        #print "Run selected alignments!\n";
        my ($results, $jobs) = $funcs->run_alignments($chunk_o_reads);
        
        # $alignment_results is an array reference, containing file paths.
        # $alignment_jobs is a hash of arrays reference, organised by aligner.
        # Sort these out, so we can get all result files for a given aligner.
        foreach my $i (1..@$results) {
            $i--;
            my $file = $results->[$i];
            
            my $aligner_used = $funcs->get_aligner_used($file);
            push @{$alignment_results{$aligner_used}}, $file;
            
        }
        
        foreach my $aligner (@$aligners) {
            if ($jobs->{$aligner}) { push @{$alignment_jobs{$aligner}}, @{$jobs->{$aligner}}; }
        }
    }
    
    # Prepare to run MEGAN on each aligner's output
    # The jobs will be submitted, but won't commence until the respective aligner jobs finish.
    my $megan_sampledir = $funcs->directory_check("$megan_dir/$subset");
    foreach my $aligner (@$aligners) {
        my $megan_alignerdir = $funcs->directory_check("$megan_sampledir/$aligner");
        my $results = $alignment_results{$aligner};
        my $jobs = $alignment_jobs{$aligner};
        my $log = "$megan_sampledir/MEGAN_log.lsf";
        
        $funcs->run_megan($results, \@fasta_files, $megan_alignerdir, $log, $jobs);
    }
}

# Do something to demonstrate conclusively that the pipeline is complete.
`echo date > $log_path/complete.txt`;