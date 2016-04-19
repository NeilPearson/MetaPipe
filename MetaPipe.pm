package MetaPipe;
use strict;
use warnings;
# This is the directory where config files, software license files etc. all go.
# Adjust as necessary.
use lib "/usr/users/ga002/pearsonn/Scripts/metapipe";
use List::Compare;
use List::Util qw(shuffle min max);
use File::Basename;
use Cwd;

$|++;
my %jobs_check_buffer = ();

sub new {
    my $class = shift;
    my $self = bless {}, $class;
    return $self;
}

sub assign_parameters {
    # Adds the user-supplied parameters to the object, so I can get at them anywhere.
    my $self = shift;
    
    $self->{param}{config} = shift;
    $self->{param}{data} = shift;
    $self->{param}{sub_start_size} = shift;
    
    $self->{param}{sub_step} = shift;
    $self->{param}{sub_num} = shift;
    $self->{param}{sub_unique} = shift;
    
    $self->{param}{num_chunks} = shift;
    $self->{param}{exclude_all} = shift;
    $self->{param}{output_prefix} = shift;
    
    $self->{param}{overwrite} = shift;
    $self->{param}{log_path} = shift;
    $self->{param}{blastn} = shift;
    
    $self->{param}{blastx} = shift;
    $self->{param}{rapsearch} = shift;
    $self->{param}{diamond} = shift;
    
    $self->{param}{halt_at} = shift;
    
    # Check for non-optional inputs here.
    if (!$self->{param}{config}) { die "ERROR: Config file input missing. Use the --config flag to supply one.\n"; }
    $self->does_file_exist($self->{param}{config});
    if (!$self->{param}{data}) { die "ERROR: Data file/directory input missing. Use the --data flag to supply one.\n"; }
    if (!$self->{param}{output_prefix}) { die "ERROR: Output prefix input missing. Use the --output_prefix flag to supply one.\n"; }
    
    if ($self->{param}{num_chunks}) {
        if ($self->{param}{num_chunks} =~ /^[0-9]+$/) {
            $self->{param}{num_chunks} = int $self->{param}{num_chunks};
        }
        else {
            die "ERROR: Parameter --num_chunks must be numeric!\n";
        }
    }
    else {
        $self->{param}{num_chunks} = 1;
    }
    
    # Lowercase and strip non-alphnum characters from this halt_at input
    if ($self->{param}{halt_at}) {
        chomp $self->{param}{halt_at};
        $self->{param}{halt_at} = lc $self->{param}{halt_at};
        $self->{param}{halt_at} =~ s/[^a-z0-9]//g;
    }
    
    $self->read_config_file;
    
    # Do a check here: if no aligner parameters are supplied, fill them in based on the defaults supplied in the config file.
    if ((!$self->{param}{blastn}) && (!$self->{param}{blastx}) && (!$self->{param}{rapsearch}) && (!$self->{param}{diamond})) {
        if ($self->{config}{use_blastn}) {
            if ($self->{config}{use_blastn} eq 'yes') { $self->{param}{blastn} = 1; }
        }
        if ($self->{config}{use_blastx}) {
            if ($self->{config}{use_blastx} eq 'yes') { $self->{param}{blastx} = 1; }
        }
        if ($self->{config}{use_rapsearch}) {
            if ($self->{config}{use_rapsearch} eq 'yes') { $self->{param}{rapsearch} = 1; }
        }
        if ($self->{config}{use_diamond}) {
            if ($self->{config}{use_diamond} eq 'yes') { $self->{param}{diamond} = 1; }
        }
    }
    
    if ($self->{param}{num_chunks} == 1) {
        # If this is set to 1, it should be treated exactly the same as not being set at all. 
        $self->{param}{num_chunks} = ();
    }
}

sub submit_job {
    # TGAC is moving to a different job manager (slurm, from LSF); metapipe should be able to support either.
    # To facilitate that, I'll have all jobs passed through this function.
    # $cmd is the string that we want to pass to the system to execute.
    # $opts is a hash reference, holding the LSF/slurm options. If not set, I'll just use defaults from the config file. 
    my $self = shift;
    my $cmd = shift;
    my $opts = shift;
    my $scheduler = $self->{config}{jobsys};
    my $output_prefix = $self->{param}{output_prefix};
    
    if (!$opts->{threads}) { $opts->{threads} = 1 }                      # Some processes are only single-threaded, so should request only 1 thread. Assume 1 if unfilled.
    if (!$opts->{memory}) { $opts->{memory} = $self->{config}{memory} }  # Memory is more likely to require adjustment for particular purposes.
    if (!$opts->{queue}) { $opts->{queue} = $self->{config}{queue} }     # Queue should be left to the default setting in just about all cases.
    if (!$opts->{log}) { $opts->{log} = '/dev/null' }                    # Different logs will be specified in almost all cases.
    if (!$opts->{hosts}) { $opts->{hosts} = 1 }                          # Only a single host will be requested in almost all cases.
    
    # Some extra things can be present in $opts, though I won't set them to any defaults because they're optional:
    # Array of job IDs to set a job as dependent to
    # Array of things to source
    # Job name
    
    if ($scheduler eq 'lsf') {
        my $str = ();
        # First, source relevant software.
        if ($opts->{source}) {
            foreach my $src (@{$opts->{source}}) { $str .= "source $src; " }
        }
        # Start command - add threads, queue, log
        $str .= "bsub -n ".$opts->{threads}." -q ".$opts->{queue}." -oo ".$opts->{log}." ";
        # Add job name, if one has been set
        if ($opts->{jobname}) { $str .= "-J ".$opts->{jobname}." " }
        # Add memory
        # Note that we also add a restriction here, in which we set the number of hosts. It's adjustable, but usually easier to leave it alone.
        $str .= '-R "rusage[mem='.$opts->{memory}.'] span[hosts='.$opts->{hosts}.']" ';
        # Add dependencies
        if ($opts->{depend}) {
            my $dependency_string = ();
            if (ref $opts->{depend} eq 'ARRAY') {
                my @dep = ();
                foreach my $job (@{$opts->{depend}}) { push @dep, "ended($job)"; }
                $str .= '-w "'.join " && ", @dep.'" ';
            }
            else { $str .= '-w "ended('.$opts->{depend}.')" '; }
        }
        # Add command string (note, it's single-quoted; I don't have to worry about interpolating into $cmd here, and it will prevent interference from special characters etc.)
        $str .= "'".$cmd."'";
        
        # Submit and get the job ID
        my $r1jobs = `$str`;
        my $jobs = $self->extract_list_of_jobs($r1jobs);
        
        # Return the full command constructed, and the job ID.
        return ($str, $jobs);
    }
    elsif ($scheduler eq 'slurm') {
        # Slurm jobs can be set up in much the same way as LSF jobs, but there are enough differences between the two to justify a separate implementation here.
        my $str = ();
        # First, source relevant software.
        if ($opts->{source}) {
            foreach my $src (@{$opts->{source}}) { $str .= "source $src; " }
        }
        # Start command - add threads, queue, log
        $str .= "sbatch -c ".$opts->{threads}." -p ".$opts->{queue}." -o ".$opts->{log}." -e ".$opts->{log}." "; # Note that output and error go to the same log file
        # Add job name, if one has been set
        if ($opts->{jobname}) { $str .= "-J ".$opts->{jobname}." " }
        # Add memory
        $str .= '--mem '.$opts->{memory}." ";
        # Add dependencies
        if ($opts->{depend}) {
            my $dependency_string = ();
            if (ref $opts->{depend} eq 'ARRAY') {
                $str .= '--dependency afterok:'.join ":", @{$opts->{depend}}.' ';
            }
            else { $str .= '--dependency afterok:'.$opts->{depend}.' '; }
        }
        # Add command string (note, it's single-quoted; I don't have to worry about interpolating into $cmd here, and it will prevent interference from special characters etc.)
        $str .= "--wrap '".$cmd."'";
        
        # Submit and get the job ID
        my $r1jobs = `$str`;
        my $jobs = $self->extract_list_of_jobs($r1jobs);
        
        # Return the full command constructed, and the job ID.
        return ($str, $jobs);
    }
    else {
        die "ERROR: Unrecognised job scheduler '$scheduler'\n";
    }
}

sub read_config_file {
    # Reads the config file and adds the parameters within to the object, so I can get at them anywhere.
    my $self = shift;
    # These can all go in $self->{config}{name}
    # It looks, reasonably enough, as though there are distinct databases for each alignment tool. Make sure they're distinguished.
    # Also include version numbers of various software! This makes switching versions that bit easier.
    my $configfile = $self->{param}{config};
    $self->does_file_exist($configfile);
    
    open(CONFIG, "<", $configfile) or die "ERROR: Cannot open config file\n  $configfile\n$!\n";
    my @config = <CONFIG>;
    close CONFIG;
    
    #print "Found the following config values:\n";
    LINE: foreach my $line (@config) {
        chomp $line;
        $line =~ s/ +/ /g;
        # This removes any comment strings
        $line =~ s/\#.*//g;
        
        # Skip over empty/whitespace only lines
        if (length $line == 0) { next LINE; }
        if ($line =~ /^[\s]+$/)  { next LINE; }
        
        my ($key, $value) = split /\s/, $line;
        $self->{config}{$key} = $value;
        #print "$key:\t$value\n";
    }
}

sub find_files {
    # Returns a list of files in a given directory matching a pattern
    my $self = shift;
    my $dir = shift;
    my $pattern = shift;
    
    my @files = `ls -l $dir/$pattern | awk '{print \$9}'`;
    my @actual_files = ();
    foreach my $file (@files) {
        $file =~ s/\n//g;
        if ($file !~ /No such file or directory/) {
            push @actual_files, $file;
        }
    }
    return \@actual_files;
}

sub does_file_exist {
    # Checks if a file exists as it ought to
    my $self = shift;
    my $file = shift;
    
    if (!$file) {
        die "ERROR: Can't even find a filename. Come on.\n";
    }
    
    unless (-e $file) {
        die "ERROR: File\n  $file\ndoes not exist.\n";
    }
    
    # Also check size of file; die if it's zero bytes.
    if (-z $file) {
        die "ERROR: File\n  $file\nis empty (0 bytes).\n";
    }
}

sub halt {
    my $self = shift;
    my $step = shift;
    my $halt_at = $self->{param}{halt_at};
    
    if ($halt_at) {
        if ($halt_at eq $step) { die "Pipeline halted at $step\n"; }
    }
}

sub directory_check {
    # Checks if a directory exists; tries to make it if not; halts the program if it can't make it.
    my $self = shift;
    my $dir = $_[0];
    if (!-d $dir) {
        `mkdir $dir`;
        if (!-d $dir) {
            # If it still doesn't exist, throw an error
            die "ERROR: Cannot create output directory\n  $dir\n";
        }
    }
    return $dir;
}

sub copy_raw_reads {
    # Copies the supplied reads files to a specific location within the data directory. Protects the master copy from accidental deletion/alteration, and
    # gets around permissions issues.
    my $self = shift;
    my $files = shift;
    my $output_prefix = $self->{param}{output_prefix};
    my $log_path = $self->{param}{log_path};
    my $queue = $self->{config}{queue};
    my $overwrite = $self->{param}{overwrite};
    
    my $file_output_prefix = $self->directory_check("$output_prefix/reads/raw");
    my $newfiles = ();
    my @jobs = ();
    print "Copying raw reads to\n  $file_output_prefix\n";
    foreach my $file (@$files) {
        # Check if it exists
        $self->does_file_exist($file);
        # Check if it's in the output location already
        my $bn = basename($file);
        if ((-e "$file_output_prefix/$bn") && (!$overwrite)) { print "File\n $file_output_prefix/$bn\nexists; skipping copy.\n"; }
        
        my $cmd = "cp $file $file_output_prefix/$bn";
        my ($outcmd, $jobs) = $self->submit_job($cmd);
        push (@jobs, @$jobs);
        push @$newfiles, "$file_output_prefix/$bn";
    }
    $self->done_when_its_done(\@jobs);
    foreach my $file (@$newfiles) { $self->does_file_exist($file); } 
    return $newfiles;
}

sub dechunk_and_unzip {
    my $self = shift;
    
    my $files = $_[0];
    my $read = $_[1];
    my $readfile = ();
    
    # OK, time to use some recursion in anger again.
    # I want to pass all my reads files in here, whether they're paired end or not. I want this to be able to figure it out.
    # As such, I won't pass in read numbers initially (at least, in some cases). If no filter input ($read) is set, I'll check the
    # config option that tells us if this is a paired end or single end run, and recursively call this sub with the correct filters.
    if ($read) {
        # Be aware that the correct file may already exist, amongst a set of chunked and/or compressed files.
        # Try to pick it out, if possible, before going any further.
        my @readfiles = ();
        FILES: foreach my $file (@$files) {
            chomp $file;
            $self->does_file_exist($file);
            if ($file =~ /_R$read\.fastq/) { $readfile = $file; last FILES;  }
            elsif ($file =~ /_R$read\_[0-9]+\.fastq/) { push @readfiles, $file; }
        }
        
        # At this point, the correct file may already have been positively ID'd - so we only need carry on with this other stuff if it hasn't.
        if (!$readfile) {
            if (@readfiles == 0) {
                die "ERROR: Cannot find Read $read data\n";
            }
            elsif (@readfiles == 1) {
                $readfile = $readfiles[0];
                chomp $readfile;
                # Is it compressed?
                if ($readfile =~ /fastq.gz/) {
                    $readfile =~ s/fastq.gz/fastq/;
                    # Unzip it. See if I can keep the compressed version (use -c)
                    print "Unzipping reads to file\n  $readfile\n";
                    my $cmd = "gunzip -c $readfile.gz 1> $readfile";
                    my ($outcmd, $jobs) = $self->submit_job($cmd);
                    $self->done_when_its_done($jobs);
                }
            }
            else {
                # Unlike with a single file, an output filename has to be assembled here. It should match the rest, but without the chunk number.
                # To get that, split the first thing in the list on _, remove numbers from the last part, and re-join.
                # For real thoroughness, I ought to check all files, but this will work well enough for now. 
                $readfile = $readfiles[0];
                chomp $readfile;
                my @parts = split /_/, $readfile;
                my $end = pop @parts;   # The end of the reads name includes the extension as well as the chunk number
                $end =~ s/[0-9]+//g;
                $parts[-1] .= $end;
                $readfile = join '_', @parts;
                $readfile =~ s/fastq.gz/fastq/;
                
                # One thing to double-check: does the file we're after already exist?
                # There's no point re-extracting data if it already exists. 
                unless (-e $readfile) {
                    # Check if these files are compressed or not; it determines what we do with them.
                    # If decompressing, we should only put compressed files in the list.
                    my $compressed = 0;
                    my @compressed_files = ();
                    foreach my $file (@readfiles) {
                        if ($file =~ /fastq.gz/) { $compressed = 1; push @compressed_files, $file; }
                    }
                    
                    # If compressed, feed all files (restricted to compressed files) to gunzip -c; if not, use cat.
                    my $call = ();
                    if ($compressed == 1) { $call = "gunzip -c "; @readfiles = @compressed_files; }
                    else                  { $call = "cat "; }
                    
                    foreach my $file (@readfiles) {
                        chomp $file;
                        $call .= "$file ";
                    }
                    $call .= "1> $readfile";
                    print "Reads file appears to be chunked. Concatenation into a single file requested.\nCommand:\n  $call\n";
                    `$call`;
                }
            }
        }
        else {
            # May have correctly identified the single read file, but it may be compressed. Deal with that here.
            print "Found single valid reads file\n  $readfile\n";
            if ($readfile =~ /fastq.gz/) {
                $readfile =~ s/fastq.gz/fastq/;
                # Unzip it. See if I can keep the compressed version (use -c)
                print "Unzipping reads to file\n  $readfile\n";
                my $cmd = "gunzip -c $readfile.gz 1> $readfile";
                my ($outcmd, $jobs) = $self->submit_job($cmd);
                $self->done_when_its_done($jobs);
            }
        }
        my $outfiles = ();
        push @$outfiles, $readfile;
        foreach my $outfile (@$outfiles) { $self->does_file_exist($outfile); }
        return $outfiles;
    }
    else {
        # Do the recursive thing I describe above
        my $r1files = $self->dechunk_and_unzip($files, 1);
        if ($self->{config}{reads_type} eq 'paired end') {
            my $r2files = $self->dechunk_and_unzip($files, 2);
            push @$r1files, @$r2files;
        }
        foreach my $outfile (@$r1files) { $self->does_file_exist($outfile); }
        return $r1files;
    }
}

sub list_subset_sizes {
    # Returns a list of the subset sizes we use (as dictated by inputs).
    my $self = shift;
    return $self->{subset_sizes};
}

sub count_fastq_reads {
    # Counts the number of reads in a fastq file.
    # Simple, really - counts lines, divides by 4.
    # Bear in mind that for large files - and we always have large files - this can take quite a long time.
    my $self = shift;
    my $fastq_file = shift;
    
    $self->does_file_exist($fastq_file);
    
    my $count = `wc -l $fastq_file`;
    
    # Error checking
    
    $count = $count / 4;
    return $count;
}

sub get_fastq_read_ids {
    # Returns a big array of all the read IDs in a FastQ file. 
    my $self = shift;
    my $file = shift;
    
    $self->does_file_exist($file);
    my @ids = ();
    
    open INFILE, "<", $file or die "ERROR: Could not open fastq file\n  $file\n$!\n";
    my @buffer = ();
    my $read_count = 0; my $outfile = ();
    while (my $line = <INFILE>) {
        chomp $line;
        # Keep track of previous lines; that will help us keep pace. (Store them in an array, and pop n' push them down).
        # This is a kind of sliding window approach. When the window can be determined to be in the right place, dump the buffer to the output file.
        push @buffer, $line;
        if (@buffer == 4) {
            my $header = $buffer[0];
            push @ids, $header;
            $read_count ++;
            @buffer = ();
        }
    }
    close INFILE;
    close OUT;
    
    return \@ids;
}

sub make_subsamples {
    # Try again, but keep it simpler. Get rid of this random stuff and just pick the first n reads, since they're not in any particular order anyway.
    # Make a consistent output structure so I always know what's happening
    # Chunking will remain a similar process, though it becomes a lot easier to see how it will work after looking at this.
    # Pro mode: make the file recognised as 'all', when requested, actually be the file that comes out of flash. Use a symlink?
    
    my $self = shift;
    my $start_size = $self->{param}{sub_start_size};
    my $step_size = $self->{param}{sub_step};
    if (!$step_size) { $step_size = 0; }
    my $num_subsamples = $self->{param}{sub_num};
    my $exclude_all = $self->{param}{exclude_all};
    my $overwrite = $self->{param}{overwrite};
    
    #my $sample_id = $self->{param}{sample_id};
    # Bear in mind that these last two inputs may not be supplied, or may be null.
    # In that case, make only one subsample.
    my $readfile = shift;
    my $outdir = shift;
    
    # Take care of some stuff to make sure we're getting an accurate, nonsense-free look at the number of reads in the input file.
    # Adjusting this lets this sub handle either fastq or fasta input. 
    #$self->remove_trailing_newlines($readfile);
    my $lines_per_read = 4;
    my ($basename, $parentdir, $extension) = fileparse($readfile, qr/\.[^.]*$/);
    if ($readfile =~ /\.fasta$/) { $lines_per_read = 2; }
    
    my %subsample_files = ();
    
    # Get number of reads in the file
    # Choose appropriate division
    my $number_of_reads = `wc -l $readfile`;
    my @sp = split /\s/, $number_of_reads;
    $number_of_reads = $sp[0];
    $number_of_reads /= $lines_per_read;
    # Store this figure; it's useful.
    $self->{param}{number_of_reads} = $number_of_reads;
   
    # Set up subsample numbers
    # Check that they don't overrun number of available reads; react appropriately if they do
    # Store the numbers in $self somewhere
    # Modify list_subset_sizes to simply recall that
    # Maybe try to also handle the case where no subsetting parameters are supplied at all.
    my %subsample_sizes = ();
    
    # Allow the option of excluding the 'all' setting - but only if a subsample start size is given.
    # If no subsample sizes are specified at all, we definitely want the 'all' option in there.
    if ($start_size) {
        if (!$exclude_all) { $subsample_sizes{all} = 1; }
    }
    
    my $subsize = $start_size;
    if ($num_subsamples) {
        # If a number of subsamples is set ($start_size will also be set), we want multiple subsamples, stepping up by $step_size each time.
        for (1..$num_subsamples) {
            if ($subsize < $number_of_reads) { $subsample_sizes{$subsize} = 1; }
            else {
                $subsample_sizes{all} = 1;
                print "WARN: Not enough reads to meet planned subsampling scheme!\n  Proposed sample:\t$subsize\n  Available:\t$number_of_reads\n";
            }
            $subsize += $step_size;
        }
    }
    elsif ($start_size) {
        # If number of subsamples is not set, we want a single subsample of $start_size. 
        if ($subsize < $number_of_reads) { $subsample_sizes{$subsize} = 1; }
        else {
            $subsample_sizes{all} = 1;
            print "WARN: Not enough reads to meet planned subsampling scheme!\n  Proposed sample:\t$subsize\n  Available:\t$number_of_reads\n";
        }
    }
    else {
        # If neither are set, we want all reads, regardless of size.
        $subsample_sizes{all} = 1;
    }
    
    my @subset_sizes = keys %subsample_sizes;
    $self->{subset_sizes} = \@subset_sizes;
    
    # OK, let's do this properly this time. No fancy rubbish - just pull n lines out at a time until we have the right number.
    my @subsample_files = ();
    print "Making these subsets:\n";
    foreach my $subset_size (@subset_sizes) {
        print "    $subset_size\n";
        my $subset_file = "$outdir/$subset_size"."$extension";
        if ($subset_size eq 'all') {
            # Make a symlink $infile->$outdir/all.fastq
            `ln -s $readfile $subset_file`;
            print "  Making 'all' subset (symlink from reads file)\n";
        }
        elsif ((!-e $subset_file) || ($overwrite)) {
            #print "  Making '$subset_size' subset\n";
            # Important: if $overwrite is set, we should delete the existing file, if present.
            if ($overwrite) { `rm $subset_file`; }
            
            # Now, let's extract that number of reads and write them to the right file.
            open INFILE, "<", $readfile or die "ERROR: Could not open fastq file\n  $readfile\n$!\n";
            my @buffer = ();    my @output_reads = ();
            my $read_count = 0;
            READIT: while (my $line = <INFILE>) {
                chomp $line;
                push @buffer, $line;
                if (@buffer >= $lines_per_read) {
                    my $read = join "\n", @buffer;
                    push @output_reads, $read;
                    
                    # Apparently this incredibly simple thing just won't work. Wrong number of reads output. 
                    $read_count ++;
                    
                    if (($read_count >= $subset_size) || ($read_count >= $number_of_reads)) {
                        print "        Subset $subset_size: writing ".@output_reads." reads (readcount $read_count) of $number_of_reads in dataset\n";
                        open (OUT, ">", $subset_file) or die "ERROR: Cannot open reads subsample file\n  $subset_file\n$!\n"; 
                        foreach my $read (@output_reads) { print OUT "$read\n"; }
                        close OUT;
                        my $number_of_output_reads = `wc -l $subset_file`;
                        my @sp = split /\s/, $number_of_output_reads;
                        $number_of_output_reads = $sp[0];
                        $number_of_output_reads /= $lines_per_read;
                        print "          ($number_of_output_reads reads actually written to the file)\n";
                        @output_reads = ();
                        $read_count = 0;
                        last READIT;
                    }
                    @buffer = ();
                }
            }
            close INFILE;
        }
        else {
            print "  Subset '$subset_size' already exists; using that\n";
        }
        push @subsample_files, $subset_file;
    }
    
    print "final subsample files:\n";
    foreach my $file (@subsample_files) {
        print "  $file\n";
    }
    return \@subsample_files;
}

sub rechunk {
    # This is meant to split read files back up into a set number of subdivisions, after the whole subsetting business has been handled. It will, hopefully, speed up our chronically slow aligners a bit.
    # I think I overcomplicated this a bit before, by randomly assigning reads to files. Don't do that; just keep plucking 4 (or 2, if fasta) lines out at a time, and switch the output over whenever necessary.
    my $self = shift;
    my $readfile = shift;
    my $num_chunks = $self->{param}{num_chunks};
    my $output_prefix = $self->{param}{output_prefix};
    my $overwrite = $self->{param}{overwrite};
    
    my @chunk_files = ();
    # We may not wish to do any rechunking (in which case we should just return $readfile)
    if (!$num_chunks) {
        push @chunk_files, $readfile;
    }
    else {
        # Take care of some stuff to make sure we're getting an accurate, nonsense-free look at the number of reads in the input file.
        # Adjusting this lets this sub handle either fastq or fasta input. 
        $self->remove_trailing_newlines($readfile);
        my $lines_per_read = 4;
        my ($basename, $parentdir, $extension) = fileparse($readfile, qr/\.[^.]*$/);
        if ($readfile =~ /\.fasta$/) { $lines_per_read = 2; }
        
        # Get number of reads in the file
        # Choose appropriate division
        my $number_of_reads = `wc -l $readfile`;
        my @sp = split /\s/, $number_of_reads;
        $number_of_reads = $sp[0];
        $number_of_reads /= $lines_per_read;
        
        # Calculate a number; this will help us decide which chunk file each read should go into.
        # Round up by 1, in order to not get a file at the end with only one or two reads in it.
        my $reads_per_chunk = int ($number_of_reads / $num_chunks) + 1;
        
        # Make the right output directory for the subset file
        my $dir = dirname($readfile);
        my $subset = basename($readfile);
        $subset =~ s/\.fastq//;
        
        my $chunkdir = $self->directory_check("$dir/$subset");
        
        # First, do a check; do these rechunked files already exist? We want to avoid repeating this is we can help it.
        # Criteria for skipping:
        # Must be correct number of rechunked files
        # All must have the correct number of reads (with a bit of variance)
        my $pattern = "*$extension*";
        my $existing_chunk_files = $self->find_files($chunkdir,$pattern);
        
        print "Rechunking check: can we skip this step?\n";
        # $overwrite must not be set...
        unless ($overwrite) {
            #print "  Overwrite not requested\n";
            # Files must exist...
            unless ($existing_chunk_files->[0] =~ /No such file or directory/) {
                #print "  Files exist in output dir\n";
                # Number of files must match...
                if (@$existing_chunk_files == $num_chunks) {
                    
                    # Number of reads must match...
                    # Check them all
                    my $num_reads_is_fine = 1;
                    foreach my $file (@$existing_chunk_files) {
                        my $number_of_subset_reads =  `wc -l $file`;
                        $number_of_subset_reads /= $lines_per_read;
                        # Bear in mind that the last read file in the set will most likely have fewer reads.
                        # Give a pass if we're on the last file and the number of subset reads equals total number of reads modulus number of chunks.
                        unless (($number_of_subset_reads == $reads_per_chunk) || ((basename($file) eq "$num_chunks.fasta") && ($number_of_subset_reads == ($number_of_reads % $num_chunks)))) {
                            $num_reads_is_fine = 0;
                            #print "    File $file has wrong number of reads! ($number_of_subset_reads vs. $reads_per_chunk)\n";
                        }
                    }
                    
                    # If that's all OK, feel free to skip.
                    if ($num_reads_is_fine == 1) {
                        print "  Yes; existing output with correct number of reads found.\n";
                        my @chunk_files = ();
                        foreach my $i (1..$num_chunks) {
                            push @chunk_files, "$chunkdir/$i.fastq";
                        }
                        return \@chunk_files;
                    }
                    else { print "  No; one or more files has the wrong number of reads\n    (indicating it was probably made under a different subsampling regime).\n"; }
                }
                else { print "  No; incorrect number of files in $chunkdir\n    (".@$existing_chunk_files." vs. $num_chunks)\n"; }
            }
            else { print "  No; cannot find any existing output files in directory\n    $chunkdir\n"; }
        }
        else { print "  No; overwrite requested\n"; }
        
        # There is a minor complication here. Since this sub adds data to files mostly through appending, we should avoid writing to files that already exist.
        # That can be achieved simply by deleting the chunks directory and recreating it every time.
        if ($subset) {
            if (-d "$dir/$subset") { `rm -r $dir/$subset`; }
        }
        $chunkdir = $self->directory_check("$dir/$subset");
        
        # Now, I need to pick out the reads from the input file, and direct them to an output file.
        # These are potentially very big files, so they should be read line by line. That would need a little read buffer type thing, since we want
        # chunks of 4 lines at a time.
        # Make sure to store the names of the chunk files created, so that they can be passed back. 
        my %chunk_files = ();
        open INFILE, "<", $readfile or die "ERROR: Could not open fastq file\n  $readfile: $!\n";
        my @buffer = ();    my @output_reads = ();
        my $read_count = 0; my $outfile = 1;    my $total_read_count = 0;
        while (my $line = <INFILE>) {
            chomp $line;
            push @buffer, $line;
            if (@buffer >= $lines_per_read) {
                push @output_reads, @buffer;
                @buffer = ();
                $read_count ++; $total_read_count ++;
                
                if (($read_count >= $reads_per_chunk) || ($total_read_count >= $number_of_reads)) {
                    $chunk_files{"$chunkdir/$outfile"."$extension"} = 1;
                    open (OUT, ">", "$chunkdir/$outfile"."$extension") or die "ERROR: Cannot open reads subsample file\n  $chunkdir/$outfile"."$extension\n"; 
                    foreach my $l (@output_reads) { print OUT "$l\n"; }
                    close OUT;
                    @output_reads = ();
                    $outfile ++;
                    $read_count = 0;
                }
            }
        }
        close INFILE;
        
        @chunk_files = keys %chunk_files;
        $self->remove_trailing_newlines(\@chunk_files);
    }
    foreach my $file (@chunk_files) { $self->does_file_exist($file); }
    return \@chunk_files;
}

sub remove_trailing_newlines {
    # A useful thing to be able to do when relying on the number of lines being useful information.
    # Make it so this can take an array of files too.
    my $self = shift;
    my $file = shift;
    
    my @files = ();
    if (ref $file eq 'ARRAY') { @files = @$file; }
    else                      { push @files, $file; }
    foreach my $file (@files) { $self->does_file_exist($file); }
    foreach my $f (@files) {
        my $gibberish = '/./,$!d;/^\n*$/{$d;N;};/\n$/ba';
        `sed -e :a -e '$gibberish' $f`;
    }
}

sub extract_list_of_jobs {
    # When given a load of STDOUT data, this will return a list of job IDs that were started within that session.
    my $self = shift;
    my $ramble = shift;
    
    my @lines = split /\n/, $ramble;
    my @jobs = ();
    foreach my $line (@lines) {
        if ($line =~ /Job \<[0-9]+\> is submitted to queue /) {
            # LSF job submitted string
            my @sp = split /[\<\>]/, $line;
            push @jobs, $sp[1];
        }
        elsif ($line =~ /Submitted batch job [0-9]+/) {
            # SLURM job submitted string
            $line =~ s/[^0-9]+//g;
            push @jobs, $line;
        }
    }
    
    if (@jobs == 0) {
        die "PROBABLE ERROR: No jobs appear to have been submitted in the above step.\n";
    }
    else {
        print "Found ".@jobs." job(s) launched by the above command.\n";
    }
    
    return \@jobs;
}

sub done_when_its_done {
    # IMPORTANT: The various processes of this pipeline are set off as separate LSF jobs. If we
    # don't wait for a job to finish before starting the next, the dependent jobs will all
    # just crash. I therefore want to make each process wait until the jobs it has kicked off finish
    # before carrying on to error-checking and the next process.
    # When given a list of job IDs, this sub will wait until they all finish (periodically
    # checking in on them all) before allowing the script to move on.
    my $self = shift;
    my $check_jobs = shift;
    my $scheduler = $self->{config}{jobsys};
    
    # A list of jobs should always be supplied.
    if (!$check_jobs) { die "ERROR: Failed to find a list of job IDs where expected.\n"; }
    
    # Check the list of current jobs every minute or so, until all the input jobs have disappeared
    # from the list.
    my @intersection = ();
    my $checks = ();
    my $waited = 0;
    my $failed = 0;
    my $failure_attempts = 20;
    do {
        $checks ++;
        # Get running job IDs using the relevant command for the job scheduler, and awk
        my $running_jobs = ();
        if ($scheduler eq 'lsf')      { $running_jobs = `bjobs | awk '{print \$1}'`; }
        elsif ($scheduler eq 'slurm') { $running_jobs = `squeue -u \$USER | awk '{print \$1}'`; }
        
        # This should never, ever be empty, because if it is, it means either the job scheduler is not working as it ought to be, or much more likely, something has gone wrong with this code. The pipeline is always run as a job, so if nothing else, it should find itself.
        #if (!$running_jobs) { die "ERROR: Job scheduler cannot find a list of running jobs.\n"; }
        # WRONG. The job scheduler will intermittently fail to return a list of running jobs, for various and sundry reasons that ultimately make no difference to this pipeline (because the jobs continue to run regardless). It always sorts itself out when left alone for a bit, but if the pipeline happens to check for jobs in that time, it will find one - and will wrongly conclude that all jobs have finished and the pipeline can proceed. We need to make it... not do that. That's what $failed is for.
        # To counteract the possibility of getting stuck in an infinite loop, though, the number of failed attempts to access the job list should be capped - albeit pretty generously. 
        if ($running_jobs) {
            $failed = 0;
            my @running_jobs = split /\n/, $running_jobs;
            my $line = shift @running_jobs;
            
            # Compare using List::Compare
            # BUT only if the scheduler has indeed returned a list of jobs!
            # (Don't know what slurm does when it fails to return a list yet. Update once I do...)
            unless ($line =~ /Please wait/) {
                my $lc = List::Compare->new($check_jobs, \@running_jobs);
                @intersection = $lc->get_intersection;
            }
        }
        else {
            $failed ++;
            if ($failed >= $failure_attempts) { die "ERROR: Unable to access list of running jobs after $failure_attempts attempts.\n"; }
        }
        
        # If any jobs in the checklist are still active, then wait a short time before checking again.
        # No point checking ten billion times per minute.
        if ((@intersection >= 1) || ($failed > 0)) {
            if ($checks <= 10)    { sleep 5; $waited += 5; }
            elsif ($checks <= 20) { sleep 10; $waited += 10; }
            else                  { sleep 15; $waited += 15; }
        }
    }
    until ((@intersection == 0) && ($failed == 0));
    
    # Wait another 60 seconds; perl/the cluster apparently needs this in order to sort out newly created files.
    sleep 60;
    $waited += 60;
    
    my $time_now = localtime();
    print "All pending jobs for this step completed by $time_now\n(Waited for $waited seconds)\n";
    return "Done";
}

sub single_job_check {
    # This is used to check that submitted jobs have completed successfully (and if they haven't, cause the pipeline
    # to stop immediately rather than try to continue).
    # Pass in a command and a job ID, and this sub will figure it out from there.
    # I've included an option on whether to die or not because we might want to test a set of jobs at once, and I want to
    # know what proportion of them have failed.
    my $self = shift;
    my $command = shift;
    my $job_input = shift;
    my $die = shift;
    my $dont_wait = shift;
    
    # $jobs may be an array reference or a string. I'm going to assume that if it's an array reference, it only contains
    # one job.
    my $jobID = ();
    if (ref $job_input eq 'ARRAY') {
        $jobID = $job_input->[0]
    }
    
    # First, we need the log file path from the command.
    # (This may not be supplied; if so, there's not a lot we can do).
    my @bits = split / /, $command;
    my $logpath = ();
    my $c = 0;
    foreach my $bit (@bits) {
        if ($bit eq '-oo') {
            $logpath = [$c+1]
        }
        $c ++;
    }
    # If no logpath found, just continue silently.
    my $error = 0;
    if ($logpath) {
        if (-e $logpath) {
            # Open the file, check the status
            open(LOG, "<", $logpath) or die "ERROR: Could not open expected log file for the following job:\n  $logpath\n$!\n";
            my @log = <LOG>; close LOG;
            # We're looking for one particular line.
            # When we reach a line with "Resource usage summary in, we know we've gone past it."
            LINES: foreach my $line (@log) {
                if ($line =~ /Resource usage summary/) { last LINES; }
                if ($line =~ /Exited with exit code/)  {
                    print "ERROR: Log reports failure for the following job:\n  $logpath\n";
                    if ($die) { die ""; }
                    else { $error ++; }
                }
            }
        }
        else {
            # If the logpath is missing, this job has probably failed.
            # Wait a few seconds and try again, and if it still fails, kill the pipeline.
            # We can achieve that with a bit of recursion.
            if ($dont_wait) {
                print "ERROR: Could not find expected log file for the following job:\n  $logpath\n";
                if ($die) { die ""; }
                else { $error ++; }
            }
            else {
                sleep 10;
                $self->jobs_check($command, $job_input, $die, 1);
            }
        }
    }
    return $error;
}

sub jobs_check {
    # The above sub reads the log for a single job, and kills the pipeline if the log reports failure.
    # But parts of this pipeline have jobs run in parallel - most notably, in the alignment stage.
    # I currently check that those jobs have all finished before continuing, but I want to check that they were all
    # successful too.
    # The currently active jobs (that have built up since the start of the pipeline, or since the last time this sub was called)
    # are in $jobs_check_buffer - a hash reference, with the job IDs as the key and the commands as the value.
    my $self = shift;
    my $jobs = shift;
    
    my $error = 0;
    foreach my $job (keys %jobs_check_buffer) {
        my $e = single_job_check($jobs_check_buffer{job},$job);
        if ($e) { $error ++; }
        
    }
    if ($error > 0) { die "Found $error failed jobs out of ".@$jobs." submitted\n"; }
    
    # Clear the jobs buffer out to prevent needless double-checking of jobs that we already know have been successful.
    %jobs_check_buffer = ();
}

sub run_nextclip {
    my $self = shift;
    my $readsfiles = shift;
    my $output_path = shift;
    my $output_path_slash =$output_path;
    $output_path_slash =~ s!/*$!/!; # Add a trailing slash if none is present
    my $log_path = $self->{param}{log_path};
    my $queue = $self->{config}{queue};
    my $run_nextclip = $self->{config}{run_nextclip};
    my $nextclip_version = $self->{config}{nextclip_version};
    my $remove_pcr_duplicates = $self->{config}{remove_pcr_duplicates};
    my $overwrite = $self->{param}{overwrite};
    my $scheduler = $self->{config}{jobsys};
    
    if (!$readsfiles) { die "ERROR: Reads files not supplied.\n"; }
    if (@$readsfiles == 0) { die "ERROR: Empty reads file array.\n"; }
    foreach my $file (@$readsfiles) { $self->does_file_exist($file); }
    
    my $readsfiles_out = ();
    # Are we actually going to run NextClip? If not, just return the input files and don't do anything.
    if ($run_nextclip eq 'no') {
        $readsfiles_out = $readsfiles;
    }
    # We may have a single-end run, which nextclip can't handle. That means we have to do something else.
    elsif ($self->{config}{reads_type} eq 'single end') {
        print "\tThis is a single-end run, so NextClip cannot be run.\n";
        # If the user has requested removal of PCR duplicates, then we need to run another sub that can do it.
        # Otherwise, we just need to return a warning, pass the filepath back, and continue merrily on our way.
        if ($remove_pcr_duplicates eq 'yes') {
            print "PCR duplicates will be removed by Metapipe's remove_pcr_duplicates function instead.\n";
            my $outfile = $self->remove_pcr_duplicates($readsfiles->[0]);
            push @$readsfiles_out, $outfile;
        }
        else {
            print "\tNo further action taken.\n";
            $readsfiles_out = $readsfiles;
        }
    }
    else {
        foreach my $in (@$readsfiles) {
            push @$readsfiles_out, "$output_path/".basename($in);
        }
        
        # Got to pick out the read 1 and 2 files from the inputs. (I put them in in this order, for definite)
        my $read1file = $readsfiles->[0];
        my $read2file = $readsfiles->[1];
        my $read1file_out = $readsfiles_out->[0];
        my $read2file_out = $readsfiles_out->[1];
        
        # How can I check if this step has been done already?
        # Look for output files in the output path, of course.
        # Check if NextClip files exist. If not, run NextClip; if they do, skip it.
        my $files = $self->find_files($output_path, "results_[ABCD]_*.fastq");
        if ((@$files > 0) && (!$overwrite)) {
            print "\tFound existing NextClip output files:\n";
            foreach my $ncfile (@$files) { print "  $ncfile\n"; }
            print "\tSkipping NextClip!\n";
        }
        else {
            # Set up the job
            my $jobname = basename($self->{param}{output_prefix});
            
            my $cmd = "nextclip -i $read1file -j $read2file -o $output_path_slash"."results ";
            if ($remove_pcr_duplicates eq 'yes') { $cmd .= "--remove_duplicates "; }
            
            my $opts->{source} = ["nextclip-$nextclip_version"];
            $opts->{jobname} = "NextClip_$jobname";
            $opts->{log} = "$log_path/nextclip_run_log.$scheduler";
            $opts->{threads} = 1;
            
            my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
            print "NextClip command:\n  $outcmd\n";
            $self->done_when_its_done($jobs);
        }
        
        # Then, after that, check if catted reads files exist. If not, make them; if they do, just set the names as the outputs.
        $files = $self->find_files($output_path, basename($read1file));
        push @$files, @{$self->find_files($output_path, basename($read1file))};
        if ((@$files > 0) && (!$overwrite)) {
            print "Found NextClip output files (merged reads, all types - A, B, C, D):\n";
            foreach my $ncfile (@$files) { print "  $ncfile\n"; }
            print "Skipping concatenation step!\n";
        }
        else {
            # Note on removal of PCR duplicates:
            # Nextclip dumps all its output into the place specified in the output prefix (obviously). The actual read files are named results_[A,B,C,D]_R[1,2], though. The majority of the reads will usually end up in category D, because we usually won't be doing metagenomics on Nextera long mate pair libraries.
            # We want all of those reads, so we have to use cat and pipe them off to a single file, or something. Nonetheless, we should then have the full set of reads, minus PCR duplicates. (If we want a specific set of reads pulled out by NextClip in the future, this is the best place to get it).
            print "Nextclip reads output files:\n  $read1file_out\n  $read2file_out\n";
            
            my $cmd = "cat $output_path/results_*_R1.fastq > $read1file_out";
            my $opts->{source} = ["nextclip-$nextclip_version"];
            $opts->{jobname} = "Cat_R1";
            $opts->{threads} = 1;
            my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
            print "Concatenation of NextClip results command R1:\n  $outcmd\n";
            $self->done_when_its_done($jobs);
            
            $cmd = "cat $output_path/results_*_R2.fastq > $read2file_out";
            $opts->{source} = ["nextclip-$nextclip_version"];
            $opts->{jobname} = "Cat_R2";
            $opts->{threads} = 1;
            ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
            print "Concatenation of NextClip results command R2:\n  $outcmd\n";
            $self->done_when_its_done($jobs);
        }
    }
    if (!$readsfiles_out) { die "ERROR: Output reads files not supplied.\n"; }
    if (@$readsfiles_out == 0) { die "ERROR: Empty output reads file array.\n"; }
    foreach my $outfile (@$readsfiles_out) { $self->does_file_exist($outfile); }
    $self->halt('nextclip');
    return $readsfiles_out;
}

sub remove_pcr_duplicates {
    my $self = shift;
    my $file = shift;
    my $output_prefix = $self->{param}{output_prefix};
    my $log_path = $self->{param}{log_path};
    my $overwrite = $self->{param}{overwrite};
    my $check_n_characters = $self->{config}{check_first_n_characters};
    my $scriptsdir = $self->{config}{scriptsdir};
    my $scheduler = $self->{config}{jobsys};
    
    if (!$file) { die "ERROR: Reads files not supplied.\n"; }
    $self->does_file_exist($file);
    
    my $lines_per_read = 4;
    my ($basename, $parentdir, $extension) = fileparse($file, qr/\.[^.]*$/);
    if ($extension =~ /\.fasta$/) { $lines_per_read = 2; }
    
    # I want to keep the output of this read separate, so I need to make an output directory.
    my $output_dir = $self->directory_check("$output_prefix/reads/nextclip/pcr_duplicate_removal");
    my $outfile = "$output_dir/$basename"."$extension";
    
    # To make a more effective use of cluster memory when running this crude and fairly memory-intensive function, I've moved it off to a separate script.
    # That script is in $scriptsdir.
    print "Removing PCR duplicates from file\n  $file\n";
    if ((!-e $outfile) || ($overwrite)) {
        my $cmd = "perl $scriptsdir/remove_pcr_duplicates.pl $file ";
        if ($check_n_characters) { $cmd .= "$check_n_characters "; }
        $cmd .= "> $outfile";
        my $opts->{jobname} = "Remove_PCR_duplicates";
        $opts->{threads} = 1;
        $opts->{log} = "$log_path/deduplicate_$basename.$scheduler";
        
        my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
        print "PCR duplicate removal command:\n  $outcmd\n";
        $self->done_when_its_done($jobs);
        $self->remove_trailing_newlines($outfile);
    }
    else {
        print "Output file exists; skipping this step!\n";
    }
    
    #print "PCR duplicates removed in reads file\n  $outfile\n";
    $self->does_file_exist($outfile);
    
    return $outfile;
}

sub run_fastqc {
    my $self = shift;
    my $readsfiles = shift;
    my $output_path = shift;
    my $log_path = $self->{param}{log_path};
    my $queue = $self->{config}{queue};
    my $fastqc_version = $self->{config}{fastqc_version};
    my $overwrite = $self->{param}{overwrite};
    my $scheduler = $self->{config}{jobsys};
    
    if (!$readsfiles) { die "ERROR: Reads files not supplied.\n"; }
    if (@$readsfiles == 0) { die "ERROR: Empty reads file array.\n"; }
    foreach my $file (@$readsfiles) { $self->does_file_exist($file); }
    
    # Since FastQC gets run twice, we need to set two different log file names to prevent the later ones overwriting the earlier ones. 
    my $fastqc_log_path = "$log_path/fastqc_";
    if ($self->{config}{trimming_done}) { $fastqc_log_path .= "trimmed.$scheduler"; }
    else                                { $fastqc_log_path .= "untrimmed.$scheduler"; }
    
    my $outfiles = ();
    READS: foreach my $readfile (@$readsfiles) {
        my $fastqcfile = basename($readfile);
        $fastqcfile =~ s/\.fastq/_fastqc\.zip/g;
        
        # How can I check if this step has been done already?
        # Look for output files in the output path, of course.
        unless ($overwrite) {
            my $zippedfiles = $self->find_files($output_path, '*fastqc.zip');
            my $unzippedfiles = $self->find_files($output_path, '*fastqc.html');
            my @files = @$zippedfiles;
            push @files, @$unzippedfiles;
            if (@files > 0) {
                print "Found existing FastQC output files:\n";
                foreach my $fqfile (@files) { print "  $fqfile\n"; }
                print "Skipping FastQC!\n";
                push @$outfiles, "$output_path/$fastqcfile";
                next READS;
            }
            else { print "No existing output file found; proceeding with FastQC.\n"; }
        }
        
        my $jobname = basename($self->{param}{output_prefix});
        # Problem reported here; sometimes FastQC appears to hang, asking if we want to "replace $some_directory/Icons/fastqc_icon.png". Options are yes, no, all and none.
        # Let's go with 'all', which we activate by typing 'A'; hence the echo statement in the command below.
        # Note: may or may not remain necessary now that I've tried to fix the way logs used to overwrite each other. 
        my $cmd = "fastqc $readfile --outdir $output_path; ".'echo "A"';
        my $opts->{source} = ["fastqc-$fastqc_version"];
        $opts->{jobname} = "FastQC_$jobname";
        $opts->{log} = "$fastqc_log_path";
        $opts->{threads} = 1;
        
        my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
        print "FastQC command:\n  $outcmd\n";
        $self->done_when_its_done($jobs);
        
        push @$outfiles, "$output_path/$fastqcfile"
    }
    foreach my $file (@$outfiles) { $self->does_file_exist($file); }
    $self->halt('fastqc');
    return ($outfiles);
}

sub examine_fastqc_results {
    my $self = shift;
    my $fastqc_results_files = shift;
    my $trimming_done = shift;
    
    if (!$fastqc_results_files) { die "ERROR: FastQC results files not supplied.\n"; }
    if (@$fastqc_results_files == 0) { die "ERROR: Empty FastQC results file array.\n"; }
    foreach my $file (@$fastqc_results_files) { $self->does_file_exist($file); }
    
    foreach my $fastqc_results_file (@$fastqc_results_files) {
        # Unzip $fastqc_results_file
        # Read fastqc_data.txt in directory created from that
        # Locate the adapter content section of that file
        # Prepare warnings, if necessary, based on what's in there
        my $fastqc_results_dir = $self->unzip_fastqc_results($fastqc_results_file);
        
        # Get the name of that directory (it'll be used for printing output to the user)
        my @dirsplit = split /\//, $fastqc_results_dir;
        my $results_dirname = $dirsplit[-1];
        
        open(RESULTS, "<", "$fastqc_results_dir/fastqc_data.txt") or die "ERROR: Cannot open unzipped FastQC data\n  $fastqc_results_dir/fastqc_data.txt\n";
        my @data = <RESULTS>;
        close RESULTS;
        
        my $line = ();
        do {
            $line = shift @data;
        }
        until ($line =~ /Adapter Content/);
        chomp $line;
        #my @linesplit = split /\t/, $line;
        #my $test_status = $linesplit[1];
        #chomp $test_status;
        #if ($test_status =~ /pass/) {
        if ($line =~ /pass/) {
            print "FastQC test: adapter content for dataset\n  $results_dirname\nwithin acceptable limits!!\n";
        }
        else {
            if ($trimming_done) {
                # If trimming has been done and adapter content is still too high, throw a big wobbly.
                print "ERROR: Adapter content too high after trimming in dataset\n  $results_dirname!\n";
                print "$line\n";
                my $ac_section = ();
                do {
                    $ac_section = shift @data;
                    print $ac_section;
                }
                until ($ac_section =~ /END_MODULE/);
            }
            else {
                # If trimming hasn't been done, throw a minor warning.
                print "WARN: Adapter content too high (before trimming) in dataset\n  $results_dirname\n";
            }
        }
    }
}

sub unzip_fastqc_results {
    my $self = shift;
    my $fastqc_file = shift;
    my $overwrite = $self->{param}{overwrite};
    
    $self->does_file_exist($fastqc_file);
    # Get the directory alone
    my $fastqc_dir = dirname($fastqc_file);
    # Get the basename of the file; we'll use it to check if this has already been unzipped
    my $fastqc_file_basename = basename($fastqc_file);
    # If it's been unzipped, there will be a directory that matches the name of the file, except for its extension.
    my $fastqc_results_dir = $fastqc_file_basename;
    $fastqc_results_dir =~ s/\.zip//g;
    
    unless ($overwrite) {
        if (-d "$fastqc_dir/$fastqc_results_dir") {
            return "$fastqc_dir/$fastqc_results_dir";
        }
    }
    
    # We have to change directory to $fastqc_dir, then do the unzip.
    # (Then change it back).
    my $pwd = cwd();
    chdir($fastqc_dir);
    my $unzip = `unzip $fastqc_file`;
    chdir($pwd);
    
    return "$fastqc_dir/$fastqc_results_dir";
}

sub run_trimming {
    my $self = shift;
    my $readsfiles = shift;
    my $run_trimming = $self->{config}{run_trimming};
    my $log_path = $self->{param}{log_path};
    my $queue = $self->{config}{queue};
    my $adaptersfile = $self->{config}{adaptersfile};
    my $seed_mismatches = $self->{config}{seed_mismatches};
    my $palindrome_clip_threshold = $self->{config}{palindrome_clip_threshold};
    my $simple_clip_threshold = $self->{config}{simple_clip_threshold};
    my $trimming_sliding_window = $self->{config}{trimming_sliding_window};
    my $trimming_min_length = $self->{config}{trimming_min_length};
    my $threads = $self->{config}{threads};
    my $trimming_output_path = shift;
    my $trimmomatic_version = $self->{config}{trimmomatic_version};
    my $overwrite = $self->{param}{overwrite};
    my $memory = $self->{config}{trimming_memory};
    my $readtype = $self->{config}{reads_type};
    my $scheduler = $self->{config}{jobsys};
    
    if (!$readsfiles) { die "ERROR: Reads files not supplied.\n"; }
    if (@$readsfiles == 0) { die "ERROR: Empty reads file array.\n"; }
    foreach my $file (@$readsfiles) { $self->does_file_exist($file); }
    
    my $outfiles = ();
    if ($run_trimming eq 'no') {
        print "\tConfig file requests no trimming; skipping this step\n";
        push @$outfiles, @$readsfiles;
    }
    else {
        $self->does_file_exist($adaptersfile);
        
        foreach my $readfile (@$readsfiles) {
            $self->does_file_exist($readfile);
        }
        
        # Trimmomatic actually uses somewhat different commands for paired-end and single-end runs, so I'll just have to arrange it into
        # different code blocks.
        
        my ($readsfiles_trimmed, $readsfiles_trimmed_single) = ();
        foreach my $readfile (@$readsfiles) {
            my $readfile_trimmed = "$trimming_output_path/".basename($readfile);
            push @$readsfiles_trimmed, $readfile_trimmed;
            my $readfile_trimmed_single = $readfile_trimmed;
            $readfile_trimmed_single =~ s/\.fastq/_single.fastq/;
            push @$readsfiles_trimmed_single, $readfile_trimmed;
            push @$outfiles, $readfile_trimmed;
        }
        
        unless ($overwrite) {
            my $exists = 1;
            foreach my $file (@$readsfiles_trimmed) {
                unless (-e $file) { $exists = 0; }
            }
            if ($exists == 1) {
                print "Found all expected trimming output files.\nSkipping trimming!\n";
               $self->halt('trimming');
                return $readsfiles_trimmed;
            }
            
            else { print "--No existing file found; proceeding with trimming.\n"; }
        }
        
        my $hot_air = "-R \"rusage[mem=$memory] span[hosts=1]\"";
        
        my $jobname = basename($self->{param}{output_prefix});
        
        # Due to the slightly different commands used for single-end and paired-end runs, this command will need a bit of perl's magic at
        # string handling.
        my $cmd = "java -jar /tgac/software/testing/trimmomatic/0.30/x86_64/bin/trimmomatic-0.30.jar ";
        if ($readtype eq 'paired end') { $cmd .= "PE"; }
        else                           { $cmd .= "SE"; }
        $cmd .= " -phred33 -threads $threads";
        foreach my $readfile (@$readsfiles) { $cmd .= " $readfile"; }
        foreach my $i (1..@$readsfiles) {
            $i --; $cmd .= " ".$readsfiles_trimmed->[$i];
            if ($readtype eq 'paired end') { $cmd .= " ".$readsfiles_trimmed_single->[$i]; }
        }
        # I don't know what the :2:30:10 bit after $adaptersfile does, but it seems to be important.
        # (It's a bunch of trimming parameters).
        $cmd .= " ILLUMINACLIP:$adaptersfile:$seed_mismatches:$palindrome_clip_threshold:$simple_clip_threshold SLIDINGWINDOW:$trimming_sliding_window MINLEN:$trimming_min_length ";
        
        my $opts->{source} = ["source trimmomatic-$trimmomatic_version"];
        $opts->{jobname} = "Trimming_$jobname";
        $opts->{log} = "$log_path/trim.$scheduler";
        $opts->{threads} = $threads;
        $opts->{memory} = $memory;
        
        my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
        print "TRIMMOMATIC command:\n  $outcmd\n";
        $self->done_when_its_done($jobs);
    }
    foreach my $file (@$outfiles) { $self->does_file_exist($file); }
    $self->halt('trimming');
    return $outfiles;
}

sub run_kontaminant {
    my $self = shift;
    my $readsfiles = shift;
    my $filtering_output_dir = shift;
    my $reference = shift;
    my $run_filtering = $self->{config}{run_filtering};
    my $log_path = $self->{param}{log_path};
    my $database =$self->{config}{kontaminant_database};
    my $queue = $self->{config}{queue};
    my $kontaminant_version = $self->{config}{kontaminant_version};
    my $overwrite = $self->{param}{overwrite};
    my $mem_width = $self->{config}{kontaminant_mem_width};
    my $mem_height = $self->{config}{kontaminant_mem_height};
    my $memory = $self->{config}{filtering_memory};
    my $scheduler = $self->{config}{jobsys};
    
    if (!$readsfiles) { die "ERROR: Reads files not supplied.\n"; }
    if (@$readsfiles == 0) { die "ERROR: Empty reads file array.\n"; }
    foreach my $readfile (@$readsfiles) { $self->does_file_exist($readfile); }
    
    my $outfiles = ();
    my ($r1outfile, $r2outfile) = ();
    if ($run_filtering =~ /yes/) {
        $self->directory_check("$log_path/kontaminant");
        
        foreach my $readfile (@$readsfiles) {
            push @$outfiles, "$filtering_output_dir/filtered_".basename($readfile);
        }
        
        unless ($overwrite) {
            if ((-e $r1outfile) && (-e $r2outfile)) {
                print "Found \n$r1outfile \nand \n$r2outfile\nSkipping filtering!\n";
                $self->halt('kontaminant');
                return ($r1outfile, $r2outfile);
            }
            #else { print "--No existing file found; proceeding with filtering.\n"; }
        }
        
        my $jobname = basename($self->{param}{output_prefix});
        my $cmd = "kontaminant -f";
        my $c = 0;
        foreach my $readfile (@$readsfiles) {
            $c++;
            $cmd .= " -$c $readfile";
        }
        $cmd .= " -c $reference -d $database -k 21 -o $filtering_output_dir/filtered_ -r $filtering_output_dir/removed_ -p $log_path/kontaminant -n $mem_height -b $mem_width\" ";
        my $opts->{source} = ["kontaminant-$kontaminant_version"];
        $opts->{jobname} = "Filter_$jobname";
        $opts->{log} = "$log_path/kontaminant.$scheduler";
        $opts->{memory} = $memory;
        $opts->{threads} = 1;
        
        my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
        print "Kontaminant command:\n  $outcmd\n";
        $self->done_when_its_done($jobs);
    }
    else {
        print "\tConfig file requests no filtering; skipping this step\n";
        $outfiles = $readsfiles;
    }
    foreach my $file (@$outfiles) { $self->does_file_exist($file); }
    $self->halt('kontaminant');
    return $outfiles;
}

sub run_flash {
    my $self = shift;
    my $readsfiles = shift;
    my $output_dir = shift;
    my $output_file = shift;
    my $log_path = $self->{param}{log_path};
    my $queue = $self->{config}{queue};
    my $flash_version = $self->{config}{flash_version};
    my $overwrite = $self->{param}{overwrite};
    my $readstype = $self->{config}{reads_type};
    my $scheduler = $self->{config}{jobsys};
    
    if (!$readsfiles) { die "ERROR: Reads files not supplied.\n"; }
    if (@$readsfiles == 0) { die "ERROR: Empty reads file array.\n"; }
    foreach my $readfile (@$readsfiles) { $self->does_file_exist($readfile); }
    
    # If running on single-end data, we can't (and don't need to) use flash.
    if ($readstype eq 'single end') {
        print "\tFLASH is unnecessary for single-end runs.\n";
        $self->halt('flash');
        return $readsfiles->[0];
    }
    
    my $read1file = $readsfiles->[0];
    my $read2file = $readsfiles->[1];
    $self->does_file_exist($read1file);
    $self->does_file_exist($read2file);
    
    # The string in $output_file isn't actually the full name of the output file. FLASH creates it by adding this:
    my $full_output_file = "$output_file.extendedFrags.fastq";
    
    unless ($overwrite) {
        if (-e "$output_dir/$full_output_file") {
            print "--Found \n$output_dir/$full_output_file \nSkipping FLASH!\n";
            $self->halt('flash');
            return "$output_dir/$full_output_file";
        }
        else { print "--No existing file found; proceeding with FLASH.\n"; }
    }
    
    # FLASH apparently doesn't like absolute file paths; we have to change directory, then set a relative file path.
    # (Then change it back).
    my $pwd = cwd();
    chdir($output_dir);
    print "Changed working directory:\n$pwd to\n$output_dir\n";
    
    my $jobname = basename($self->{param}{output_prefix});
    my $cmd = "flash $read1file $read2file -M 150 -o $output_file";
    my $opts->{source} = ["flash-$flash_version"];
    $opts->{jobname} = "Flash_$jobname";
    # NOTE: Log for FLASH will be set in the output directory, because there's no guarantee that the $log_path is still directly accessible from where we've changed dir to.
    $opts->{log} = "flash.$scheduler";
    $opts->{threads} = 1;
    
    my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
    print "FLASH command:\n  $outcmd\n";
    $self->done_when_its_done($jobs);
    
    # Check for errors
    # Output files should have the value supplied as $output_file, but with a couple of different extensions indicating what they are.
    
    # Remember to change back to original working directory.
    chdir($pwd);
    print "Changed working directory back to\n  $pwd\n";
    
    $self->halt('flash');
    # Return full file path, so that other stuff can use it
    return "$output_dir/$full_output_file";
}

sub convert_fastq_to_fasta {
    # Converts supplied fastq file (assume path included) into fasta in set output location. 
    my $self = shift;
    my $fastq_filepath = shift;
    my $outdir = shift;
    my $log_path = $self->{param}{log_path};
    my $queue = $self->{config}{queue};
    my $scriptsdir = $self->{config}{scriptsdir};
    my $overwrite = $self->{param}{overwrite};
    
    $self->does_file_exist($fastq_filepath);
    
    my $fastq_filename = basename($fastq_filepath);
    my $fasta_filename = $fastq_filename;
    $fasta_filename =~ s/\.fastq/\.fasta/;
    my $fasta_filepath = "$outdir/$fasta_filename";
    
    unless ($overwrite) {
        if (-e $fasta_filepath) {
            print "Found existing FASTA output files in\n  $fasta_filepath \nSkipping FASTQ->FASTA coversion!\n";
            return $fasta_filepath;
        }
        #else { print "--No existing file found; proceeding with FASTQ->FASTA conversion.\n"; }
    }
    
    # This function is one of a few homebrewed scripts, which I'll need to include a path to as a config file option.
    my $jobname = basename($self->{param}{output_prefix});
    my $cmd = "perl $scriptsdir/fastq_to_fasta.pl $fastq_filepath > $fasta_filepath";
    my $opts->{jobname} = "FASTQ_to_FASTA";
    $opts->{threads} = 1;
    
    # Note that we just kick the jobs off and return job IDs here, since it's very likely we'll want to do this to multiple files at the same time.
    my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);    
    return ($fasta_filepath, $jobs);
}

sub list_aligners {
    # Returns a list of all the aligners we use (as dictated by inputs). Annoyingly, I can't think of a way to fully automate this, so it has to be hard-coded.
    my $self = shift;
    my @aligners_used = ();
    if($self->{param}{blastn})      { push @aligners_used, "blastn"; }
    if($self->{param}{blastx})      { push @aligners_used, "blastx"; }
    if($self->{param}{rapsearch})   { push @aligners_used, "rapsearch"; }
    if($self->{param}{diamond})     { push @aligners_used, "diamond"; }
    return \@aligners_used;
}


sub get_aligner_used {
    # In a couple of situations, we want to know which aligner was used to create a particular data file. That information is contained within the path; it's easily accessible enough, but more convenient to get it using a function. 
    my $self = shift;
    my $path = shift;
    
    # Going to make the assumption that if multiple paths are supplied, they're all from the same aligner.
    if (ref $path eq 'ARRAY') { $path = $path->[0]; }
    
    my $aligners = $self->list_aligners();
    my @path = split /\//, $path;
    my $lc = List::Compare->new($aligners, \@path);
    my @intersection = $lc->get_intersection;
    return $intersection[0];
}

sub run_alignments {
    # Run BLAST on a single sample/subset.
    # I'm going to assume flashed input for now, but handling paired-end data can be added if necessary.
    # This can be thought of as a kind of wrapper. We may want to call multiple different alignments; this collects the results of everything we run.
    my $self = shift;
    #my $db = $self->{config}{blast_database};
    my $query = shift;  # MUST be FASTA format!
    my $blastn = $self->{param}{blastn};
    my $blastx = $self->{param}{blastx};
    my $rapsearch = $self->{param}{rapsearch};
    my $diamond = $self->{param}{diamond};
    #my $output_prefix = $self->{param}{output_prefix};
    
    $self->does_file_exist($query);
    
    # It's possible to make the aligners run in parallel here.
    # Set them off using the appropriate subs. They'll return their output file paths, and also a list of jobs.
    # Add all the jobs to a big list, then call the done_when_its_done function on them.
    # I've also used a hash, which allows us to separate the jobs out by the aligner they represent.
    my (@results, @jobs) = ();
    my %jobs = ();
    if ($blastn)    {
        my ($results, $jobs) = $self->run_blastn($query);
        push @results, $results;
        if ($jobs) {
            push @jobs, @$jobs;
            push @{$jobs{blastn}}, @$jobs;
        }
    }
    if ($blastx)    {
        my ($results, $jobs) = $self->run_blastx($query);
        push @results, $results;
        if ($jobs) {
            push @jobs, @$jobs;
            push @{$jobs{blastx}}, @$jobs;
        }
    }
    if ($diamond)   { 
        my ($results, $jobs) = $self->run_diamond($query);
        push @results, $results;
        if ($jobs) {
            push @jobs, @$jobs;
            push @{$jobs{diamond}}, @$jobs;
        }
    }
    if ($rapsearch) {
        my ($results, $jobs) = $self->run_rapsearch($query);
        push @results, $results;
        if ($jobs) {
             push @jobs, @$jobs;
            push @{$jobs{rapsearch}}, @$jobs;
        }
    }
    
    if (@jobs > 0) { my $done = $self->done_when_its_done(\@jobs); }
    $self->halt('alignment');
    return (\@results, \%jobs);
}

sub run_blastn {
    my $self = shift;
    my $query = shift;
    my $alignment_base_dir = $self->{config}{alignment_base_dir};    #my $log_path = $self->{param}{log_path};
    my $log_path = $self->{param}{log_path};
    my $db = $self->{config}{blastn_database};
    my $sample_id = $self->{param}{sample_id};
    my $queue = $self->{config}{queue};
    my $threads = $self->{config}{threads};
    my $blast_version = $self->{config}{blast_version};
    my $overwrite = $self->{param}{overwrite};
    my $memory = $self->{config}{blast_memory};
    my $evalue = $self->{config}{blast_evalue};
    my $scheduler = $self->{config}{jobsys};
    
    $log_path = $self->directory_check("$log_path/aligners");
    $log_path = $self->directory_check("$log_path/blastn");
    
    print "    Run BLASTn\n";
    $self->does_file_exist($query);
    
    # Get filename alone (no path)
    my $query_filename = basename($query);
    
    # Make an output directory for this aligner
    my $outdir = $self->directory_check("$alignment_base_dir/blastn");
    my $outfile = "$outdir/$query_filename.txt";
    
    unless ($overwrite) {
        if (-e $outfile) {
            print "    Found existing BLASTn output\n      $outfile\n    Skipping BLASTn!\n";
            return $outfile;
        }
        #else { print "--No existing file found; proceeding with BLASTn.\n"; }
    }
    
    my $jobname = basename($self->{param}{output_prefix});
    my $queryname = $query_filename;
    $queryname =~ s/[^0-9]//g;
    $jobname .= "_$queryname";
    
    my $cmd = "blastn -num_threads $threads -db $db -query $query -evalue $evalue -out $outfile";
    my $opts->{source} = ["blast-$blast_version"];
    $opts->{jobname} = "BLASTN_$jobname";
    $opts->{threads} = $threads;
    $opts->{memory} = $memory;
    $opts->{log} = "$log_path/".basename($query).".$scheduler";
    
    my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
    print "      $outcmd\n";
    
    # Check for errors
    
    # Return results
    return ($outfile, $jobs);
}

sub run_blastx {
    my $self = shift;
    my $query = shift;
    my $alignment_base_dir = $self->{config}{alignment_base_dir};
    my $log_path = $self->{param}{log_path};
    my $db = $self->{config}{blastx_database};
    my $sample_id = $self->{param}{sample_id};
    my $queue = $self->{config}{queue};
    my $threads = $self->{config}{threads};
    my $blast_version = $self->{config}{blast_version};
    my $overwrite = $self->{param}{overwrite};
    my $memory = $self->{config}{blast_memory};
    my $evalue = $self->{config}{blast_evalue};
    my $scheduler = $self->{config}{jobsys};
    
    $log_path = $self->directory_check("$log_path/aligners");
    $log_path = $self->directory_check("$log_path/blastx");
    
    print "    Run BLASTx\n";
    $self->does_file_exist($query);
    
    # Get filename alone (no path)
    my $query_filename = basename($query);
    
    # Make an output directory for this aligner
    my $outdir = $self->directory_check("$alignment_base_dir/blastx");
    my $outfile = "$outdir/$query_filename.txt";
    
    unless ($overwrite) {
        if (-e $outfile) {
            print "    Found existing BLASTx output\n      $outfile\n    Skipping BLASTx!\n";
            return $outfile;
        }
        #else { print "--No existing file found; proceeding with BLASTx.\n"; }
    }
    
    my $jobname = basename($self->{param}{output_prefix});
    my $queryname = $query_filename;
    $queryname =~ s/[^0-9]//g;
    $jobname .= "_$queryname";
    
    my $cmd = "blastx -num_threads $threads -db $db -query $query -evalue $evalue -out $outfile";
    my $opts->{source} = ["blast-$blast_version"];
    $opts->{jobname} = "BLASTX_$jobname";
    $opts->{threads} = $threads;
    $opts->{memory} = $memory;
    $opts->{log} = "$log_path/".basename($query).".$scheduler";
    
    my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
    print "      $outcmd\n";
    
    # Check for errors
    
    # Return results
    return ($outfile, $jobs);
}

sub run_rapsearch {
    my $self = shift;
    my $query = shift;
    my $alignment_base_dir = $self->{config}{alignment_base_dir};
    my $log_path = $self->{param}{log_path};
    my $db = $self->{config}{rapsearch_database};
    my $sample_id = $self->{param}{sample_id};
    my $queue = $self->{config}{queue};
    my $threads = $self->{config}{threads};
    my $rapsearch_version = $self->{config}{rapsearch_version};
    my $overwrite = $self->{param}{overwrite};
    my $memory = $self->{config}{rapsearch_memory};
    my $evalue = $self->{config}{rapsearch_evalue};
    my $scheduler = $self->{config}{jobsys};
    
    $log_path = $self->directory_check("$log_path/aligners");
    $log_path = $self->directory_check("$log_path/rapsearch");
    
    print "    Run RapSearch\n";
    $self->does_file_exist($query);
    
    # Get filename alone (no path)
    my $query_filename = basename($query);
    
    # Make an output directory for this aligner
    my $outdir = $self->directory_check("$alignment_base_dir/rapsearch");
    my $outfile = "$outdir/$query_filename.txt";
    
    unless ($overwrite) {
        if (-e $outfile) {
            print "    Found existing RapSearch output\n      $outfile\n    Skipping RapSearch!\n";
            return $outfile;
        }
        #else { print "--No existing file found; proceeding with RapSearch.\n"; }
    }
    
    my $jobname = basename($self->{param}{output_prefix});
    my $queryname = $query_filename;
    $queryname =~ s/[^0-9]//g;
    $jobname .= "_$queryname";
    
    my $cmd = "rapsearch -q $query -d $db -o $outfile -z $threads -e $evalue";
    my $opts->{source} = ["rapsearch-$rapsearch_version"];
    $opts->{jobname} = "RapSearch_$jobname";
    $opts->{threads} = $threads;
    $opts->{memory} = $memory;
    $opts->{log} = "$log_path/".basename($query).".$scheduler";
    
    my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
    print "      $outcmd\n";
    
    # Check for errors
    
    # Return results
    # Due to the way RapSearch outputs files, we need to append '.m8' to $outfile here in order to get a file that exists.
    $outfile .= ".m8";
    return ($outfile, $jobs);
}

sub run_diamond {
    my $self = shift;
    my $query = shift;
    my $alignment_base_dir = $self->{config}{alignment_base_dir};
    my $log_path = $self->{param}{log_path};
    my $db = $self->{config}{diamond_database};
    my $sample_id = $self->{param}{sample_id};
    my $queue = $self->{config}{queue};
    my $threads = $self->{config}{threads};
    my $diamond_version = $self->{config}{diamond_version};
    my $overwrite = $self->{param}{overwrite};
    my $memory = $self->{config}{diamond_memory};
    my $scheduler = $self->{config}{jobsys};
    
    $log_path = $self->directory_check("$log_path/aligners");
    $log_path = $self->directory_check("$log_path/diamond");
    
    print "    Run Diamond\n";
    $self->does_file_exist($query);
    
    # Get filename alone (no path)
    my $query_filename = basename($query);
    
    # Make an output directory for this aligner
    my $outdir = $self->directory_check("$alignment_base_dir/diamond");
    my $outfile = "$outdir/$query_filename.txt";
    
    # Diamond needs a temp directory as well, apparently
    my $tempdir = $self->directory_check("$alignment_base_dir/diamond/temp");
    
    unless ($overwrite) {
        if (-e $outfile) {
            print "    Found existing Diamond output\n      $outfile\n    Skipping Diamond!\n";
            return $outfile;
        }
        #else { print "--No existing file found; proceeding with Diamond.\n"; }
    }
    
    my $jobname = basename($self->{param}{output_prefix});
    my $queryname = $query_filename;
    $queryname =~ s/[^0-9]//g;
    $jobname .= "_$queryname";
    
    my $cmd = "diamond blastx -d $db -q $query -o $outfile -t $tempdir --threads $threads";
    my $opts->{source} = ["diamond-$diamond_version"];
    $opts->{jobname} = "Diamond_$jobname";
    $opts->{threads} = $threads;
    $opts->{memory} = $memory;
    $opts->{log} = "$log_path/".basename($query).".$scheduler";
    
    my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
    print "      $outcmd\n";
    
    # Check for errors
    
    # Return results
    return ($outfile, $jobs);
}

sub set_up_megan_command_file {
    # The 'command file' used by MEGAN seems to be dynamically generated. Do that here.
    # Cunningly, these commands include importing data, etc. 
    my $self = shift;
    my $commandfile = shift;
    my $alignment_file = shift;
    my $fasta_file = shift;
    my $megan_file = shift;
    my $aligner_used = $self->get_aligner_used($alignment_file);
    my $maxmatches = $self->{config}{megan_maxmatches};
    my $maxexpected = $self->{config}{megan_maxexpected};
    my $minsupport = $self->{config}{megan_minsupport};
    my $mincomplexity = $self->{config}{megan_mincomplexity};
    my $megan_taxafile = $self->{config}{megan_taxafile};
    
    if (ref $alignment_file eq 'ARRAY') { $alignment_file = join ', ', @$alignment_file; }
    if (ref $fasta_file eq 'ARRAY')     { $fasta_file = join ', ', @$fasta_file; }
    
    open(COMFILE, ">", $commandfile) or die "ERROR: Cannot open command file for aligner $aligner_used\n  $commandfile\n$!\n";
    if (($aligner_used eq 'diamond') || ($aligner_used eq 'rapsearch')) {
        print COMFILE "load taxGIFile='$megan_taxafile';\n";
        print COMFILE "import blastFile=$alignment_file fastaFile=$fasta_file meganFile=$megan_file maxMatches=$maxmatches maxExpected=$maxexpected minSupport=$minsupport minComplexity=$mincomplexity blastFormat=BlastTAB mapping='Taxonomy:GI_MAP=true';\n";
        print COMFILE "quit;";
    }
    if ($aligner_used =~ /blast/) {
        print COMFILE "import blastFile=$alignment_file fastaFile=$fasta_file meganFile=$megan_file maxMatches=$maxmatches maxExpected=$maxexpected minSupport=$minsupport minComplexity=$mincomplexity;\n";
        print COMFILE "quit;";
        
    }
    close COMFILE;
    $self->does_file_exist($commandfile);
}

sub run_megan {
    # Calls MEGAN on the supplied file
    # MEGAN wants an alignment file and a FASTA file. Also, a command file - a list of instructions, since you obviously can't mouse around menus in command-line mode. 
    # It creates a "meganfile"? Unsure, but I imagine it holds the data for generating plots etc.
    # Annoyingly, the commands are slightly different when different aligners are used. set_up_megan_command_file deals with that.
    # Import multiple alignment and FASTA files in a comma-separated list in the command file. Needs a bit of string modification here, though, and also some array
    # niftiness to identify which aligner is being used from the paths. 
    
    my $self = shift;
    my $alignment_file = shift;
    my $fasta_file = shift;
    my $output_dir = shift;
    my $log = shift;
    my $prerequisite_jobs = shift;
    my $megan_version = $self->{config}{megan_version};
    my $megan_license_file = $self->{config}{megan_license_file};
    #my $log_path = $self->{param}{log_path};
    my $queue = $self->{config}{queue};
    my $aligner_used = $self->get_aligner_used($alignment_file);
    my $memory = $self->{config}{megan_memory};
    my $display = $self->{config}{megan_display};
    
    #$self->does_file_exist($fasta_file);
    
    my $meganfile = "$output_dir/megan.rma";
    my $commandfile = "$output_dir/commands.cmds";
    $self->set_up_megan_command_file($commandfile, $alignment_file, $fasta_file, $meganfile);
    
    my $jobname = basename($self->{param}{output_prefix});
    my $cmd = "MEGAN -c $commandfile -L $megan_license_file --commandLineMode";
    my $opts->{source} = ["MEGAN-$megan_version"];
    #$opts->{depend} = $prerequisite_jobs;
    $opts->{jobname} = "MEGAN_$jobname";
    $opts->{threads} = 1;
    $opts->{memory} = $memory;
    $opts->{log} = $log;
    
    my ($outcmd, $jobs) = $self->submit_job($cmd, $opts);
    print "    MEGAN command:\n      $outcmd\n";
    
    $self->halt('megan');
    return $meganfile;
}

1;
