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

sub new {
    my $class = shift;
    my $self = bless {}, $class;
    return $self;
}

sub assign_parameters {
    # Adds the user-supplied parameters to the object, so I can get at them anywhere.
    my $self = shift;
    
    $self->{param}{config} = shift;
    $self->{param}{datadir} = shift;
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
    
    if (($self->{param}{exclude_all}) && (!$self->{param}{sub_start_size})) {
        die "ERROR: Option --exclude_all cannot be set if no subsampling size is specified!\nUse --subsample to specify a number of reads to use.\n";
    }
    
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
    
    # Print the aligners we'll use, because it doesn't seem thay're all firing right now
    if ($self->{param}{blastn}) { print "Use BLASTn\n"; }
    if ($self->{param}{blastx}) { print "Use BLASTx\n"; }
    if ($self->{param}{rapsearch}) { print "Use RapSearch\n"; }
    if ($self->{param}{diamond}) { print "Use Diamond\n"; }
    
}

sub read_config_file {
    # Reads the config file and adds the parameters within to the object, so I can get at them anywhere.
    my $self = shift;
    # These can all go in $self->{config}{name}
    # It looks, reasonably enough, as though there are distinct databases for each alignment tool. Make sure they're distinguished.
    # Also include version numbers of various software! This makes switching versions that bit easier.
    my $configfile = $self->{param}{config};
    
    open(CONFIG, "<", $configfile) or die "ERROR: Cannot open config file $configfile:\n$!\n";
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
        die "ERROR: File $file does not exist\n";
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
            die "ERROR: Cannot create output directory $dir\n";
        }
    }
    return $dir;
}

sub dechunk_and_unzip {
    my $self = shift;
    my $read = $_[0];
    my $files = $_[1];
    my $readfile = ();
    
    # Be aware that the correct file may already exist, amongst a set of chunked and/or compressed files.
    # Try to pick it out, if possible, before going any further.
    my @readfiles = ();
    FILE: foreach my $file (@$files) {
        chomp $file;
        if ($file =~ /_R$read\.fastq$/) { $readfile = $file; last FILE;  }
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
            if ($readfile =~ /fastq\.gz$/) {
                $readfile =~ s/fastq\.gz/fastq/;
                # Unzip it. See if I can keep the compressed version (use -c)
                # I don't think this needs to be bsubbed - this script itself should be bsubbed.
                `gunzip -c $readfile.gz 1> $readfile`;
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
            $readfile =~ s/fastq\.gz/fastq/;
            
            # One thing to double-check: does the file we're after already exist?
            # There's no point re-extracting data if it already exists. 
            unless (-e $readfile) {
                # Check if these files are compressed or not; it determines what we do with them.
                # If decompressing, we should only put compressed files in the list.
                my $compressed = 0;
                my @compressed_files = ();
                foreach my $file (@readfiles) {
                    if ($file =~ /fastq\.gz/) { $compressed = 1; push @compressed_files, $file; }
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
                print "Join chunked files:\n$call\n-----\n";
                `$call`;
            }
        }
    }
    return $readfile;
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
    
    my $reads = `wc -l $fastq_file`;
    
    # Error checking
    
    $reads = $reads / 4;
    return $reads;
}

sub get_fastq_read_ids {
    # Returns a big array of all the read IDs in a FastQ file. 
    my $self = shift;
    my $file = shift;
    
    $self->does_file_exist($file);
    my @ids = ();
    
    open INFILE, "<", $file or die "ERROR: Could not open fastq file $file: $!\n";
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
    # Try again, but keep it simpler
    # Make a consistent output structure so I always know what's happening
    # Try it in here first.
    # Chunking will remain a similar process, though it becomes a lot easier to see how it will work after looking at this.
    # Pro mode: make the file recognised as 'all', when requested, actually be the file that comes out of flash. Use a symlink?
    
    my $self = shift;
    my $start_size = $self->{param}{sub_start_size};
    my $step_size = $self->{param}{sub_step};
    my $num_subsamples = $self->{param}{sub_num};
    my $exclude_all = $self->{param}{exclude_all};
    my $overwrite = $self->{param}{overwrite};
    
    #my $sample_id = $self->{param}{sample_id};
    # Bear in mind that these last two inputs may not be supplied, or may be null.
    # In that case, make only one subsample.
    my $infile = shift;
    my $outdir = shift;
    
    my %subsample_files = ();
    
    # Set up the right directories - we want outdir/sample_id
    #my $sampledir = $self->directory_check("$outdir/$sample_id");
    
    # Get all the read IDs in the input file
    my $read_ids = $self->get_fastq_read_ids($infile);
    my $number_of_reads = @$read_ids;
    
    ## Get the numbers of reads to get in each subsample
    #my $subset_sizes = $self->list_subset_sizes();
    #
    ## Get a max of the number of reads expected to be used across our subsampling regimen.
    #my $max = max @$subset_sizes;
    #
    
    # Set up subsample numbers
    # Check that they don't overrun number of available reads; react appropriately if they do
    # Store the numbers in $self somewhere
    # Modify list_subset_sizes to simply recall that
    my %subsample_sizes = ();
    
    # Allow the option of excluding the 'all' setting - but only if no subsample start size is given.
    if (!$start_size) {
        if (!$exclude_all) { $subsample_sizes{all} = 1; }
    }
    #else {
    #    if ($start_size >= $number_of_reads) { push @subsample_sizes, 'all'; }
    #}
    
    my $subsize = $start_size;
    if ($num_subsamples) {
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
        if ($subsize < $number_of_reads) { $subsample_sizes{$subsize} = 1; }
        else {
            $subsample_sizes{all} = 1;
            print "WARN: Not enough reads to meet planned subsampling scheme!\n  Proposed sample:\t$subsize\n  Available:\t$number_of_reads\n";
        }
    }
    
    my @subset_sizes = keys %subsample_sizes;
    $self->{subset_sizes} = \@subset_sizes;
    
    # Assign read IDs to files, as necessary
    # Easiest way is probably to shuffle the array, then pick the first x ($start_size) IDs for the first subsample, then repeat for the next x + y for the second, and so on.
    # Store the subsample each sample has been assigned to in a hash for easy lookup shortly.
    # This is where I need to assign reasonable output structure.
    # The 'all' subset is excluded from this assignment process on the grounds that a file containing all reads is already present, at least in the form of a symlink.
    # I can also exclude (more precisely, decline to assign reads to) already-written files here unless the --overwrite parameter is supplied.
    my %subsamples = ();
    print "Making these subsets:\n";
    foreach my $subset_size (@subset_sizes) {
        print "    $subset_size\n";
        my $subset_file = ();
        if ($subset_size eq 'all') {
            # Make a symlink $infile->$outdir/all.fastq
            # If that doesn't work, copy-paste it.
            `ln -s $infile $outdir/all.fastq`;
            # If the 'all' subset is specified, add that link as a key to $subsample_files too. That will put it in the output file list.
            $subset_file = "$outdir/all.fastq";
            print "  Making 'all' subset (symlink from reads file)\n";
        }
        elsif ((!-e "$outdir/$subset_size.fastq") || ($overwrite)) {
            print "  Making '$subset_size' subset\n";
            @$read_ids = shuffle @$read_ids;
            my $i = 0;
            while ($i < $subset_size) {
                my $read_id = $read_ids->[$i];
                chomp $read_id;
                $subset_file = "$outdir/$subset_size.fastq";
                push @{$subsamples{$read_id}}, $subset_file;
                $i++;
                # Important: if $overwrite is set, we should delete the existing file, if present.
                if (($overwrite) && (-e $subset_file)) { `rm $subset_file`; }
            }
        }
        else {
            print "  Subset '$subset_size' already exists; using that\n";
            $subset_file = "$outdir/$subset_size.fastq";
        }
        $subsample_files{$subset_file} = 1;
    }
    my @uniq = keys %subsamples;
    
    # Now, I need to pick out the reads from the input file, and direct them to the subsample files they've been assigned to.
    # These are potentially very big files, so they should be read line by line. That would need a little read buffer type thing, since we want
    # chunks of 4 lines at a time.
    # Make sure to store the names of the subsample files created, so that they can be passed back. 
    # Of course, if no new subsample files are actually going to be created (@uniq is empty), we can skip this.
    my $lines = 0;
    if (@uniq) {
        open INFILE, "<", $infile or die "ERROR: Could not open fastq file $infile: $!\n";
        my @buffer = ();;  my $store_lines = 0;    my $outfile = ();   my $j = 0;
        while (my $line = <INFILE>) {
            # If line is a read ID, and has been assigned to a subsample, prepare to store it
            chomp $line;
            # Wow, that is some gnarly regex. It means a bunch of stuff begining with an @, then some stuff with semicolons at regular intervals.
            if ($line =~ /^\@.*:.*:.*:.*:.*:.*:.*/) {
                my $read_id = $line;
                if ($subsamples{$read_id}) {
                    $store_lines = 1;
                    push @buffer, $line;
                    # This sets the appropriate output paths for this read - the output path, the subsample ID we assigned earlier, and the source filename as a tail.
                    # $outfile is an array reference containing one or more of those.
                    $outfile = $subsamples{$read_id};
                }
            }
            elsif ($store_lines == 1) {
                push @buffer, $line;
                if (@buffer >= 4) {
                    # $outfile is an array reference to all the files this read should be printed to. Loop them, and print the stuff to each.
                    foreach my $file (@$outfile) {
                        $j++;
                        # If 4 lines have built up in the buffer, write them out to the right file and clear the buffer and counter
                        if (-e $file) { open (OUT, ">>", $file) or die "ERROR: Cannot open reads subsample file $file\n"; }
                        else          { open (OUT, ">", $file) or die "ERROR: Cannot open reads subsample file $file\n"; }
                        foreach my $l (@buffer) { print OUT "$l\n"; $lines ++; }
                        close OUT;
                    }
                    
                    $store_lines = 0;
                    @buffer = ();
                }
            }
        }
        close INFILE;
    }
    print "Wrote $lines lines to subset files\n";
    
    #$subsample_files{"remainder"} = ();
    my @subsample_files = keys %subsample_files;
    print "final subsample files:\n@subsample_files\n";
    return \@subsample_files;
}

sub rechunk {
    # This is meant to split read files back up into a set number of subdivisions, after the whole subsetting business has been handled. It will, hopefully, speed up our chronically slow aligners a bit.
    my $self = shift;
    my $readfile = shift;
    my $num_chunks = $self->{param}{num_chunks};
    my $output_prefix = $self->{param}{output_prefix};
    my $overwrite = $self->{param}{overwrite};
    
    # Get number of reads in the file
    my $read_ids = $self->get_fastq_read_ids($readfile);
    my $number_of_reads = @$read_ids;
    
    # Calculate a number; this will help us decide which chunk file each read should go into
    my $reads_per_chunk = int ($number_of_reads / $num_chunks);
    
    # Make the right output directory for the subset file
    my $dir = dirname($readfile);
    my $subset = basename($readfile);
    $subset =~ s/\.fastq//;
    
    my $chunkdir = $self->directory_check("$dir/$subset");
    
    # First, do a check; do these rechunked files already exist? We want to avoid repeating this is we can help it.
    # Criteria for skipping:
    # Must be correct number of rechunked files
    # All must have the correct number of reads (with a bit of variance)
    my $pattern = '*.fastq*';
    my $existing_chunk_files = $self->find_files($chunkdir,$pattern);
    
    print "Rechunking check: can we skip this step?\n";
    # $overwrite must not be set...
    unless ($overwrite) {
        print "  Overwrite not requested\n";
        # Files must exist...
        unless ($existing_chunk_files->[0] =~ /No such file or directory/) {
            print "  Files exist in output dir\n";
            # Number of files must match...
            if (@$existing_chunk_files == $num_chunks) {
                
                # Number of reads must (roughly) match...
                # Check them all
                my $var = int ($reads_per_chunk / 100);
                my $num_reads_is_fine = 1;
                foreach my $file (@$existing_chunk_files) {
                    my $read_ids = $self->get_fastq_read_ids($file);
                    my $number_of_subset_reads = @$read_ids;
                    unless (($number_of_subset_reads <= $reads_per_chunk + $var) && ($number_of_subset_reads >= $reads_per_chunk - $var)) {
                        $num_reads_is_fine = 0;
                        print "    File $file has wrong number of reads! ($number_of_subset_reads vs. ".($reads_per_chunk - $var)."-".($reads_per_chunk + $var).")\n";
                    }
                }
                
                # If that's all OK, feel free to skip.
                if ($num_reads_is_fine == 1) {
                    my @chunk_files = ();
                    foreach my $i (1..$num_chunks) {
                        push @chunk_files, "$chunkdir/$i.fastq";
                    }
                    return \@chunk_files;
                }
            }
            else { print "  Incorrect number of files in $chunkdir\n  (".@$existing_chunk_files." vs. $num_chunks)\n"; }
        }
        else { print "  No files exist in directory $chunkdir\n"; }
    }
    else { print "  Overwrite requested\n"; }
    
    # There is a minor complication here. Since this sub adds data to files mostly through appending, we should avoid writing to files that already exist.
    # That can be achieved simply by deleting the chunks directory and recreating it every time.
    if ($subset) {
        if (-d "$dir/$subset") { `rm -r $dir/$subset`; }
    }
    $chunkdir = $self->directory_check("$dir/$subset");
    
    
    # Now, I need to pick out the reads from the input file, assign them to a chunk file, and direct them there.
    # These are potentially very big files, so they should be read line by line. That would need a little read buffer type thing, since we want
    # chunks of 4 lines at a time.
    # Make sure to store the names of the chunk files created, so that they can be passed back. 
    my %chunk_files = ();
    open INFILE, "<", $readfile or die "ERROR: Could not open fastq file $readfile: $!\n";
    my @buffer = ();
    my $store_lines = 0;    my $read_count = 0; my $outfile = ();
    while (my $line = <INFILE>) {
        chomp $line;
        # If line is a read ID, and has been assigned to a subsample, prepare to store it
        if ($line =~ /^\@.*:.*:.*:.*:.*:.*:.*/) {
            my $read_id = $line;
            chomp $read_id;
            $store_lines = 1;
            push @buffer, $line;
            
            # This works out which chunk to assign this read to. It will split the file up pretty much evenly.
            # I increment $chunk because it will be inted down to 0 at first. The unless condition stops $chunk going to $num_chunks + 1 at the last read.
            my $chunk = int($read_count / $reads_per_chunk);
            $chunk ++;
            if ($chunk > $num_chunks) { $chunk = $num_chunks; }
            
            $outfile = "$chunkdir/$chunk.fastq";
            $chunk_files{$outfile} = 1;
        }
        elsif ($store_lines == 1) {
            push @buffer, $line;
            if (@buffer >= 4) {
                # If 4 lines have built up in the buffer, write them out to the right file and clear the buffer and counter
                if (-e $outfile) { open (OUT, ">>", $outfile) or die "ERROR: Cannot open reads subsample file $outfile\n"; }
                else             { open (OUT, ">", $outfile) or die "ERROR: Cannot open reads subsample file $outfile\n"; }
                foreach my $l (@buffer) { print OUT "$l\n"; }
                close OUT;
                
                $store_lines = 0;
                @buffer = ();
                $read_count ++;
            }
        }
    }
    close INFILE;
    
    my @chunk_files = keys %chunk_files;
    return \@chunk_files;
}

sub extract_list_of_jobs {
    # When given a load of STDOUT data, this will return a list of LSF job IDs that were started within that session.
    my $self = shift;
    my $ramble = shift;
    
    my @lines = split /\n/, $ramble;
    my @jobs = ();
    foreach my $line (@lines) {
        if ($line =~ /Job \<[0-9]+\> is submitted to queue /) {
            my @sp = split /[\<\>]/, $line;
            push @jobs, $sp[1];
        }
    }
    
    if (@jobs == 0) {
        die "PROBABLE ERROR: No LSF jobs found for the above step.\n";
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
    
    # Check the list of current jobs every minute or so, until all the input jobs have disappeared
    # from the list.
    my @intersection = ();
    my $checks = ();
    do {
        $checks ++;
        # Get running job IDs using bjobs and awk
        my $running_jobs = `bjobs | awk '{print \$1}'`;
        my @running_jobs = split /\n/, $running_jobs;
        my $line = shift @running_jobs;
        
        # Compare using List::Compare
        # BUT only if LSF has indeed returned a list of jobs!
        unless ($line =~ /Please wait/) {
            my $lc = List::Compare->new($check_jobs, \@running_jobs);
            @intersection = $lc->get_intersection;
        }
        
        # If any jobs in the checklist are still active, then wait a short time before checking again.
        # No point checking ten billion times per minute.
        if (@intersection >= 1) {
            if ($checks <= 10)    { sleep 5; }
            elsif ($checks <= 20) { sleep 15; }
            else                  { sleep 30; }
        }
    }
    until (@intersection == 0);
    
    # Wait another 60 seconds; perl/the cluster apparently needs this in order to sort out newly created files.
    sleep 60;
    
    print "All pending jobs for this step have completed.\n";
    return "Done";
}

sub run_nextclip {
    my $self = shift;
    my $read1file = shift;
    my $read2file = shift;
    my $output_path = shift;
    my $output_path_slash =$output_path;
    $output_path_slash =~ s!/*$!/!; # Add a trailing slash if none is present
    my $log_path = $self->{param}{log_path};
    my $queue = $self->{config}{queue};
    my $run_nextclip = $self->{config}{run_nextclip};
    my $nextclip_version = $self->{config}{nextclip_version};
    my $remove_pcr_duplicates = $self->{config}{remove_pcr_duplicates};
    
    # Are we actually going to run NextClip? If not, just return the input files and don't do anything.
    unless ($run_nextclip eq 'yes') {
        return ($read1file, $read2file);
    }
    
    # Set up the job
    my $jobname = basename($self->{param}{output_prefix});
    my $bsub = "source nextclip-$nextclip_version; bsub -J NextClip_$jobname -q $queue -oo $log_path/nextclip_run_log.lsf \"nextclip -i $read1file -j $read2file -o $output_path_slash"."results \" ";
    if ($remove_pcr_duplicates eq 'yes') {
        $bsub .= "--remove_duplicates ";
    }
    
    # Run it and wait for it to finish
    print "NextClip command:\n$bsub\n";
    my $r1jobs = `$bsub`;
    
    my $jobs = $self->extract_list_of_jobs($r1jobs);
    my $done = $self->done_when_its_done($jobs);
    
    # Note on removal of PCR duplicates:
    # Nextclip dumps all its output into the place specified in the output prefix (obviously). The actual read files are named results_[A,B,C,D]_R[1,2], though. The majority of the reads will usually end up in category D, because we usually won't be doing metagenomics on Nextera long mate pair libraries.
    # We want all of those reads, so we have to use cat and pipe them off to a single file, or something. Nonetheless, we should then have the full set of reads, minus PCR duplicates. (If we want a specific set of reads pulled out by NextClip in the future, this is the best place to get it).
    # Try to keep the basename of the original input files, just in case anything else needs it.
    my $read1file_out = "$output_path/".basename($read1file);
    my $read2file_out = "$output_path/".basename($read2file);
    print "Nextclip reads output files:\n$read1file_out\n$read2file_out\n";
    $bsub = "bsub -J CatR1 -oo /dev/null \"cat $output_path/results_*_R1.fastq > $read1file_out \" ";
    $r1jobs = `$bsub`;
    $jobs = $self->extract_list_of_jobs($r1jobs);
    $done = $self->done_when_its_done($jobs);
    $bsub = "bsub -J CatR2 -oo /dev/null \"cat $output_path/results_*_R2.fastq > $read2file_out \" ";
    $r1jobs = `$bsub`;
    $jobs = $self->extract_list_of_jobs($r1jobs);
    $done = $self->done_when_its_done($jobs);
    return ($read1file_out, $read2file_out);
}

sub run_fastqc {
    my $self = shift;
    my $read1file = shift;
    my $read2file = shift;
    my $output_path = shift;
    my $log_path = $self->{param}{log_path};
    my $queue = $self->{config}{queue};
    my $fastqc_version = $self->{config}{fastqc_version};
    my $overwrite = $self->{param}{overwrite};
    my $halt_at = $self->{param}{halt_at};
    
    if ($halt_at) {
        if ($halt_at eq 'fastqc') { die "Pipeline halted at FastQC\n"; }
    }
    
    my $fastqcfileR1 = basename($read1file);
    $fastqcfileR1 =~ s/\.fastq/_fastqc\.zip/g;
    my $fastqcfileR2 = basename($read2file);
    $fastqcfileR2 =~ s/\.fastq/_fastqc\.zip/g;
    
    # How can I check if this step has been done already?
    # Look for output files in the output path, of course.
    unless ($overwrite) {
        my $files = $self->find_files($output_path, '*fastqc.html');
        if (@$files > 0) {
            print "--Found FastQC files\n@$files\nSkipping FastQC!\n";
            return  ("$output_path/$fastqcfileR1", "$output_path/$fastqcfileR2");
        }
        else { print "--No existing file found; proceeding with trimming.\n"; }
    }
    
    my $jobname = basename($self->{param}{output_prefix});
    my $bsub = "source fastqc-$fastqc_version; bsub -J FastQC_$jobname -q $queue -oo $log_path/fastqc_R1.lsf \"fastqc $read1file --outdir $output_path\" ";
    print "FASTQC commands:\n$bsub\n";
    my $r1jobs = `$bsub`;
    my $r2jobs = $r1jobs;
    $bsub = "source fastqc-0.11.2; bsub -q $queue -oo $log_path/fastqc_R2.lsf \"fastqc $read2file --outdir $output_path\" ";
    print "$bsub\n";
    $r2jobs .= `$bsub`;
    # (That concatenates all the jobs into a single string)
    
    my $jobs = $self->extract_list_of_jobs($r2jobs);
    my $done = $self->done_when_its_done($jobs);
    
    # Check for errors
    
    return ("$output_path/$fastqcfileR1", "$output_path/$fastqcfileR2");
}

sub examine_fastqc_results {
    my $self = shift;
    my $fastqc_results_file = shift;
    my $trimming_done = shift;
    
    # Unzip $fastqc_results_file
    # Read fastqc_data.txt in directory created from that
    # Locate the adapter content section of that file
    # Prepare warnings, if necessary, based on what's in there
    my $fastqc_results_dir = $self->unzip_fastqc_results($fastqc_results_file);
    
    # Get the name of that directory (it'll be used for printing output to the user)
    my @dirsplit = split /\//, $fastqc_results_dir;
    my $results_dirname = $dirsplit[-1];
    
    open(RESULTS, "<", "$fastqc_results_dir/fastqc_data.txt") or die "ERROR: Cannot open unzipped FastQC data\n$fastqc_results_dir/fastqc_data.txt\n";
    my @data = <RESULTS>;
    close RESULTS;
    
    my $line = ();
    do {
        $line = shift @data;
    }
    until ($line =~ /Adapter Content/);
    chomp $line;
    my @linesplit = split /\t/, $line;
    my $test_status = $linesplit[1];
    chomp $test_status;
    if ($test_status =~ /pass/) {
        print "FastQC test: adapter content for dataset $results_dirname OK\n";
    }
    else {
        if ($trimming_done) {
            # If trimming has been done and adapter content is still too high, throw a big wobbly.
            print "ERROR: Adapter content too high after trimming in dataset $results_dirname!\n";
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
            print "WARN: Adapter content too high (before trimming) in dataset $results_dirname\n";
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
    my $read1file = shift;
    my $read2file = shift;
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
    my $halt_at = $self->{param}{halt_at};
    my $memory = $self->{config}{trimming_memory};
    
    if ($halt_at) {
        if ($halt_at eq 'trimming') { die "Pipeline halted at trimming\n"; }
    }
    
    $self->does_file_exist($read1file);
    $self->does_file_exist($read2file);
    
    my ($read1file_trimmed, $read2file_trimmed) = ();
    if ($run_trimming =~ /yes/) {
        $read1file_trimmed = "$trimming_output_path/".basename($read1file);
        $read2file_trimmed = "$trimming_output_path/".basename($read2file);
        $read1file_trimmed =~ s/\.fastq/_trimmed.fastq/;
        $read2file_trimmed =~ s/\.fastq/_trimmed.fastq/;
        
        my $read1file_trimmed_single = $read1file_trimmed;
        my $read2file_trimmed_single = $read2file_trimmed;
        $read1file_trimmed_single =~ s/\.fastq/_single.fastq/;
        $read2file_trimmed_single =~ s/\.fastq/_single.fastq/;
        
        unless ($overwrite) {
            if ((-e $read1file_trimmed) && (-e $read2file_trimmed)) {
                print "--Found \n$read1file_trimmed \nand \n$read2file_trimmed\nSkipping trimming!\n";
                return ($read1file_trimmed, $read2file_trimmed);
            }
            else { print "--No existing file found; proceeding with trimming.\n"; }
        }
        
        my $hot_air = "-R \"rusage[mem=$memory] span[hosts=1]\"";
        # I don't know what the :2:30:10 bit after $adaptersfile does, but it seems to be important.
        my $jobname = basename($self->{param}{output_prefix});
        my $bsub = "source trimmomatic-$trimmomatic_version; bsub -J Trimming_$jobname -q $queue -oo $log_path/trim.lsf $hot_air -n $threads \"java -jar /tgac/software/testing/trimmomatic/0.30/x86_64/bin/trimmomatic-0.30.jar PE -phred33 -threads $threads $read1file $read2file $read1file_trimmed $read1file_trimmed_single $read2file_trimmed $read2file_trimmed_single ILLUMINACLIP:$adaptersfile:$seed_mismatches:$palindrome_clip_threshold:$simple_clip_threshold SLIDINGWINDOW:$trimming_sliding_window MINLEN:$trimming_min_length\" ";
        print "TRIMMOMATIC command:\n$bsub\n";
        my $r1jobs = `$bsub`;
        
        my $jobs = $self->extract_list_of_jobs($r1jobs);
        my $done = $self->done_when_its_done($jobs);
        
        # Check for errors
    }
    else {
        print "    Config file requests no trimming; skipping this step\n";
        $read1file_trimmed = $read1file;
        $read2file_trimmed = $read2file;
    }
    
    return ($read1file_trimmed, $read2file_trimmed);
}

sub run_kontaminant {
    my $self = shift;
    my $read1file = shift;
    my $read2file = shift;
    my $filtering_output_dir = shift;
    my $reference = shift;
    my $run_filtering = $self->{config}{run_filtering};
    my $log_path = $self->{param}{log_path};
    my $database =$self->{config}{kontaminant_database};
    my $queue = $self->{config}{queue};
    my $kontaminant_version = $self->{config}{kontaminant_version};
    my $overwrite = $self->{param}{overwrite};
    my $halt_at = $self->{param}{halt_at};
    my $mem_width = $self->{config}{kontaminant_mem_width};
    my $mem_height = $self->{config}{kontaminant_mem_height};
    my $memory = $self->{config}{filtering_memory};
    
    if ($halt_at) {
        if ($halt_at eq 'filtering') { die "Pipeline halted at filtering\n"; }
    }
    
    $self->does_file_exist($read1file);
    $self->does_file_exist($read2file);
    
    my ($r1outfile, $r2outfile) = ();
    if ($run_filtering =~ /yes/) {
        $self->directory_check("$log_path/kontaminant");
        
        $r1outfile = "$filtering_output_dir/filtered_".basename($read1file);
        $r2outfile = "$filtering_output_dir/filtered_".basename($read2file);
        
        unless ($overwrite) {
            if ((-e $r1outfile) && (-e $r2outfile)) {
                print "--Found \n$r1outfile \nand \n$r2outfile\nSkipping filtering!\n";
                return ($r1outfile, $r2outfile);
            }
            else { print "--No existing file found; proceeding with filtering.\n"; }
        }
        
        my $jobname = basename($self->{param}{output_prefix});
        my $bsub = "source kontaminant-$kontaminant_version; bsub -J Filter_$jobname -R \"rusage[mem=$memory]\" -q $queue -oo $log_path/kontaminant.lsf \"kontaminant -f -1 $read1file -2 $read2file -c $reference -d $database -k 21 -o $filtering_output_dir/filtered_ -r $filtering_output_dir/removed_ -p $log_path/kontaminant -n $mem_height -b $mem_width\" ";
        print "KONTAMINANT command:\n$bsub\n";
        my $r1jobs = `$bsub`;
        
        my $jobs = $self->extract_list_of_jobs($r1jobs);
        my $done = $self->done_when_its_done($jobs);
        
        # Check for errors
    }
    else {
        print "    Config file requests no filtering; skipping this step\n";
        $r1outfile = $read1file;
        $r2outfile = $read2file;
    }
    
    return ($r1outfile, $r2outfile);
}

sub run_flash {
    my $self = shift;
    my $read1file = shift;
    my $read2file = shift;
    my $output_dir = shift;
    my $output_file = shift;
    my $log_path = $self->{param}{log_path};
    my $queue = $self->{config}{queue};
    my $flash_version = $self->{config}{flash_version};
    my $overwrite = $self->{param}{overwrite};
    my $halt_at = $self->{param}{halt_at};
    
    if ($halt_at) {
        if ($halt_at eq 'flash') { die "Pipeline halted at FLASH\n"; }
    }
    
    $self->does_file_exist($read1file);
    $self->does_file_exist($read2file);
    
    # The string in $output_file isn't actually the full name of the output file. FLASH creates it by adding this:
    my $full_output_file = "$output_file.extendedFrags.fastq";
    
    unless ($overwrite) {
        if (-e "$output_dir/$full_output_file") {
            print "--Found \n$output_dir/$full_output_file \nSkipping FLASH!\n";
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
    my $bsub = "source flash-$flash_version; bsub -J Flash_$jobname -q $queue -o $log_path/flash.lsf \"flash $read1file $read2file -M 150 -o $output_file\" ";
    print "FLASH command:\n$bsub\n";
    my $r1jobs = `$bsub`;
    
    my $jobs = $self->extract_list_of_jobs($r1jobs);
    my $done = $self->done_when_its_done($jobs);
    
    # Check for errors
    # Output files should have the value supplied as $output_file, but with a couple of different extensions indicating what they are.
    
    # Remember to change back to original working directory.
    chdir($pwd);
    print "Changed working directory back to\n$pwd\n";
    
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
    my $fastx_version = $self->{config}{fastx_toolkit_version};
    my $overwrite = $self->{param}{overwrite};
    
    $self->does_file_exist($fastq_filepath);
    
    my $fastq_filename = basename($fastq_filepath);
    my $fasta_filename = $fastq_filename;
    $fasta_filename =~ s/\.fastq/\.fasta/;
    my $fasta_filepath = "$outdir/$fasta_filename";
    
    unless ($overwrite) {
        if (-e $fasta_filepath) {
            print "--Found \n$fasta_filepath \nSkipping FASTQ->FASTA coversion!\n";
            return $fasta_filepath;
        }
        else { print "--No existing file found; proceeding with FASTQ->FASTA conversion.\n"; }
    }
    
    #my $bsub = "source fastx_toolkit-$fastx_version; bsub -q $queue -oo $log_path/convert_fasta.lsf \"fastq_to_fasta -Q 33 -n -i $fastq_filepath -o $fasta_filepath\" ";
    #print "CONVERT FASTQ TO FASTA command:\n$bsub\n";
    #my $r1jobs = `$bsub`;
    #
    #my $jobs = $self->extract_list_of_jobs($r1jobs);
    ##my $done = $self->done_when_its_done($jobs);
    #
    ## Check for errors
    #
    #return ($fasta_filepath, $jobs);
    
    # First, count read IDs so we can check the right number were added
    #my $read_ids = $self->get_fastq_read_ids($fastq_filepath);
    #my $number_of_reads = @$read_ids;
    
    open INFILE, "<", $fastq_filepath or die "ERROR: Could not open fastq file $fastq_filepath: $!\n";
    open (OUT, ">", $fasta_filepath) or die "ERROR: Cannot open FASTA output file $fasta_filepath: $!\n";
    my @buffer = ();
    my $read_count = 0; my $outfile = ();
    while (my $line = <INFILE>) {
        chomp $line;
        # Keep track of previous lines; that will help us keep pace. (Store them in an array, and pop n' push them down).
        # This is a kind of sliding window approach. When the window can be determined to be in the right place, dump the buffer to the output file.
        push @buffer, $line;
        if (@buffer == 4) {
            my $header = $buffer[0];
            my $seq = $buffer[1];
            $header =~ s/^\@/\>/;
            print OUT "$header\n$seq\n"; 
            $read_count ++;
            @buffer = ();
        }
    }
    close INFILE;
    close OUT;
    
    #unless ($read_count == $number_of_reads) {
    #    print "WARN: Number of reads in FASTA file ($read_count) does not match number in input FASTQ! ($number_of_reads)\n";
    #}
    
    return ($fasta_filepath);
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
    my $halt_at = $self->{param}{halt_at};
    
    if ($halt_at) {
        if ($halt_at eq 'alignment') { die "Pipeline halted at alignment\n"; }
    }
    
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
    
    #if (@jobs > 0) { my $done = $self->done_when_its_done(\@jobs); }
    
    return (\@results, \%jobs);
}

sub run_blastn {
    my $self = shift;
    my $query = shift;
    my $alignment_base_dir = $self->{config}{alignment_base_dir};
    #my $log_path = $self->{param}{log_path};
    my $log_path = dirname($query);
    my $db = $self->{config}{blastn_database};
    my $sample_id = $self->{param}{sample_id};
    my $queue = $self->{config}{queue};
    my $threads = $self->{config}{threads};
    my $blast_version = $self->{config}{blast_version};
    my $overwrite = $self->{param}{overwrite};
    my $memory = $self->{config}{blast_memory};
    my $evalue = $self->{config}{blast_evalue};
    
    print "Run BLASTn\n";
    $self->does_file_exist($query);
    
    # Get filename alone (no path)
    my $query_filename = basename($query);
    
    # Make an output directory for this aligner
    my $outdir = $self->directory_check("$alignment_base_dir/blastn");
    my $outfile = "$outdir/$query_filename.txt";
    
    unless ($overwrite) {
        if (-e $outfile) {
            print "--Found \n$outfile\nSkipping BLASTn!\n";
            return $outfile;
        }
        else { print "--No existing file found; proceeding with BLASTn.\n"; }
    }
    
    my $jobname = basename($self->{param}{output_prefix});
    my $queryname = $query_filename;
    $queryname =~ s/[^0-9]//g;
    $jobname .= "_$queryname";
    my $command = "blastn -num_threads $threads -db $db -query $query -evalue $evalue -out $outfile";
    my $bsub = "source blast-$blast_version; bsub -J BLASTN_$jobname -n $threads -q $queue -R \"rusage[mem=$memory] span[hosts=1]\" -oo $log_path/blastn.lsf \"$command\" ";
    print "BLASTN command:\n$bsub\n";
    my $r1jobs = `$bsub`;

    my $jobs = $self->extract_list_of_jobs($r1jobs);
    
    # Check for errors
    
    # Return results
    return ($outfile, $jobs);
}

sub run_blastx {
    my $self = shift;
    my $query = shift;
    my $alignment_base_dir = $self->{config}{alignment_base_dir};
    #my $log_path = $self->{param}{log_path};
    my $log_path = dirname($query);
    my $db = $self->{config}{blastx_database};
    my $sample_id = $self->{param}{sample_id};
    my $queue = $self->{config}{queue};
    my $threads = $self->{config}{threads};
    my $blast_version = $self->{config}{blast_version};
    my $overwrite = $self->{param}{overwrite};
    my $memory = $self->{config}{blast_memory};
    my $evalue = $self->{config}{blast_evalue};
    
    print "Run BLASTx\n";
    $self->does_file_exist($query);
    
    # Get filename alone (no path)
    my $query_filename = basename($query);
    
    # Make an output directory for this aligner
    my $outdir = $self->directory_check("$alignment_base_dir/blastx");
    my $outfile = "$outdir/$query_filename.txt";
    
    unless ($overwrite) {
        if (-e $outfile) {
            print "--Found \n$outfile\nSkipping BLASTx!\n";
            return $outfile;
        }
        else { print "--No existing file found; proceeding with BLASTx.\n"; }
    }
    
    my $jobname = basename($self->{param}{output_prefix});
    my $queryname = $query_filename;
    $queryname =~ s/[^0-9]//g;
    $jobname .= "_$queryname";
    my $command = "blastx -num_threads $threads -db $db -query $query -evalue $evalue -out $outfile";
    my $bsub = "source blast-$blast_version; bsub -J BLASTX_$jobname -n $threads -q $queue -R \"rusage[mem=$memory] span[hosts=1]\" -oo $log_path/blastx.lsf \"$command\" ";
    print "BLASTX command:\n$bsub\n";
    my $r1jobs = `$bsub`;

    my $jobs = $self->extract_list_of_jobs($r1jobs);
    
    # Check for errors
    
    # Return results
    return ($outfile, $jobs);
}

sub run_rapsearch {
    my $self = shift;
    my $query = shift;
    my $alignment_base_dir = $self->{config}{alignment_base_dir};
    #my $log_path = $self->{param}{log_path};
    my $log_path = dirname($query);
    my $db = $self->{config}{rapsearch_database};
    my $sample_id = $self->{param}{sample_id};
    my $queue = $self->{config}{queue};
    my $threads = $self->{config}{threads};
    my $rapsearch_version = $self->{config}{rapsearch_version};
    my $overwrite = $self->{param}{overwrite};
    my $memory = $self->{config}{rapsearch_memory};
    my $evalue = $self->{config}{rapsearch_evalue};
    
    print "Run RapSearch\n";
    $self->does_file_exist($query);
    
    # Get filename alone (no path)
    my $query_filename = basename($query);
    
    # Make an output directory for this aligner
    my $outdir = $self->directory_check("$alignment_base_dir/rapsearch");
    my $outfile = "$outdir/$query_filename.txt";
    
    unless ($overwrite) {
        if (-e $outfile) {
            print "--Found \n$outfile\nSkipping RapSearch!\n";
            return $outfile;
        }
        else { print "--No existing file found; proceeding with RapSearch.\n"; }
    }
    
    my $jobname = basename($self->{param}{output_prefix});
    my $queryname = $query_filename;
    $queryname =~ s/[^0-9]//g;
    $jobname .= "_$queryname";
    my $command = "rapsearch -q $query -d $db -o $outfile -z $threads -e $evalue";
    my $bsub = "source rapsearch-$rapsearch_version; bsub -J RapSearch_$jobname -n $threads -q $queue -R \"rusage[mem=$memory] span[hosts=1]\" -oo $log_path/rapsearch.lsf \"$command\" ";
    print "RAPSEARCH command:\n$bsub\n";
    my $r1jobs = `$bsub`;
    
    my $jobs = $self->extract_list_of_jobs($r1jobs);
    
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
    #my $log_path = $self->{param}{log_path};
    my $log_path = dirname($query);
    my $db = $self->{config}{diamond_database};
    my $sample_id = $self->{param}{sample_id};
    my $queue = $self->{config}{queue};
    my $threads = $self->{config}{threads};
    my $diamond_version = $self->{config}{diamond_version};
    my $overwrite = $self->{param}{overwrite};
    my $memory = $self->{config}{diamond_memory};
    
    print "Run Diamond\n";
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
            print "--Found \n$outfile\nSkipping Diamond!\n";
            return $outfile;
        }
        else { print "--No existing file found; proceeding with Diamond.\n"; }
    }
    
    my $jobname = basename($self->{param}{output_prefix});
    my $queryname = $query_filename;
    $queryname =~ s/[^0-9]//g;
    $jobname .= "_$queryname";
    my $command = "diamond blastx -d $db -q $query -o $outfile -t $tempdir --threads $threads";
    my $bsub = "source diamond-$diamond_version; bsub -J Diamond_$jobname -n $threads -q $queue -R \"rusage[mem=$memory] span[hosts=1]\" -oo $log_path/diamond.lsf \"$command\" ";
    print "DIAMOND command:\n$bsub\n";
    my $r1jobs = `$bsub`;
    
    my $jobs = $self->extract_list_of_jobs($r1jobs);
    
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
    
    open(COMFILE, ">", $commandfile) or die "ERROR: Cannot open command file for aligner $aligner_used\n: $!\n";
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
}

sub run_megan {
    # Calls MEGAN on the supplied file
    # MEGAN wants an alignment file and a FASTA file. Also, a command file - a list of instructions, since you obviously can't mouse around menus in command-line mode. 
    # It creates a "meganfile"? Unsure, but I imagine it holds the data for generating plots etc.
    # Annoyingly, the commands are slightly different when different aligners are used. set_up_megan_command_file deals with that.
    # Import multiple alignment and FASTA files in a comma-separated list in the command file. Needs a bit of string modification here, though, and also some array niftiness to identify which aligner is being used from the paths. 
    
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
    my $halt_at = $self->{param}{halt_at};
    my $memory = $self->{config}{megan_memory};
    my $display = $self->{config}{megan_display};
    
    if ($halt_at) {
        if ($halt_at eq 'megan') { die "Pipeline halted at MEGAN\n"; }
    }
    
    #$self->does_file_exist($fasta_file);
    
    my $meganfile = "$output_dir/megan.rma";
    my $commandfile = "$output_dir/commands.cmds";
    $self->set_up_megan_command_file($commandfile, $alignment_file, $fasta_file, $meganfile);
    
    # Job dependency string
    
    my $dependency_string = ();
    if (ref $prerequisite_jobs eq 'ARRAY') {
        my @dep = ();
        foreach my $job (@$prerequisite_jobs) { push @dep, "ended($job)"; }
        my $dependency = join " && ", @dep;
        $dependency_string = " -w \"$dependency\" ";
    }
    elsif ($prerequisite_jobs) {
        $dependency_string = " -w \"$prerequisite_jobs\" ";
    }
    else {
        $dependency_string = " ";
    }
    
    my $jobname = basename($self->{param}{output_prefix});
    my $bsub = "source MEGAN-$megan_version; bsub".$dependency_string."-q $queue -J MEGAN_$jobname -R \"rusage[mem=$memory]\" -oo $log \"export DISPLAY=$display; MEGAN -g -c $commandfile -L $megan_license_file\" ";
    print "MEGAN command:\n$bsub\n";
    my $r1jobs = `$bsub`;
    my $jobs = $self->extract_list_of_jobs($r1jobs);
    my $done = $self->done_when_its_done($jobs);
    
    return $meganfile;
}

sub calculate_max_diversity {
    # Over a set of subsamples of varying size, species diversity will tend towards an asymptote. This sub predicts the number of species likely to exist based on that curve.
    # In future, I would like to expand this to estimate x depth of coverage required to get y depth of coverage on every species (so, y depth of least abundant species?)
    # This is kind of a tough one. I might only be able to get this information from MEGAN, and I don't know how to do that yet. Keep it in mind anyway.
    
}

1;
