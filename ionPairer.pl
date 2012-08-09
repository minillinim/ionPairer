#!/usr/bin/env perl
###############################################################################
#
#   ionPairer.pl
#
#   Work out paired mappings based on information in a sam file
#
#   Copyright (C) Michael Imelfort, Ben Woodcroft
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;
use Carp;

#CPAN modules
use File::Spec;
use File::Basename;
use Bio::SeqIO;
use Data::Dumper;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
my $global_options = checkParams();
if(!exists $global_options->{'silent'}) {
    printAtStart();
}
######################################################################
# CODE HERE
######################################################################
# get a working dir
my $global_working_dir = File::Spec->rel2abs($global_options->{'working_dir'});
unless(-e $global_working_dir or mkdir $global_working_dir) {
    croak sprintf("Unable to create output directory %s: %s\n", $global_working_dir, $!);
}

# globals
my %global_con_2_int = ();         # contigIDs to integers
my %global_int_2_con = ();         # integers to contigIDs
my $global_conInt = 1;             # UID
my %global_con_2_len = ();         # contig integers versus lengths
my %global_reads_2_map = ();       # reads vs mapping info
my $global_minimum_links = overrideDefault(3, 'min_links'); # minimum number of links accepted during graph creation
my $reference_fasta = $global_options->{'reference_fasta'};

# get output file names and handles
my ($file_root, undef, undef) = fileparse($global_options->{'sam1'});

my $unpaired_file = File::Spec->catfile($global_options->{'working_dir'}, $file_root.".unpaired.csv");
my $pcr_duplicate_file = File::Spec->catfile($global_options->{'working_dir'}, $file_root.".pcr_duplicates.csv");
my $error_paired_file = File::Spec->catfile($global_options->{'working_dir'}, $file_root.".error_paired.csv");
my $paired_file = File::Spec->catfile($global_options->{'working_dir'}, $file_root.".paired.csv");
my $short_links_file = File::Spec->catfile($global_options->{'working_dir'}, $file_root.".short_links.csv");
my $all_links_file = File::Spec->catfile($global_options->{'working_dir'}, $file_root.".unique_links.csv");

my $pcr_duplicate_fh = openWrite($pcr_duplicate_file);
my $unpaired_fh = openWrite($unpaired_file);
my $error_paired_fh = openWrite($error_paired_file);
my $paired_fh = openWrite($paired_file);
my $short_links_fh = openWrite($short_links_file);
my $all_links_fh = openWrite($all_links_file);

# make headers # Array indicies: conID, pos, prob, strand, conID, pos, prob, strand
my $delimeter = "\t";
my $unpaired_header = "#".join($delimeter, ("contig_name","contig_length","read_id","position","strand"))."\n";
my $paired_header = "#".join($delimeter, ("contig_name","contig_length","read_id","position1","strand1","position2","strand2"))."\n";
my $all_links_header = "#".join($delimeter, ("contig1_name","contig1_position","contig1_probability","contig1_strand","contig2_name","contig2_position","contig2_probability","contig2_strand"))."\n";
my $links_header = "#".join($delimeter, ("contig1_name","contig1_length","contig2_name","contig2_length","read_id","contig1_position","contig1_strand","contig2_position","contig2_strand"))."\n";
my $pcr_header = "#".join($delimeter, ("contig1_name","contig1_position","contig1_strand","contig2_name","contig2_position","contig2_strand"))."\n";

print $unpaired_fh $unpaired_header;
print $paired_fh $paired_header;
print $error_paired_fh $paired_header;
print $all_links_fh $all_links_header;
print $short_links_fh $links_header;
print $pcr_duplicate_fh $pcr_header;

# parse the sam files
my @samfiles = ($global_options->{'sam1'},$global_options->{'sam2'});
my $first_sam = 1;
foreach my $sam_fn (@samfiles) {
    if(!exists $global_options->{'silent'}) { print "Parsing: $sam_fn\n"; }
    open my $sam_fh, "<", $sam_fn or die "**ERROR: could not open sam file $!\n";
    while(<$sam_fh>) {
        chomp $_;
        if($_ =~ /^@/) {
            # work out and store contig lengths
            next if(1 != $first_sam);
            my @fields = split(/:/, $_);
            $fields[1] =~ s/\tLN//;
            my $con_id = $fields[1];
            my $con_int = 0;
            if(!exists $global_con_2_int{$con_id}) {
                $con_int = $global_conInt;
                $global_conInt++;
                $global_con_2_int{$con_id} = $con_int;
                $global_int_2_con{$con_int} = $con_id;
            } else {
                $con_int = $global_con_2_int{$con_id};
            }
            $global_con_2_len{$con_int} = $fields[2];         
        } else {
            my @sam_fields = split(/\t/, $_);
            my @mapping_flags = split(//, dec2bin($sam_fields[1]));
            
            # make sure it mapped
            next if($mapping_flags[5] eq "1");
            
            # get a contig ID (int)
            my $con_id = $sam_fields[2];
            if(!exists $global_con_2_int{$con_id}) {
                croak "Unknown contig ID: $con_id\n";
            }
            $con_id = $global_con_2_int{$con_id};
            
            # make sure it hit one of our contigs
            if(exists $global_con_2_len{$con_id}) {
                # make sure there's a container for our read ID
                my $read_id = $sam_fields[0];
                $read_id =~ s/_.$//;
                $read_id =~ s/\.[rf]$//;
                my $position = int($sam_fields[3]);
                my $prob = int($sam_fields[4]);
                if(!exists $global_reads_2_map{$read_id}) {
                    my @tmp = ();
                    $global_reads_2_map{$read_id} = \@tmp;
                }
                
                # add it to the map
                push @{$global_reads_2_map{$read_id}}, $con_id;
                push @{$global_reads_2_map{$read_id}}, $position;
                push @{$global_reads_2_map{$read_id}}, $prob;
                if($sam_fields[1] eq '0') {
                    push @{$global_reads_2_map{$read_id}}, 0;   
                } elsif($sam_fields[1] eq '16') {
                    push @{$global_reads_2_map{$read_id}}, 1; 
                } else {
                    croak "We are expecting sam mapping flags to be either '16' or '0'\n";
                }
            }
        }
    }
    $first_sam = 0;
    close $sam_fh;
}
if(!exists $global_options->{'silent'}) { print "Read ".scalar(keys %global_con_2_len)." contigs from the SAM files.\n"; }




#=================================================================================================
# Remove PCR duplicates. A duplicate is defined as:
# 1. both ends map to the same positions on the contigs
if(!exists $global_options->{'silent'}) {
    print "\n";
    print "Before PCR deduplication, there was ".keys(%global_reads_2_map)." reads (possibly PCR duplicates, single end mapped or chimeric) in \%global_reads_2_map\n"; 
    print "Removing PCR duplicates where both ends are mapped...\n";
}
my %already_mapped_contigs_positions = ();
my $number_of_duplicates_removed = 0;
my $original_number_of_pairs = 0;
foreach my $read_id (keys %global_reads_2_map) {
    my @array = @{$global_reads_2_map{$read_id}};
    if($#array == 7)# if both ends mapped, and no chimeras detected
    {
        $original_number_of_pairs += 1;
        
        my $contig1 = $array[0];
        my $contig2 = $array[4];
        my $position1 = $array[1];
        my $position2 = $array[5];
        my $direction1 = $array[3];
        my $direction2 = $array[7];
        my @key_parts = (
                 $contig1.'_'.$position1.'_'.$direction1,
                 $contig2.'_'.$position2.'_'.$direction2
                 );
        @key_parts = sort @key_parts; #sort so that prdering of the pairs is irrelevant
        my $deduplication_key = $key_parts[0].'__'.$key_parts[1];
        #print $deduplication_key."\n"; #debug
        
        # it is a duplicate if both ends map to the same positions in the same contigs
        if(exists $already_mapped_contigs_positions{$deduplication_key}) {
            $number_of_duplicates_removed += 1;
            print $pcr_duplicate_fh join($delimeter, ($contig1, $position1, $direction1, $contig2, $position2, $direction2))."\n";
            delete $global_reads_2_map{$read_id};
        } else {
            # else this is the first time that it has been mapped
            $already_mapped_contigs_positions{$deduplication_key} = 1;
        }
    }
}
my $number_after_deduplication = $original_number_of_pairs-$number_of_duplicates_removed;
my $percent_duplicate = ($number_of_duplicates_removed/$original_number_of_pairs)*100;
if(!exists $global_options->{'silent'}) {
    print "Removed $number_of_duplicates_removed PCR duplicates, leaving ".$number_after_deduplication." reads pairs with both ends mapped in \%global_reads_2_map.\n";
    print $percent_duplicate."% of sequences were PCR duplicates\n\n";
}
close $pcr_duplicate_fh;





#=================================================================================================
# now go through and make up the results files
if(!exists $global_options->{'silent'}) { print "Determining insert size and orientation...\n"; }
      # stats holders
my @diffs = ();
my @type_array = ();
$type_array[0] = 0;   # <--- --->
$type_array[1] = 0;   # ---> --->
$type_array[2] = 0;   # ---> <---
# conID, pos, prob, strand
foreach my $read_id (keys %global_reads_2_map) {
    my @array = @{$global_reads_2_map{$read_id}};
    if($#array == 7) {
        # paired mapper
        if($array[0] == $array[4]) {
            # mapped onto self
            my $diff = abs($array[1] - $array[5]);
            push @diffs, $diff;
            if($array[3] == $array[7]) {
              $type_array[1]++;
            } else {
                if(($array[1] < $array[5]) ^ ($array[3] == 1)) {
                    $type_array[2]++;
                } else {
                    $type_array[0]++;
                }
            }
        }
    }
#    elsif($#array != 3)
#    {
    # something fishy happening with this hit
    # we only want to keep the two most probably matches
    # ignore for now
#    }
}

# get the type of the read insert
my $true_type = 0;
for my $i (0..2) {
    if($type_array[$i] > $type_array[$true_type]) {
       $true_type = $i;
    }
}

# calculate the stdev on only the middle 80% of the sorted differences, otherwise it can be become very inflated. Can it still do this?
@diffs = sort @diffs;
my $length = scalar @diffs;
my $bottom = int($length*0.1);
my $top = int($length*0.9);
croak "not enough data mapped to within contigs to determine insert size empirically" if $top-1 <= $bottom;

my $cum_diff = 0.0;
my $num_counted = 0;
for my $index ($bottom..$top){
    $cum_diff += $diffs[$index];
    $num_counted += 1;
}

my $mean = $cum_diff / $num_counted;
my $stdev = 0.0;
for my $index ($bottom..$top)
{
    $stdev += ($diffs[$index] - $mean)**2;
}
$stdev /= ($num_counted-1);
$stdev = $stdev ** 0.5;

# we only care about reads which match the given type and
# have an insert size comparable to the distribution
# calculated above. ALSO...
# use look up tables to save on if statements
my $tol = 2;
my $lower_limit = 0; #no need for this 
my $upper_limit = overrideDefault(($mean + $stdev * $tol), 'max_insert');
if(!exists $global_options->{'silent'}) {
    print "Stats estimate:\n";
    print "Mean: $mean\nStdev: $stdev\n";
    print "Not accepting insert greater than $upper_limit or less than $lower_limit.\n";
    print "Found $type_array[2] read pairs facing inwards on the same contig (type 2 read pairs). This is what you want for IonTorrent mate pair data.\n";
    print "Found $type_array[0] read pairs facing outwards on the same contig (type 0 read pairs). This is not what you want for mate pair data.\n";
    print "Found $type_array[1] read pairs facing the same direction on the same contig (type 1 read pairs). This is not what you want for IonTorrent mate pair data.\n";
    print "\nProceding using type $true_type as the expected mate pair type\n*****\n";
}
##
# KEYS:
#
# ([1 at end],[1 rev],[2 at end],[2 rev])
# 0 False, 1 True
##

##
# READ $types:
#
##
my %type_table = ();
# S: start
# E: end
# F: forward
# R: reverse
$type_table{'0000'} = 0; # 1SF 2SF
$type_table{'0001'} = 1; # 1SF 2SR
$type_table{'0010'} = 1; # 1SF 2EF
$type_table{'0011'} = 0; # 1SF 2ER
$type_table{'0100'} = 1; # 1SR 2SF
$type_table{'0101'} = 2; # 1SR 2SR
$type_table{'0110'} = 2; # 1SR 2EF
$type_table{'0111'} = 1; # 1SR 2ER
$type_table{'1000'} = 1; # 1EF 2SF
$type_table{'1001'} = 2; # 1EF 2SR
$type_table{'1010'} = 2; # 1EF 2EF
$type_table{'1011'} = 1; # 1EF 2ER
$type_table{'1100'} = 0; # 1ER 2SF
$type_table{'1101'} = 1; # 1ER 2SR
$type_table{'1110'} = 1; # 1ER 2EF
$type_table{'1111'} = 0; # 1ER 2ER

if(!exists $global_options->{'silent'}) { print "Munging it all together...\n"; }
foreach my $read_id (keys %global_reads_2_map) {
    # conID, pos, prob, strand
    my @array = @{$global_reads_2_map{$read_id}};
    my $contig1_name = $global_int_2_con{$array[0]};
    my $contig1_length = $global_con_2_len{$array[0]}; 
    if($#array == 3) {
        # single mapper
        print $unpaired_fh join($delimeter, ($contig1_name,$contig1_length,$read_id,$array[1],$array[3]))."\n";
    } elsif($#array == 7) {
        # paired mapper
        # 0      1    2     3       4      5    6     7
        # conID, pos, prob, strand, conID, pos, prob, strand
        my $contig2_name = $global_int_2_con{$array[4]};
        my $contig2_length = $global_con_2_len{$array[4]}; 
        
        if($contig1_name eq $contig2_name) {
            # mapped onto self
            my $pair_print = join($delimeter, ($contig1_name,$contig1_length,$read_id,$array[1],$array[3],$array[5],$array[7]));
            if(($array[1] < $array[5]) ^ ($array[3] == 1) ^ ($true_type == 0)) {
                # seems OK, check the insert
                my $diff = abs($array[1] - $array[5]);
                if(($diff > $lower_limit) and ($diff < $upper_limit)) {
                    # all good
                    print $paired_fh $pair_print."\n";
                } else {
                    # insert mismatch
                    print $error_paired_fh $pair_print.$delimeter."I\n";
                }
            } else {
                # type mismatch
                print $error_paired_fh $pair_print.$delimeter."T\n";
            }
        } else {
            # possible linker
            
            # make the print statements consistant
            my $all_print = join $delimeter, ($contig2_name, $array[5], $array[6], $array[7],$contig1_name, $array[1], $array[2], $array[3]);
            my $link_print = join($delimeter, ($contig2_name,$contig2_length,$contig1_name,$contig1_length,$read_id,$array[5],$array[7],$array[1],$array[3]));
            if($contig1_name lt $contig2_name) {
                $all_print = join $delimeter, ($contig1_name, $array[1], $array[2], $array[3],$contig2_name, $array[5], $array[6], $array[7]);
                $link_print = join($delimeter, ($contig1_name,$contig1_length,$contig2_name,$contig2_length,$read_id,$array[1],$array[3],$array[5],$array[7]));
            }

            # print to the all links file for making the dot file
            print $all_links_fh $all_print."\n";
            
            # work out whether the reads map in the right way
            if(($contig1_length < $mean) or ($contig2_length < $mean)) {
                # one contig or another is too short.
                print $short_links_fh $link_print."\n";
            } 
        }
    }
}

close $paired_fh;
close $unpaired_fh;
close $error_paired_fh;
close $short_links_fh;
close $all_links_fh;

# Call the ruby script that does the creation of the graphviz file

my $path_to_contig_linker = File::Spec->catfile(dirname($0), 'contig_linker.rb');
checkAndRunCommand($path_to_contig_linker, [{
                                 -l => $all_links_file,
                                 -L => $global_minimum_links,
                                 -m => $upper_limit,
                                 -f => $reference_fasta,
                                 -u => $mean,
                                 "--trace" => 'info',
                                 -u => $upper_limit
                                 }], DIE_ON_FAILURE);

######################################################################
# CUSTOM SUBS
######################################################################
sub dec2bin { my ($dec) = @_; return sprintf "%08b", $dec; }

sub revcompl {
    my ($seq) = @_;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return scalar reverse $seq;
}


######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ("sam1|1:s", "sam2|2:s", "reference_fasta|f:s", "working_dir|w:s", "help|h+", "max_insert|m:i", "min_links|l:i", "silent|s+");
    my %options;

    # Add any other command line options, and the code to handle them
    #
    GetOptions( \%options, @standard_options);

    # Default to current directory if no directory is given
    $options{'working_dir'} = '.' unless defined($options{'working_dir'});

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    if(!(exists $options{'sam1'})) { printParamError ("No forward SAM file supplied."); }
    if(!(exists $options{'sam2'})) { printParamError ("No reverse SAM file supplied."); }
    if(!(exists $options{'reference_fasta'})) { printParamError ("No reference fasta file supplied."); }

    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #
    my ($error) = @_;
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) {
        return $global_options->{$option_name};
    }
    return $default_value;
}

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
        unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    #
    my ($cmd, $params, $failure_type) = @_;
    
    isCommandInPath($cmd, $failure_type);
    
    # join the parameters to the command
    my $param_str = join " ", map {formatParams($_)} @{$params};
    
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    # make sure that all went well
    if (system($cmd_str)) {
        handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to
    # checkAndRunCommand
    #
    my $ref = shift;
    
    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
    'only contain references to arrays or hashes\n';
}


sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : " . $! . "\n";
        }
    }
}


######################################################################
# MISC

sub printAtStart {
print<<"EOF";
----------------------------------------------------------------
$0
Copyright (C) Michael Imelfort, Ben Woodcroft
This program comes with ABSOLUTELY NO WARRANTY;
This is free software, and you are welcome to redistribute it
under certain conditions: See the source for more details.
----------------------------------------------------------------
EOF
}

__DATA__

=head1 NAME

__Script__Name__

=head1 COPYRIGHT

    copyright (C) Michael Imelfort, Ben Woodcroft
    
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

    Parse a pair of sam files and determine read pairs which link contigs 

=head1 SYNOPSIS

    ionPairer.pl -sam1|1 SAMFILE1 -sam2|2 SAMFILE2 -reference_fasta|f REFERENCE_FASTA

    -sam1 -1 SAMFILE1                    Sam file of forward read
    -sam2 -2 SAMFILE2                    Sam file of reverse read
    -reference_fasta -f REFERENCE_FASTA  Fasta file of sequences to be scaffolded (used during creation of sam files)
    [-working_dir -w]                    Somewhere to write all the files
    [-min_links -l]                      Minimum number of links to report in the dot file [default: 3]
    [-max_insert -m]                     Maximum insert size accepted [default: mean + 2 * std measured empirically (middle 80% paired reads)]
    [-silent -s]                         Suppress print statements

    Produces output files:
    
     SAMFILE1.unpaired.csv                           - Reads where only one end mapped to a contig
     SAMFILE1.pcr_duplicates.csv                     - Reads removed from further analysis as they were judged to be PCR duplicates
     SAMFILE1.paired.csv                             - Reads where both ends mapped onto one contig
     SAMFILE1.error_paired.csv                       - Reads where both ends mapped, but erroneously due to insert size or relative orientation
     SAMFILE1.unique_links.csv                       - All pairs of mate pairs that span between two contigsand pass the pcr filter  
     SAMFILE1.short_links.csv                        - Pairs where one of the contigs was shorter than the insert size
     SAMFILE1.unique_links.csv.error_links.csv       - Pairs which link two contigs, but erroneously due to insert size, position or relative orientation
     SAMFILE1.unique_links.csv.filtered_links.csv    - Subset of unique pairs which are neither short nor erroneous
     
=cut
