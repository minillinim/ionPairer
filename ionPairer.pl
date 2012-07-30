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
printAtStart();
my $global_options = checkParams();

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
my %global_con_2_len = ();         # contig integerss versus lengths
my %global_reads_2_map = ();       # reads vs mapping info

# get output file names and handles
my $l_file = File::Spec->catfile($global_options->{'working_dir'}, $global_options->{'sam1'}.".links.csv");
my $p_file = File::Spec->catfile($global_options->{'working_dir'}, $global_options->{'sam1'}.".paired.csv");
my $u_file = File::Spec->catfile($global_options->{'working_dir'}, $global_options->{'sam1'}.".unpaired.csv");
my $ep_file = File::Spec->catfile($global_options->{'working_dir'}, $global_options->{'sam1'}.".perrors.csv");
my $el_file = File::Spec->catfile($global_options->{'working_dir'}, $global_options->{'sam1'}.".lerrors.csv");
my $al_file = File::Spec->catfile($global_options->{'working_dir'}, $global_options->{'sam1'}.".lamb.csv");
my $s_file = File::Spec->catfile($global_options->{'working_dir'}, $global_options->{'sam1'}.".summary.gv");
my $pcr_duplicate_file = File::Spec->catfile($global_options->{'working_dir'}, $global_options->{'sam1'}.".pcr_duplicates.csv");

open my $l_fh, ">", $l_file or die "**ERROR: could not open link output file $l_file $!\n";
open my $p_fh, ">", $p_file or die "**ERROR: could not open paired output file $p_file $!\n";
open my $u_fh, ">", $u_file or die "**ERROR: could not open unpaired output file $u_file $!\n";
open my $ep_fh, ">", $ep_file or die "**ERROR: could not open error_pair output file $ep_file $!\n";
open my $el_fh, ">", $el_file or die "**ERROR: could not open error_link output file $el_file $!\n";
open my $al_fh, ">", $al_file or die "**ERROR: could not open ambiguous output file $al_file $!\n";
open my $s_fh, ">", $s_file or die "**ERROR: could not open summary output file $s_file $!\n";
open my $pcr_duplicate_fh, ">", $pcr_duplicate_file or die "**ERROR: could not open pcr duplicate output file $s_file $!\n";

# parse the sam files
my @samfiles = ($global_options->{'sam1'},$global_options->{'sam2'});
my $first_sam = 1;
foreach my $sam_fn (@samfiles)
{
    print "Parsing: $sam_fn\n";
    open my $sam_fh, "<", $sam_fn or die "**ERROR: could not open sam file $!\n";
    while(<$sam_fh>)
    {
        chomp $_;
        if($_ =~ /^@/)
        {
            # work out and store contig lengths
            next if(1 != $first_sam);
            my @fields = split(/:/, $_);
            $fields[1] =~ s/\tLN//;
            my $con_id = $fields[1];
            my $con_int = 0;
            if(!exists $global_con_2_int{$con_id})
            {
                $con_int = $global_conInt;
                $global_conInt++;
                $global_con_2_int{$con_id} = $con_int;
                $global_int_2_con{$con_int} = $con_id;
            }
            else
            {
                $con_int = $global_con_2_int{$con_id};
            }
            $global_con_2_len{$con_int} = $fields[2];         
        }
        else
        {
            my @sam_fields = split(/\t/, $_);
            my @mapping_flags = split(//, dec2bin($sam_fields[1]));
            
            # make sure it mapped
            next if($mapping_flags[5] eq "1");
            
            # get a contig ID (int)
            my $con_id = $sam_fields[2];
            if(!exists $global_con_2_int{$con_id})
            {
                croak "Unknown contig ID: $con_id\n";
            }
            $con_id = $global_con_2_int{$con_id};
            
            # make sure it hit one of our contigs
            if(exists $global_con_2_len{$con_id})
            {
                # make sure there's a container for our read ID
                my $read_id = $sam_fields[0];
                $read_id =~ s/_.$//;
                $read_id =~ s/\.[rf]$//;
                my $position = int($sam_fields[3]);
                my $prob = int($sam_fields[4]);
                if(!exists $global_reads_2_map{$read_id})
                {
                    my @tmp = ();
                    $global_reads_2_map{$read_id} = \@tmp;
                }
                
                # add it to the map
                push @{$global_reads_2_map{$read_id}}, $con_id;
                push @{$global_reads_2_map{$read_id}}, $position;
                push @{$global_reads_2_map{$read_id}}, $prob;
                if($sam_fields[1] eq '0')
                {
                    push @{$global_reads_2_map{$read_id}}, 0;   
                }
                elsif($sam_fields[1] eq '16')
                {
                    push @{$global_reads_2_map{$read_id}}, 1; 
                }
                else
                {
                    croak "We are expecting sam mapping flags to be either '16' or '0'\n";
                }
            }
        }
    }
    $first_sam = 0;
    close $sam_fh;
}


#=================================================================================================
# Remove PCR duplicates. A duplicate is defined as:
# 1. both ends map to the same positions on the contigs
print "\n";
print "Before PCR deduplication, there was ".keys(%global_reads_2_map)." reads (possibly PCR duplicates, single end mapped or chimeric) in \%global_reads_2_map\n";
print "Removing PCR duplicates where both ends are mapped...\n";
my %already_mapped_contigs_positions = ();
my $number_of_duplicates_removed = 0;
my $original_number_of_seqs = keys(%global_reads_2_map);
foreach my $read_id (keys %global_reads_2_map)
{
    my @array = @{$global_reads_2_map{$read_id}};
    if($#array == 7)# if both ends mapped, and no chimeras detected
    {
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
        if(exists $already_mapped_contigs_positions{$deduplication_key}){
            $number_of_duplicates_removed += 1;
            print $pcr_duplicate_fh join("\t", ($contig1, $position1, $direction1, $contig2, $position2, $direction2))."\n";
            delete $global_reads_2_map{$read_id};
        } else {
            # else this is the first time that it has been mapped
            $already_mapped_contigs_positions{$deduplication_key} = 1;
        }
    }
}
print "Removed $number_of_duplicates_removed PCR duplicates, leaving ".keys(%global_reads_2_map)." reads in \%global_reads_2_map.\n";
my $percent_duplicate = (1.0-keys(%global_reads_2_map)/$original_number_of_seqs)*100;
print $percent_duplicate."% of sequences were PCR duplicates\n\n";





#=================================================================================================
# now go through and make up the results files
print "Determining insert size and orientation...\n";
my $cum_diff = 0;       # stats holders
my @diffs = ();
my @type_array = ();
$type_array[0] = 0;   # <--- --->
$type_array[1] = 0;   # ---> --->
$type_array[2] = 0;   # ---> <---
# conID, pos, prob, strand
foreach my $read_id (keys %global_reads_2_map)
{
    my @array = @{$global_reads_2_map{$read_id}};
    if($#array == 7)
    {
        # paired mapper
        if($array[0] == $array[4])
        {
            # mapped onto self
            my $diff = abs($array[1] - $array[5]);
            push @diffs, $diff;
            $cum_diff += $diff;
            if($array[3] == $array[7])
            {
                $type_array[1]++;
            }
            else
            {
                if(($array[1] < $array[5]) ^ ($array[3] == 1))
                {
                    $type_array[2]++;
                }
                else
                {
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
for my $i (0..2)
{
    if($type_array[$i] > $type_array[$true_type])
    {
        $true_type = $i;
    }
}

my $mean = $cum_diff / (scalar @diffs);
my $stdev = 0;
for my $val (@diffs)
{
    $stdev += ($val - $mean)**2;
}
$stdev /= (scalar @diffs-1);
$stdev = $stdev ** 0.5;

print "Stats estimate:\n";
print "Mean: $mean\nStdev: $stdev\nRAW: ";
print "Found $type_array[2] read pairs facing inwards on the same contig (type 2). This is what you want for IonTorrent mate pair data.\n";
print "Found $type_array[0] read pairs facing outwards on the same contig (type 0). This is not what you want for mate pair data.\n";
print "Found $type_array[1] read pairs facing the same direction on the same contig (type 1). This is not what you want for IonTorrent mate pair data.\n";
print "\nProceding using type $true_type as the expected mate pair type\n*****\n";

# we only care about reads which match the given type and
# have an insert size comparable to the distribution
# calculated above. ALSO...
# use look up tables to save on if statements
my $tol = 2;
my $lower_limit = $mean - $stdev * $tol; 
my $upper_limit = $mean + $stdev * $tol; 

##
# KEYS:
# 
# ([1 at end],[1 rev],[2 at end],[2 rev])
# 0 False, 1 True
##

##
# READ $typeS:
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

##
# $orientATIONS:
#
#    2   1   2
#
#       --> -->     2F1A    0
#       --> <--     2F1D    1
#   --> -->         1F2A    2
#   <-- -->         1F2D    3
#
##
my %orient_table = ();
$orient_table{'0000'} = 3;
$orient_table{'0001'} = 3;
$orient_table{'0010'} = 2;
$orient_table{'0011'} = 2;
$orient_table{'0100'} = 3;
$orient_table{'0101'} = 3;
$orient_table{'0110'} = 2;
$orient_table{'0111'} = 2;
$orient_table{'1000'} = 0;
$orient_table{'1001'} = 0;
$orient_table{'1010'} = 1;
$orient_table{'1011'} = 1;
$orient_table{'1100'} = 0;
$orient_table{'1101'} = 0;
$orient_table{'1110'} = 1;
$orient_table{'1111'} = 1;
my %print_table = ();
$print_table{'0000'} = "1F2D";
$print_table{'0001'} = "1F2D";
$print_table{'0010'} = "1F2A";
$print_table{'0011'} = "1F2A";
$print_table{'0100'} = "1F2D";
$print_table{'0101'} = "1F2D";
$print_table{'0110'} = "1F2A";
$print_table{'0111'} = "1F2A";
$print_table{'1000'} = "2F1A";
$print_table{'1001'} = "2F1A";
$print_table{'1010'} = "2F1D";
$print_table{'1011'} = "2F1D";
$print_table{'1100'} = "2F1A";
$print_table{'1101'} = "2F1A";
$print_table{'1110'} = "2F1D";
$print_table{'1111'} = "2F1D";

# hash of hashes for dertermining insert / orientation
my %post = ();
my %pre = ();
my %post_dist = ();
my %pre_dist = ();

# smallest to largest con_id int
# whoever holds the most is telling the truth
my @orients = ({},{},{},{});

print "Munging it all together...\n";
foreach my $read_id (keys %global_reads_2_map)
{
    # conID, pos, prob, strand
    my @array = @{$global_reads_2_map{$read_id}};
    if($#array == 3)
    {
        # single mapper
        my $con_id = $global_int_2_con{$array[0]};
        print $u_fh join("\t", ($con_id,$global_con_2_len{$array[0]},$read_id,$array[1],$array[3]))."\n";
    }
    elsif($#array == 7)
    {
        # paired mapper
        if($array[0] == $array[4])
        {
            # mapped onto self
            my $con_id = $global_int_2_con{$array[0]};
            if(($array[1] < $array[5]) ^ ($array[3] == 1) ^ ($true_type == 0))
            {
                # seems OK, check the insert
                my $diff = abs($array[1] - $array[5]);
                if(($diff > $lower_limit) and ($diff < $upper_limit))
                {
                    # all good
                    print $p_fh join("\t", ($con_id,$global_con_2_len{$array[0]},$read_id,$array[1],$array[3],$array[5],$array[7]))."\n";
                }
                else
                {
                    # insert mismiatch
                    print $ep_fh join("\t", ($con_id,$global_con_2_len{$array[0]},$read_id,$array[1],$array[3],$array[5],$array[7],"I"))."\n";
                }
            }
            else
            {
                # type mismatch
                print $ep_fh join("\t", ($con_id,$global_con_2_len{$array[0]},$read_id,$array[1],$array[3],$array[5],$array[7],"T"))."\n";
            }
        }
        else
        {
            # 0      1    2     3       4      5    6     7
            # conID, pos, prob, strand, conID, pos, prob, strand
            # possible linker
            # work out whether the reads map in the right way
            my $con_1_len = $global_con_2_len{$array[0]}; 
            my $con_2_len = $global_con_2_len{$array[4]}; 
            if(($con_1_len < 2*$upper_limit) or ($con_2_len < 2*$upper_limit))
            {
                # one contig or another is too short.
                if ($con_1_len < $con_2_len)
                {
                    print $al_fh join("\t", ($global_int_2_con{$array[0]},$con_1_len,$global_int_2_con{$array[4]},$con_2_len,$read_id,$array[1],$array[3],$array[5],$array[7]))."\n";
                }
                else
                {
                    print $al_fh join("\t", ($global_int_2_con{$array[4]},$con_2_len,$global_int_2_con{$array[0]},$con_1_len,$read_id,$array[5],$array[7],$array[1],$array[3]))."\n";
                }
            }
            else
            {
                # keep the distances
                my $dist_1 = 0;
                my $dist_2 = 0;
                
                # make a key!
                my @key = ('0','0','0','0');
                
                # strandedness
                if(1 == $array[3]) { $key[1] = '1'; }
                if(1 == $array[7]) { $key[3] = '1'; }
                if($array[1] < $upper_limit) { $key[0] = '0'; $dist_1 = $array[1]; }
                elsif($array[1] > ($con_1_len - $upper_limit)) {  $key[0] = '1'; $dist_1 = $con_1_len - $array[1]; }
                else
                {
                    # read one lies too far into the contig
                    print $el_fh join("\t", ($global_int_2_con{$array[0]},$con_1_len,$global_int_2_con{$array[4]},$con_2_len,$read_id,$array[1],$array[3],$array[5],$array[7],"1"))."\n";
                    next;
                }
                
                if($array[5] < $upper_limit) { $key[2] = '0'; $dist_2 = $array[5]; }
                elsif($array[5] > ($con_2_len - $upper_limit)) { $key[2] = '1'; $dist_2 = $con_2_len - $array[5];}
                else
                {
                    # read two lies too far into the contig
                    print $el_fh join("\t", ($global_int_2_con{$array[0]},$con_1_len,$global_int_2_con{$array[4]},$con_2_len,$read_id,$array[1],$array[3],$array[5],$array[7],"2"))."\n";
                    next;
                }
                my $index_key = join("", @key);
                my $k_type = $type_table{$index_key};
                if($k_type == $true_type)
                {
                    print $l_fh join("\t", ($global_int_2_con{$array[0]},$con_1_len,$global_int_2_con{$array[4]},$con_2_len,$read_id,$array[1],$array[3],$array[5],$array[7]))."\n";
                    my $orient = $orient_table{$index_key};
                    my $insert = int($mean - $dist_1 - $dist_2);
                    if(0 == $orient)           # ---1--> ---2-->
                    {
                        incOrients(\@orients, 0, $array[0],$array[4]);
                        if(!exists $pre{$array[4]}) { my %tmp = (); $pre{$array[4]} = \%tmp; }
                        if(!exists $post{$array[0]}) { my %tmp = (); $post{$array[0]} = \%tmp; }
                        ${$pre{$array[4]}}{$array[0]}++;
                        ${$post{$array[0]}}{$array[4]}++;

                        if(!exists $pre_dist{$array[4]}) { my %tmp = (); $pre_dist{$array[4]} = \%tmp; }
                        if(!exists $post_dist{$array[0]}) { my %tmp = (); $post_dist{$array[0]} = \%tmp; }
                        if(!exists ${$pre_dist{$array[4]}}{$array[0]}) { my @tmp = (); ${$pre_dist{$array[4]}}{$array[0]} = \@tmp; }
                        if(!exists ${$post_dist{$array[0]}}{$array[4]}) { my @tmp = (); ${$post_dist{$array[0]}}{$array[4]} = \@tmp; }
                        push @{${$pre_dist{$array[4]}}{$array[0]}}, $insert;
                        push @{${$post_dist{$array[0]}}{$array[4]}}, $insert;

                        #print "$global_int_2_con{$array[4]} follows $global_int_2_con{$array[0]} agrees $insert\n";
                    }
                    elsif(1 == $orient)        # ---1--> <--2---
                    {
                        incOrients(\@orients, 1, $array[0],$array[4]);
                        if(!exists $post{$array[0]}) { my %tmp = (); $post{$array[0]} = \%tmp; }
                        if(!exists $post{$array[4]}) { my %tmp = (); $post{$array[4]} = \%tmp; }
                        ${$post{$array[0]}}{$array[4]}++;
                        ${$post{$array[4]}}{$array[0]}++;

                        if(!exists $post_dist{$array[0]}) { my %tmp = (); $post_dist{$array[0]} = \%tmp; }
                        if(!exists $post_dist{$array[4]}) { my %tmp = (); $post_dist{$array[4]} = \%tmp; }
                        if(!exists ${$post_dist{$array[0]}}{$array[4]}) { my @tmp = (); ${$post_dist{$array[0]}}{$array[4]} = \@tmp; }
                        if(!exists ${$post_dist{$array[4]}}{$array[0]}) { my @tmp = (); ${$post_dist{$array[4]}}{$array[0]} = \@tmp; }
                        push @{${$post_dist{$array[0]}}{$array[4]}}, $insert;
                        push @{${$post_dist{$array[4]}}{$array[0]}}, $insert;

                        #print "$global_int_2_con{$array[4]} follows $global_int_2_con{$array[0]} disagrees $insert\n";
                    }
                    elsif(2 == $orient)        # ---2--> ---1-->
                    {
                        incOrients(\@orients, 2, $array[0],$array[4]);
                        if(!exists $pre{$array[0]}) { my %tmp = (); $pre{$array[0]} = \%tmp; }
                        if(!exists $post{$array[4]}) { my %tmp = (); $post{$array[4]} = \%tmp; }
                        ${$pre{$array[0]}}{$array[4]}++;
                        ${$post{$array[4]}}{$array[0]}++;

                        if(!exists $pre_dist{$array[0]}) { my %tmp = (); $pre_dist{$array[0]} = \%tmp; }
                        if(!exists $post_dist{$array[4]}) { my %tmp = (); $post_dist{$array[4]} = \%tmp; }
                        if(!exists ${$pre_dist{$array[0]}}{$array[4]}) { my @tmp = (); ${$pre_dist{$array[0]}}{$array[4]} = \@tmp; }
                        if(!exists ${$post_dist{$array[4]}}{$array[0]}) { my @tmp = (); ${$post_dist{$array[4]}}{$array[0]} = \@tmp; }
                        push @{${$pre_dist{$array[0]}}{$array[4]}}, $insert;
                        push @{${$post_dist{$array[4]}}{$array[0]}}, $insert;

                        #print "$global_int_2_con{$array[0]} follows $global_int_2_con{$array[4]} agrees $insert\n";
                    }
                    else #(3 == $orient)       # <--2--- ---1-->
                    {
                        incOrients(\@orients, 3, $array[0],$array[4]);
                        if(!exists $pre{$array[4]}) { my %tmp = (); $pre{$array[4]} = \%tmp; }
                        if(!exists $pre{$array[0]}) { my %tmp = (); $pre{$array[0]} = \%tmp; }
                        ${$pre{$array[4]}}{$array[0]}++;
                        ${$pre{$array[0]}}{$array[4]}++;

                        if(!exists $pre_dist{$array[4]}) { my %tmp = (); $pre_dist{$array[4]} = \%tmp; }
                        if(!exists $pre_dist{$array[0]}) { my %tmp = (); $pre_dist{$array[0]} = \%tmp; }
                        if(!exists ${$pre_dist{$array[4]}}{$array[0]}) { my @tmp = (); ${$pre_dist{$array[4]}}{$array[0]} = \@tmp; }
                        if(!exists ${$pre_dist{$array[0]}}{$array[4]}) { my @tmp = (); ${$pre_dist{$array[0]}}{$array[4]} = \@tmp; }
                        push @{${$pre_dist{$array[4]}}{$array[0]}}, $insert;
                        push @{${$pre_dist{$array[0]}}{$array[4]}}, $insert;

                        #print "$global_int_2_con{$array[0]} follows $global_int_2_con{$array[4]} disagrees $insert\n";
                    }
                }
                else
                {
                    print $el_fh join("\t", ($global_int_2_con{$array[0]},$con_1_len,$global_int_2_con{$array[4]},$con_2_len,$read_id,$array[1],$array[3],$array[5],$array[7],"T"))."\n";
                }
            } 
        }
    }
}

close $l_fh;
close $p_fh;
close $u_fh;
close $el_fh;
close $ep_fh;
close $al_fh;

my $min_hits = 3;
my %global_GV_nodes = ();
my @global_GV_edges = ();
my %global_GV_seen_edges = ();

print "Making weighted graph\n";
foreach my $con_id (keys %global_con_2_len)
{
    if(exists $pre{$con_id})
    {
        # something lies before this guy
        # make a node for this contig... ...perhaps
        makeGVNode($con_id);
        # find the guys above the lower cutoff
        my %valid_links = ();
        my %tmp_hash = %{$pre{$con_id}};
        foreach my $link (keys %tmp_hash)
        {
            if($tmp_hash{$link} >= $min_hits)
            {
                $valid_links{$link} = $tmp_hash{$link};
            }
        }
        
        # now make the edges
        %tmp_hash = %{$pre_dist{$con_id}};
        foreach my $link (keys %valid_links)
        {
            my $insert = 0;
            my $num_ins = 0;
            foreach my $ins (@{$tmp_hash{$link}})
            {
                $insert += $ins;
                $num_ins++; 
            }
            $insert = int($insert/$num_ins); 
            makeGVedge(\@orients, $con_id, $link, $insert,  $num_ins);
            #print $global_int_2_con{$link} . "\t" . $global_con_2_len{$link}. "\t" . $insert . "\t" . $num_ins . "\n";
        }
    }
    
    if(exists $post{$con_id})
    {
        # something lies after this guy
        # make a node for this contig... ...perhaps
        makeGVNode($con_id);
        my %valid_links = ();
        my %tmp_hash = %{$post{$con_id}};
        foreach my $link (keys %tmp_hash)
        {
            if($tmp_hash{$link} >= $min_hits)
            {
                $valid_links{$link} = $tmp_hash{$link};
            }
        }
        %tmp_hash = %{$post_dist{$con_id}};
        foreach my $link (keys %valid_links)
        {
            my $insert = 0;
            my $num_ins = 0;
            foreach my $ins (@{$tmp_hash{$link}})
            {
                $insert += $ins;
                $num_ins++; 
            }
            $insert = int($insert/$num_ins); 
            makeGVedge(\@orients, $con_id, $link, $insert,  $num_ins);
        }
    }
}

# now print the graph
print $s_fh "digraph auto {\tnode [shape = point];\n";
foreach my $key (keys %global_GV_nodes)
{
    print $s_fh $global_GV_nodes{$key};
}
foreach my $edge (@global_GV_edges)
{
    print $s_fh $edge;
}
print $s_fh "};\n";
close $s_fh;

######################################################################
# CUSTOM SUBS
######################################################################
sub incOrients
{
    #-----
    # increment counters for contig orientation
    #
    my ($oref, $orient, $id_1, $id_2) = @_;
    my $key = 0;
    if($id_1 < $id_2) 
    {
        $key = $id_2 * 10000000 + $id_1; 
        if(0 == $orient) { $orient = 2; }
        elsif(2 == $orient) { $orient = 0; }
    }
    else { $key = $id_1 * 10000000 + $id_2; }
    ${$oref}[$orient]{$key}++;
}

sub getOrient
{
    #-----
    # increment counters for contig orientation
    #
    my ($oref, $id_1, $id_2) = @_;
    
    # get the key
    my $switch = 0;
    my $key = 0;
    if($id_1 < $id_2) { $key = $id_2 * 10000000 + $id_1; $switch = 1; }
    else { $key = $id_1 * 10000000 + $id_2; }
    
    # we only want to return one copy of each edge
    if(exists $global_GV_seen_edges{$key})
    {
        return -1;
    }
    $global_GV_seen_edges{$key} = 1;
    
    # find the maximum one
    my $max_orient_count = 0;
    my $orient = 0;
    for my $i (0..3)
    {
        if(exists ${$oref}[$i]{$key})
        {
            if(${$oref}[$i]{$key} > $max_orient_count) 
            {
                $max_orient_count = ${$oref}[$i]{$key}; 
                $orient = $i; 
            }
        }
    }
    
    # fix any switching
    if(1 == $switch)
    {
        if(0 == $orient) { $orient = 2; }
        elsif(2 == $orient) { $orient = 0; }
    }
    return $orient;
}

sub makeGVNode
{
    my ($con_id)  = @_;
    if(!exists $global_GV_nodes{$con_id})
    {
        #
        # JUST for this CLC data
        #
        my $pretty_name = $global_int_2_con{$con_id};
        $pretty_name =~ s/pairedrawreadstrimmed\(paired\)//;
        my $node_str  = "\t\"".$global_int_2_con{$con_id}."START\" -> \"".$global_int_2_con{$con_id}."END\" [ label = \"$pretty_name"."_$global_con_2_len{$con_id}\" style=\"setlinewidth(7)\"];\n";
        $global_GV_nodes{$con_id} = $node_str;
    }
}

sub makeGVedge
{
    #-----
    # Make edges in the graph
    #
    my ($oref, $con_id, $link, $insert,  $num_ins) = @_;
    my $orient = getOrient($oref, $con_id, $link);
    # return of -1 implies the edge is alreay seen
    if($orient == 0)
    {
        my $edge_str = "\t\"".$global_int_2_con{$con_id}."END\" -> \"".$global_int_2_con{$link}."START\" [ label = \"".$insert."_".$num_ins."\"];\n"; 
        push @global_GV_edges, $edge_str;
    }
    elsif($orient == 1)
    {
        my $edge_str = "\t\"".$global_int_2_con{$con_id}."END\" -> \"".$global_int_2_con{$link}."END\" [ label = \"".$insert."_".$num_ins."\"];\n";
        push @global_GV_edges, $edge_str;
    }
    elsif($orient == 2)
    {
        my $edge_str = "\t\"".$global_int_2_con{$con_id}."START\" -> \"".$global_int_2_con{$link}."END\" [ label = \"".$insert."_".$num_ins."\"];\n";
        push @global_GV_edges, $edge_str;
    }
    elsif($orient == 3)
    {
        my $edge_str = "\t\"".$global_int_2_con{$con_id}."START\" -> \"".$global_int_2_con{$link}."START\" [ label = \"".$insert."_".$num_ins."\"];\n";
        push @global_GV_edges, $edge_str;
    }
}

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
    my @standard_options = ("sam1|1:s", "sam2|2:s", "working_dir|w:s", "help|h+");
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
    if(exists $global_options->{$option_name})
    {
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

    ionPairer.pl -sam1|1 SAMFILE1 -sam2|2 SAMFILE2

        -sam1 -1 SAMFILE1               Sam file of forward read
        -sam2 -2 SAMFILE2               Sam file of reverse read
        [-working_dir -w]               Somewhere to write all the files

    Produces output files:
    
         SAMFILE1.links.csv       - Reads which link two contigs
         SAMFILE1.paired.csv      - Reads where both ends mapped onto one contig
         SAMFILE1.lerror.csv      - Reads which link two contigs, but erroneously
         SAMFILE1.perror.csv      - Reads where both ends mapped, but erroneously
         SAMFILE1.lamb.csv        - Ambiguous links
         SAMFILE1.unpaired.csv    - Reads where only one end mapped
         SAMFILE1.summary.gv      - Graphviz file for visualising the links
=cut
