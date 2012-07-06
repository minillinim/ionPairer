#!/usr/bin/env perl
###############################################################################
##
##    scaffolder.pl
##    
##    Given a graphviz file and a contigs file scaffold them together
##
##    Copyright (C) Fauzi Haroon, Connor Skennerton
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
################################################################################
use strict;
use warnings;
#core Perl modules
use Getopt::Long;
use Data::Dumper;
use Carp;

#CPAN modules
use Bio::SeqIO;

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
# globals
my %global_start_hash = (); #contig pairs that begin with START
my %global_end_hash = (); #contig pairs that begin with END

# First read all the sequences into a hash
# use contigID -> seq
my %global_seq_hash = ();
my %global_seq_printed_hash = (); # in this hash means printed to file
my $global_contig_id = 1; #increment the number when printing, so that it is contig1 ... contig2 ...

# need to read file
# populate start and end hashes
open my $gv_fh, "<",$global_options->{'gv'}  || die $! ;
while (<$gv_fh>) {
    chomp $_;
    next if ($_ =~ /digraph/);
    next if ($_ =~ /style/);
	next if ($_ =~ /};/);
	$_ =~ s/"//g;
	$_ =~ s/\[.*\];//g;
	$_ =~ s/^\s+//;
	# Parse both pairs in .gv
	my @fields = split(/ -> /, $_);
	# Store the name of the contigs with the START/END
    my ($pair1) = $fields[0];
    my ($pair2) = $fields[1];
	# trim the name to remove the START and then store in the global start hash
    my ($pair1_trim) = $pair1 =~ /(.+)(START|END)/;
    my ($pair2_trim) = $pair2 =~ /(.+)(START|END)/;
	# if the name of the contigs on pair1 contains START
	if($pair1 =~ /START/) {
		# to check if the contig already exists in the hash
		if(exists $global_start_hash{$pair1_trim}) {
                	die  "$_ has already been assigned to a START";
                }
                $global_start_hash{$pair1_trim} = $pair2_trim;
	} else {
		# if pair1 doesn't have START, store in the global end hash
		# to check if the contig already exists in the hash	
		if(exists $global_end_hash{$pair1_trim}) {
                        die  "$_ has already been assigned to a START";
		}
		$global_end_hash{$pair1_trim} = $pair2_trim;
		}
	if($pair2 =~ /START/) {
                # to check if the contig already exists in the hash
                if(exists $global_start_hash{$pair2_trim}) {
                        die  "$_ has already been assigned to a START";
                }
                $global_start_hash{$pair2_trim} = $pair1_trim;
        } else {
                # if pair2 doesn't have START, store in the global end hash
                # to check if the contig already exists in the hash
                if(exists $global_end_hash{$pair2_trim}) {
                        die  "$_ has already been assigned to a START";
                }
                $global_end_hash{$pair2_trim} = $pair1_trim;
                }
	} 

# run through the start hash and populate the seen hash and orient hash
my %global_seen_hash = ();
my %global_orient_hash = ();

foreach my $start_key (keys %global_start_hash) {
	# if you see the key in the start hash, then put 1 in the seen hash
	$global_seen_hash{$start_key} = 1;
	# also put E into the orient hash
	$global_orient_hash{$start_key} = 'E';
}
		
# now run through the end hash and update / populate the seen hash
#print $value, "\t", $global_start_hash{$value}, "\n";
# if you see the key in the end hash, then change the 1 into the seen hash
foreach my $end_key (keys %global_end_hash) {
	if(exists $global_seen_hash{$end_key}) {
		$global_seen_hash{$end_key} = 0;
	} else	{		
		# and if not put 1 into the seen hash
		$global_seen_hash{$end_key} = 1;
	}	
	# also put S into the orient hash
	$global_orient_hash{$end_key} = 'S';
}
# open fasta file for all the contigs for reading
my $seqio = Bio::SeqIO->new(-file => $global_options->{'fasta'}, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) {
	my $string = $seq->seq;
	my $header = $seq->id;
	$global_seq_hash{$header} = $string;
}

# open file for writing contigs
my $fasta_out = $global_options->{'fasta'}.".scaffolded";
open fasta_out, ">", $fasta_out || die $!;

my $contig_id = 0;

foreach my $seen_key (keys %global_seen_hash) {
	# value of 1 equals to ends of scaffolds, so look for these guys
	if($global_seen_hash{$seen_key} eq "1") {
		$contig_id++;
		print ">New_Contig_$contig_id\n";
        my $contig = $seen_key;
        print fasta_out ">Contig_$contig_id\n";
		my $new_contig;
		my $STF = 0;
		# this is just to count the number of seen
		# now to look for orientation of these end scaffolds in the orient hash
		# populate a table to assign the contigs to STF 1 or 0. 0 -> ST, 1 -> E 
		if($global_orient_hash{$contig} eq "E") {
			$STF = 1;
		} else {
			$STF = 0;
		}
		my $scaff_done = 0;
		while(0 == $scaff_done) {
			if(exists $global_seen_hash{$contig}) {
				$global_seen_hash{$contig} = 0; # make sure we don't get this guy twice
			}
			# if the STF is 1, means we need to reverse complement the seq
			if ($STF eq "1") {
				#print revcomp contig
				my $fa_seq = revcompl($global_seq_hash{$contig});	
				print fasta_out $fa_seq , "NNNNNNNNNNNNNNNNNNNNNNNNN";
				print "$contig\treverse\n";
				# look for the next contig in start hash
				# store in new contig
				if(exists $global_start_hash{$contig}) {
					$new_contig = $global_start_hash{$contig};
				}
				else
				{
					$scaff_done = 1;
					print fasta_out "\n";
				}
			} else {
				my $fa_seq = $global_seq_hash{$contig};
				#print reg contig
				print fasta_out $fa_seq, "NNNNNNNNNNNNNNNNNNNNNNNNN";
				print "$contig\tnormal\n";
				#look for contig in end hash
                                if(exists $global_end_hash{$contig}) {
                                        $new_contig = $global_end_hash{$contig};
                                }
                                else
                                {
                                        $scaff_done = 1;
					print fasta_out "\n";
                                }
			}
			if(0 == $scaff_done) {
				my $new_set = 0;
				if (exists $global_start_hash{$new_contig}) {
					if($global_start_hash{$new_contig} eq $contig) {
						$new_set = 1;
						$STF = 0;
					}
				}
				if(0 == $new_set)
				{
					if (exists $global_end_hash{$new_contig}) {
						if($global_end_hash{$new_contig} eq $contig) {
							$STF = 1;
						}
					}
				}
				# STF is set, now we can forget the old contig
				$contig = $new_contig;
			}
		}	
	}
}

#my $lone_contigs;
#while (<gv_fh>) {
#        chomp $_;
#	next if ($_ =~ /digraph/);
#	my ($lone_contigss) =~ /"(contig\d+)"/;
	#$_ =~ s/"//g;
        #$_ =~ s/\[.*\];//g;
        #$_ =~ s/^\s+//;
        #my @fields = split(/ -> /, $_);
	#my ($lone_contigs) = $fields[0] =~ /(contig\d+)/;
#}
#print "$lone_contigs\n";

sub fastaCut {
    #-----
    # Cut up a fasta sequence
    #
    my ($string, $prot, $line_wrap) = @_;

    # translate if need be
    if(0 != $prot)
    {
        my $codon_table = Bio::Tools::CodonTable -> new ( -id => $prot );
        $string = $codon_table->translate($string);
    }

    # wrap the line if need be
    if(0 != $line_wrap)
    {
        my $return_str = "";
        my $len = length $string;
        my $start = 0;
        while($start < $len)
        {
            $return_str .= substr $string, $start, $line_wrap;
            $return_str .="\n";
            $start += $line_wrap;
        }
        return $return_str;
    }
    return "$string\n";
}

sub print_seq{
    my ($name_ref, $seq_ref, $qual_ref, $fh) = @_;
    my $seq = $$seq_ref;
    if(defined $global_options->{'wrap'})
    { 
        if(defined $global_options->{'protein'})
        {
            $seq = fastaCut($seq, $global_options->{'protein'}, $global_options->{'wrap'});
        }
        else
        {
            $seq = fastaCut($seq, 0, $global_options->{'wrap'});
        }
    }
    elsif(defined $global_options->{'protein'})
    {
        $seq = fastaCut($seq, $global_options->{'protein'}, 0);
    }
    else
    {
        $seq .= "\n";
    }

    if (defined $$qual_ref)
    {
        # fastq file
        print $fh "@".$$name_ref."\n".$seq."+".$$name_ref."\n".$$qual_ref."\n";
    }
    else
    {
        print $fh ">".$$name_ref."\n".$seq;
    }
}
sub revcompl {
	my ($seq) = @_;
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return scalar reverse $seq;
}
#####################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "gv|g:s", "fasta|f:s");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );

    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    if(!exists $options{'gv'} || !exists $options{'fasta'} ) { printParamError ("both -g and -f must be specified"); }

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
# checkAndRunCommand("ls", {
#                          -a => ""
#                          }, 
#                          WARN_ON_FAILURE);

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
Copyright (C) Michael Imelfort

This program comes with ABSOLUTELY NO WARRANTY;
This is free software, and you are welcome to redistribute it
under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

scaffolder.pl

=head1 COPYRIGHT

copyright (C) Fauzi Haroon, Connor Skennerton

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 DESCRIPTION

Given a Graphviz (.gv) file and a contig file, scaffold the contigs together

=head1 SYNOPSIS

scaffolder.pl  [-help|h] -g|gv <gv_file> -f|fasta <fasta_file>

[-help -h]                   Displays basic usage information
-g -gv FILE                  Graphviz file
-f -fasta FILE               File containing contigs
=cut
