#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;

# globals
my %global_start_hash = (); #contig pairs that begin with START
my %global_end_hash = (); #contig pairs that begin with END
my ($gv_fh, $fasta_in) = @ARGV;

# First read all the sequences into a hash
# use contigID -> seq
my %global_seq_hash = ();
my %global_seq_printed_hash = (); # in this hash means printed to file
my $global_contig_id = 1; #increment the number when printing, so that it is contig1 ... contig2 ...

# need to read file
# populate start and end hashes
open gv_fh, "<", $gv_fh || die $! ;
while (<gv_fh>) {
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
	my ($pair1) = $fields[0] =~ /(contig.*)/;
        my ($pair2) = $fields[1] =~ /(contig.*)/;
	# trim the name to remove the START and then store in the global start hash
        my ($pair1_trim) = $pair1 =~ /(contig\d+)/;
        my ($pair2_trim) = $pair2 =~ /(contig\d+)/;

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
open fasta_in, "<", $fasta_in || die $! ;
my $seqio = Bio::SeqIO->new(-file => $fasta_in, '-format' => 'Fasta');
while(my $seq = $seqio->next_seq) {
	my $string = $seq->seq;
	my $header = $seq->id;
	$global_seq_hash{$header} = $string;
}

# open file for writing contigs
my $fasta_out = "$ARGV[1].scaffolded";
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

sub revcompl {
	my ($seq) = @_;
	$seq =~ tr/ACGTacgt/TGCAtgca/;
	return scalar reverse $seq;
}
