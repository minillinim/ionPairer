          
               d8b                   8888888b.          d8b                         
               Y8P                   888   Y88b         Y8P                         
                                     888    888                                     
               888  .d88b.  88888b.  888   d88P 8888b.  888 888d888 .d88b.  888d888      
               888 d88""88b 888 "88b 8888888P"     "88b 888 888P"  d8P  Y8b 888P"   
               888 888  888 888  888 888       .d888888 888 888    88888888 888     
               888 Y88..88P 888  888 888       888  888 888 888    Y8b.     888     
               888  "Y88P"  888  888 888       "Y888888 888 888     "Y8888  888     
                                                                                                        
           
                 Scrimpy way to scaffold an assembly using ion torrent reads.

# Overview

ionPairer is a library for scaffolding assemblies somewhat manually using IonTorrent mate pair reads.

A typical workflow would be to start with an assembly generated from single-ended or paired-end (insert size ~300),
and mate pair data generated through ionTorrent. The original data files are not necessary

# Installation

It's a bit complicated. You need several things:

* Perl
* Ruby
* GraphViz
* Perl dependencies (some core + Bio::SeqIO, Bio::DB::Sam)
* Ruby dependencies

The few Ruby depencies can be installed by changing directory to the base directory
of ionPairer and running

```sh
$ bundle install
```
    
Note: those people at ACE need only ```module load ionPairer```. 

# Example usage

## First map mate pair reads against the assembly

For IonTorrent read and certain versions of BWA, it may be helpful to run bwa twice,
once for the forward reads and once for the reverse reads.

ie.

Take your existing assembly (or binned contigs from a metagenome), then
run ```bwa-sw``` twice (once for each of the forward and reverse reads) to map your
ion MP reads against your contigs. e.g.

```sh
$ bwa index my_assembly.fna    #to create database for bwasw to search against
$ bwa bwasw my_assembly.fna forward_mates.fna >forward_matesVmy_assembly.sam
$ bwa bwasw my_assembly.fna reverse_mates.fna >reverse_matesVmy_assembly.sam
```

For Illumina-type reads, you can just run BWA once with both pairs. Perhaps
makeSam.py would be helpful here: http://www.github.com/minillinim/mikesbioscripts

## Then run ionPairer:
```sh
$ ionPairer.pl -sam1 forward_matesVmy_assembly.sam -sam2 reverse_matesVmy_assembly.sam -s -w ion_pairer_outputs
```
or

```sh
$ ionPairer.pl -1 pairedReads_Vmy_assembly.bam -w ion_pairer_outputs
```

In that output directory ```ion_pairer_outputs``` there should be the following files:

* ```forward_matesVmy_assembly.sam.unpaired.csv``` Reads where only one end mapped to a contig
* ```forward_matesVmy_assembly.sam.pcr_duplicates.csv``` Reads removed from further analysis as they were judged to be PCR duplicates
* ```forward_matesVmy_assembly.sam.paired.csv``` Reads where both ends mapped onto one contig
* ```forward_matesVmy_assembly.sam.error_paired.csv``` Reads where both ends mapped, but erroneously due to insert size or relative orientation
* ```forward_matesVmy_assembly.sam.unique_links.csv``` All pairs of mate pairs that span between two contigs and pass the pcr duplicates filter  
* ```forward_matesVmy_assembly.sam.short_links.csv``` Pairs where one of the contigs was shorter than the insert size
* ```forward_matesVmy_assembly.sam.unique_links.error_links.csv``` Pairs which link two contigs, but erroneously due to insert size, position or relative orientation (these may indicate chimeras in the contigs)
* ```forward_matesVmy_assembly.sam.unique_links.filtered_links.csv``` Subset of unique pairs which are neither short nor erroneous

There'll also be three graphviz-related files:

* ```forward_matesVmy_assembly.sam.unique_links.gv``` The dot file (GraphViz file) used to specify the links between contigs
* ```forward_matesVmy_assembly.sam.unique_links.png``` A png picture representation of the dot file. The colours are meaningful. Specifically red contigs are longer than 2x 
* ```forward_matesVmy_assembly.sam.unique_links.svg``` An svg picture representation of the dot file

## The next step is manual inspection of the scaffolds. 
The fun bit! In 
```forward_matesVmy_assembly.sam.all_links.png``` (or the svg file)
is a representation of the links that have been made, and the number of mate-pairs that 
agree with that linking.

To modify the scaffolds, modify the ```gv``` file (maybe best to also copy it to a new file ```forward_matesVmy_assembly.sam.all_links.manually_modified.gv```). For instance if this link is no good, for whatever reason:

```
contig00056END -- contig00073START [label="3links", ...
```
Then simply delete the entire line. It may be advisable to save this modified version as a different file so you can use ```diff``` or ```meld``` to work out what has changed between the original and manually curated versions.

Afterwards, rerun graphviz, which will make a new png/svg file for you, like so
```sh
$ neato -Tpng forward_matesVmy_assembly.sam.all_links.manually_modified.gv >forward_matesVmy_assembly.sam.all_links.manually_modified.png
```
(or ```-Tsvg``` for svg output.)

## Once you are happy with the manual assembly and want to do the actual scaffolding
Once the ```.gv``` file has been editted to your satisfaction, you have 2 options.

1. Use the ```scaffolder.pl``` script to create the scaffolds. This script will orient the contigs (forward or reverse complement) and scaffold, based on the dot file. There is a small bug though, where single unlinked contigs have to be manually added into the final output ```.scaffolded```.

Example usage:
```sh
$ scaffolder.pl -g edited.gv -f contigs.fna > scaffold_map.fna
```

2. Use the ```dot_to_scaffolder_yaml.rb``` script to convert the dot file to a YAML file, in a particular [format](http://next.gs). This YAML file can be further modified, for example to overlap contig ends. This script is in the tools folder.

Example usage:
```sh
$ mkdir scaffolds
$ dot_to_scaffolder_yaml.rb -g edited.gv -d scaffolds
```

## Joining contigs across unresolved parts
Sometimes if two contigs reside next to each other in a scaffold, then they in reality overlap - there is no Ns between them. The tool ```blast_contig_ends.rb``` automates the process of looking for these using a simple blast-based method (blast+, that is). If the output of this script shows that, indeed, some contigs do overlap, then modify the output from ```dot_to_scaffolder_yaml.rb```, deleting the relevant ```unresolved``` bit and specify where exactly the overlap occurs by using ```start``` and ```stop```. See http://next.gs/man/scaffolder-format/ for more information. When finished use the ```scaffolder``` tool from http://next.gs to generate the final, assembled, fasta file.

Example usage:
```sh
$ blast_contig_ends.rb -d edited.gv -s contigs.fna
 INFO blast_contig_ends: Cached 117 sequences, e.g. contig00001 => CGCTTGGCCT...

Possible overlap on scaffold 3 between contigs contig00097 (start) and contig00035 (start), %ID 92.96, Length 199, 1-199 vs 199-2
   first TCATCAAATTACGTAGGAGCACGGGAGGTAATAAGACATACACAAGGAAAAAaCAGATGACGAAGGAGAAACACGCTGCGCTATTATAGTTCGGGCCTACCCGAGAATTGTGGAGGGCTCTCTACTCACTTAAGTGAAAACGAAATTGCGACTGAAAAAGCACTCCCAGTAATCTGCTAACAGCATACGGCTAGAAATTG
  second TAATTTCTAGCCGCATGCTGTTAGCAGATTACTGCGAGTGCTTTTTCAGTCGCAATTtCGTTTTCACTTAAGTGAGTAGAGGGCCTGTCACAGTTCTCGGGTAGGCCCGAGCTATAATAGCGCAGCGTGTTTCTCCTTCGTCATCTGTTTTTTCCTTGTGTATGCCTTATTACCTTCCGTGCGCTACTCAATTTGATGAT

(manually modify scaffolds/scaffold3.yml to indicate the overlap)

$ scaffolder sequence scaffolds/scaffold3.yml contigs.fna >scaffolds.scaffold3.fna
```
That last step needs to be repeated for each scaffold (potentially in a bash for loop).

# Administration

Project home page, info on the source tree, documentation, issues and how to contribute, see http://github.com/minillinim/ionPairer

This software is currently unpublished.

Copyright © 2012 Michael Imelfort, Ben Woodcroft, Fauzi Haroon. See LICENSE.txt for further details.
