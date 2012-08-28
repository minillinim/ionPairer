#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'bio'

$:.unshift File.join(File.dirname(__FILE__),'..','lib')
require 'ion_pairer'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :overlap_length => 100,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -d <dot_file> -s <fasta_file>
      
      Take an assembly, and works out if there is any overlapping ends between the contigs that got scaffolded together. \n\n"
      
    opts.on("-d", "--dot-file DOT_FILE_PATH", "assembly file in graphviz dot format [required]") do |arg|
      options[:dot_file] = arg
    end
    opts.on("-s", "--fasta-file CONTIGS_FASTA_FILE", "fasta file of contig sequences [required]") do |arg|
      options[:fasta_file] = arg
    end
    opts.on("-o", "--overlap INTEGER", "How far into the sequences to look for overlaps into [default #{options[:overlap_length]}]]") do |arg|
      options[:overlap_length] = arg.to_i
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 0 or options[:dot_file].nil? or options[:fasta_file].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
  # Read in the contigs
  contigs = {} #hash of contig names to sequences
  Bio::FlatFile.foreach(options[:fasta_file]) do |s|
    first_name = s.definition.split(/\s/)[0]
    raise if contigs[first_name]
    contigs[first_name] = s.seq
  end
  log.info "Cached #{contigs.length} sequences, e.g. #{contigs.to_a[0][0]} => #{contigs.to_a[0][1].gsub(/^(.{10}).*/,'\1')}..."
  
  # Read assembly in
  assembly = IonPairer::GraphvizReader.new.read_assembly options[:dot_file]
  
  # Blast each contig end against each other
  blast = Bio::BlastPlus::Runner.new
  scaffold_counter = 0
  assembly.each do |scaffold|
    (0...scaffold.length-2).each do |contig_index|
      contig1 = scaffold[contig_index]
      unresolved = scaffold[contig_index+1]
      contig2 = scaffold[contig_index+2]
      next if contig1.kind_of?(IonPairer::Unresolved) #only can link between contigs, duh
      
      unless contig1.kind_of?(IonPairer::Contig) and unresolved.kind_of?(IonPairer::Unresolved) and contig2.kind_of?(IonPairer::Contig)
        raise "Unexpected assembly near #{contig1.inspect}, #{unresolved.inspect}, #{contig2.inspect}, I expect contig,unresolved,contig"
      end
      
      # Look for blast overlap between the end of contig1 and the start of contig2
      contig1_seq = contigs[contig1.name]
      seq_one = contig1_seq[contig1_seq.length-options[:overlap_length]...contig1_seq.length]
      contig2_seq = contigs[contig2.name]
      seq_two = contig2_seq[0...options[:overlap_length]]
      results = blast.bl2seq(seq_one, seq_two)
      
      unless results.empty?
        puts
        contig1_start_or_stop = contig1.reverse ? 'start' : 'end'
        contig2_start_or_stop = contig2.reverse ? 'end' : 'start'
        r = results[0]
        puts "Possible overlap on scaffold #{scaffold_counter} between contigs #{contig1.name} (#{contig1_start_or_stop}) and #{contig2.name} (#{contig2_start_or_stop}), %ID #{r.pident}, Length #{r.length}, #{r.qstart}-#{r.qend} vs #{r.sstart}-#{r.subject_end}"
        print '   first '
        puts seq_one
        print '  second '
        puts seq_two
      end
    end
    scaffold_counter += 1
  end
end #end if running as a script