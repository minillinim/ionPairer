#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'graphviz'
require 'bio'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :default_distance => 25,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} <dot_file> <fasta_file>
      
      Takes a dot file originally output from ionPairer.pl/contig_linker.rb, after it has been ironed out (all scaffolds are linear).
      
      Outputs YAML suitable for input to the scaffolder rubygem (http://next.gs) which can then be used to build contigs (potentially after further manual adjustment of the scaffolding process e.g. resolving overlapping contig ends.)
      \n\n"
      
    opts.on("-f", "--reference-fasta-file PATH", "fasta file of reference sequences mapped to the sam file [required]") do |arg|
      options[:reference_fasta_file] = arg
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 1 or options[:reference_fasta_file].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  log.warn "Note that this script currently cannot handle fully circular inputs (i.e. cyclic graphs). To handle these, break it up and then modify the output file manually"
  graph = GraphViz.parse(ARGV[0])
  log.info "Parsed in graph file containing #{graph.node_count} nodes (i.e. #{graph.node_count/2} contigs) and #{graph.edge_count} edges (i.e. #{graph.edge_count-graph.node_count/2} links between contigs)"
  
  # read in the SAM file, caching the contig lengths
  contig_lengths = {} #hash of name to contig length
  Bio::FlatFile.foreach(options[:reference_fasta_file]) do |seq|
    name = seq.definition.split(/\s/)[0]
    if contig_lengths.key?(name)
      raise "Duplicate contig name in fasta #{name}, giving up"
    else
      contig_lengths[name] = seq.seq.length
    end
  end
  log.info "Cached #{contig_lengths.length} contig lengths from the fasta file, e.g #{contig_lengths.keys[0]} => #{contig_lengths[contig_lengths.keys[0]]}."
  
  # Find nodes that have no edges. These are the ones that will begin or end each scaffold
  node_edges = {}
  graph.each_edge do |edge|
    [edge.node_one, edge.node_two].each do |node|
      node_edges[node] ||= []
      node_edges[node].push edge
    end
  end
  bads = node_edges.reject{|node, edges| [1,2].include?(edges.length)}
  unless bads.empty?
    log.error "Unexpectedly found too many edges for some nodes e.g. #{bads[0].inspect}, giving up"
    exit 1
  end
  stoppers = node_edges.select{|node, edges| edges.length==1}.keys
  raise unless stoppers.length.modulo(2) == 0 #is this demonstrably true by this point regardles of input? Maybe, still being paranoid.
  
  
  class GraphViz::Edge
    def other_node(node_name)
      if node_name==node_one and node_name != node_two
        return node_two
      elsif node_name!=node_one and node_name == node_two
        return node_one
      else
        raise
      end
    end
  end
  
  class Contig
    attr_accessor :node1_name, :node2_name
    
    def contig_name
      if matches = @node1_name.match(/(.*)START$/)
        return matches[1]
      elsif matches = @node1_name.match(/(.*)END$/)
        return matches[1]
      else
        raise "Unexpected contig name: #{@node1_name}, can't assign a contig name to it and its pair"
      end
    end
    
    def orientation
      if @node1_name.match(/START$/) and @node2_name.match(/END$/)
        return '+'
      elsif @node2_name.match(/START$/) and @node1_name.match(/END$/)
        return '-'
      else
        raise "Unexpected contig pairing, #{node1_name} and #{node2_name}"
      end
    end
  end
  log.debug "Ends of scaffolds: #{stoppers.join(", ")}"
  
  while !stoppers.empty?
    puts '---'
    
    start = stoppers[0]
    stoppers.delete start
    edges = node_edges[start]
    raise unless edges.length == 1
    edge = edges[0]
    stop = edge.other_node(start)
    node_edges[stop].delete edge
    
    contig = Contig.new
    contig.node1_name = start
    contig.node2_name = stop
    
    puts [
      '  -',
      "  sequence:",
      "    source: #{contig.contig_name}",
    ].join("\n  ")
    if contig.orientation == '-'
      puts '      reverse: true'
    end
    
    # For each new contig until we've run out
    while !stoppers.include?(stop)
      puts [
        '  -',
        '    unresolved:',
        '      length: 25',
      ]
      
      start = node_edges[stop][0].other_node(stop) #node that connects to the trailing end of the last contig
      last_trailing_node = stop
      edges = node_edges[start].reject{|edge| edge.node_one == last_trailing_node or edge.node_two == last_trailing_node}
      #p edges
      raise unless edges.length == 1
      edge = edges[0]  
      stop = edge.other_node(start)
      node_edges[stop].delete edge
      
      contig = Contig.new
      contig.node1_name = start
      contig.node2_name = stop
      
      puts [
        '  -',
        "  sequence:",
        "    source: #{contig.contig_name}",
      ].join("\n  ")
      if contig.orientation == '-'
        puts '      reverse: true'
      end
      #log.debug "Finished with #{stop}, which now has #{node_edges[stop].length} edges"
    end
    stoppers.delete stop
  end
  
  
  
end #end if running as a script