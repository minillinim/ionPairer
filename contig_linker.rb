#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'graphviz'
require 'bio'
require 'pp'
require 'csv'
# Use a different CSV parser for Ruby 1.8 than 1.9+
if CSV.const_defined? :Reader
  require 'fastercsv'
  csv = FasterCSV
else
  csv = CSV
end


class ContigLink
  START_OF_CONTIG = 'contig_start'
  END_OF_CONTIG = 'contig end'
  INCONSISTENT_POSITION_AND_DISTANCE = 'inconsistent position and distance on contig ends'
  NOT_NEAR_EITHER_END = 'not near either end of the contig'
  
  attr_accessor :position1, :position2,
    :direction1, :direction2
    
  def log
    Bio::Log::LoggerPlus['contig_linker']
  end
    
  # Determine ends that the links are. They must be
  # * within the max_distance from the ends of the contig
  # * have consistent direction i.e. can't be pointing towards the inside of the contig, must be pointing out
  def predict_link_ends(max_distance, contig1_length, contig2_length)
    return [
      ContigLink.predict_link_direction(max_distance, contig1_length, @position1, @direction1),
      ContigLink.predict_link_direction(max_distance, contig2_length, @position2, @direction2),
      ]
  end
  
  # Given info about a contig and a read direction and position and max distance from the end
  # return something about the way it 
  def self.predict_link_direction(max_distance, contig_length, position, direction)
    near_start = position < max_distance
    near_end = contig_length-position < max_distance
    
    if near_start and near_end
      # a small contig. Just take it based on distance, not possible to be inconsistent
      if direction == '+'
        END_OF_CONTIG
      elsif direction == '-'
        START_OF_CONTIG
      else
        raise "Unexpected direction #{direction.inspect}"
      end
    else
      # a longer contig
      if near_start
        case direction
        when '+' then INCONSISTENT_POSITION_AND_DISTANCE
        when '-' then START_OF_CONTIG
        else raise "Unexpected direction #{direction.inspect}"
        end
      elsif near_end
        case direction
        when '+' then END_OF_CONTIG
        when '-' then INCONSISTENT_POSITION_AND_DISTANCE
        else raise "Unexpected direction #{direction.inspect}"
        end
      else
        NOT_NEAR_EITHER_END
      end
    end
  end
  
  # Given 2 contig lengths, and the joining between those contigs, estimate how many Ns are needed between them
  # Simply mean - distance_from_contig1_end - distance_from_contig2_end
  def estimate_distance_between_contigs(contig_lengths, contig_ends, mean_insert_size)
    compute_distance = lambda do |position, contig_length, contig_end|
      distance_from_end = nil
      if contig_end==START_OF_CONTIG
        distance_from_end = position
      elsif contig_end==END_OF_CONTIG
        distance_from_end = contig_length-position
      else
        raise "programming error"
      end
      distance_from_end
    end
    
    distance_from_end1 = compute_distance.call(@position1, contig_lengths[0], contig_ends[0])
    distance_from_end2 = compute_distance.call(@position2, contig_lengths[1], contig_ends[1])
    return mean_insert_size-distance_from_end1-distance_from_end2
  end
end

# The set of links between two contigs
class ContigLinkageSet
  attr_accessor :contig1_name, :contig2_name
  attr_accessor :contig1_length, :contig2_length
  attr_accessor :links

  def log
    Bio::Log::LoggerPlus['contig_linker']
  end
  
  def classified_abuttings_hash(max_distance)
    end_identifier_hash = {}
    @links.each do |link|
      # near the start and correct direction for the first contig?
      # which end is it near
      contig_length1 = @contig_lengths[contig1]
      contig_length2 = @contig_lengths[contig2]
      if contig_length1.nil?
        raise "Unexpected lack of @contig_lengths[#{contig1}] or @contig_lengths[#{contig2}]"
      end
      closest_ends = link.predict_link_ends(max_distance, @contig1_length, @contig2_length)
      
      key = closest_ends[0], closest_ends[1]
      end_identifier_hash[key] ||= []
      end_identifier_hash[key].push link
    end
    return
  end
  
  # Work out the maximal direction of the linking - i.e. which ends should be joined.
  # Return nil if there is no winner
  def predict_best_abutting(max_distance)
    # classify each link as start to start, etc.
    end_identifier_hash = classified_abuttings_hash(max_distance)
    log.debug "For #{@contig1_name} and #{@contig2_name} found the following types of links: #{end_identifier_hash.keys.inspect}"
    
    # Remove presumably spurious links
    possibilities = [
      ContigLink::START_OF_CONTIG,
      ContigLink::END_OF_CONTIG,
    ]
    end_identifier_hash.select! do |key, links|
      possibilities.include?(key[0]) and possibilities.include?(key[1])
    end
    
    # sort the possibilities in reverse order
    sorts = counts.to_a.sort{|a,b| -(a[1].length<=>b[1].length)}
    log.debug "After sorting, got sorted pairs #{sorts.inspect} for #{contig1} and #{contig2}"
    
    return nil if sorts.length == 0
    
    # look for tied winners to be anal, return no links if this is the case
    if sorts[0][1]==sorts[1][1] and
      log.warn "Contig link #{@contig1_name} and #{@contig2_name} have confusing orientation statistics (#{sorts.inspect}), not putting into the graphviz file"
      return nil
    end
      
    return sorts[0][0], sorts[0][1], abutting_links
  end
  
  def self.distance_between_contigs(contig_lengths, abuttings, abutting_links, mean_insert_size)
    # simply the mean of the individual estimates
    abutting_links.collect {|link|
      link.estimate_distance_between_contigs(contig_lengths, abuttings, mean_insert_size)
      }.inject{|total, cur| total+=cur}.to_f/abutting_links.length
  end
end

# A class for sets of sets of contig linkages
# A hash of [contig_name1, contig_name2] => Array of ContigLink's.
class ContigLinkSet < Hash
  # Hash of contig names to lengths (lengths being integers)
  attr_accessor :contig_lengths
  
  def add_contig_set(contig1_name, pair1_position, pair1_direction,
                     contig2_name, pair2_position, pair2_direction)
    
    # order by contig names, otherwise hashing is harder
    one = [contig1_name, pair1_position, pair1_direction]
    two = [contig2_name, pair2_position, pair2_direction]
    array = [one, two].sort{|a1, a2| a1[0] <=> a2[0]}
    
    relationship = ContigLink.new
    relationship.position1 = array[0][1]
    relationship.position2 = array[1][1]
    relationship.direction1 = array[0][2]
    relationship.direction2 = array[1][2]
    
    key = [array[0][0], array[1][0]]
    self[key] ||= []
    self[key].push relationship
  end
  
  def log
    Bio::Log::LoggerPlus['contig_linker']
  end
  
  # output to the graphviz file iff each side of the pair maps close enough to the ends of the contigs
  # and is in a direction that makes sense like this.
  #
  # options:
  # :mean_insert_size: mean insert size, to estimate amount of sequence between contig lengths [required]
  # :min_links: minimum number of links to bother drawing on the graph [default 3]
  def generate_graphviz(max_distance, options = {})
    options[:min_links] ||= 3
    raise unless options[:mean_insert_size]
    raise unless options[:contig_lengths]
    
    # create initial graphviz nodes
    # each contig is 2 nodes - one for the start and one for the end
    # that way it is easy to tell which end it is referring to
    get_start_node_name = lambda{|contig| "#{contig}START"}
    get_end_node_name = lambda{|contig| "#{contig}END"}
    graphviz = GraphViz.new(:G, :type => :graph)
    graphviz.node_attrs[:shape] = :point
    
    # Create contig start/end nodes and add edges between the start and the end
    @contig_lengths.each do |contig, length|
      start_name = get_start_node_name.call(contig)
      stop_name = get_end_node_name.call(contig)
      graphviz.add_nodes(start_name)
      graphviz.add_nodes(stop_name)
      color = 'blue'
      color = 'grey' if length < max_distance*2
      graphviz.add_edges(start_name, stop_name, :style => "setlinewidth(4)", :label => contig, :color => color)
    end
    
    
  
    # Iterate through contig pairs, outputing as necessary
    each do |contig_name_pair, links|
      linkset = ContigLinkageSet.new
      linkset.contig1_name = contig_name_pair[0]
      linkset.contig2_name = contig_name_pair[1]
      linkset.contig1_length = @contig_lengths[contig_name_pair[0]]
      linkset.contig2_length = @contig_lengths[contig_name_pair[1]]
      linkset.links = links
      
      # Find the most likely form of linkage between the contigs
      abutting1, abutting2, abutting_links = linkset.predict_best_abutting
      

      contig1 = contig_name_pair[0]
      contig2 = contig_name_pair[1]
      log.debug "evaluating #{contig1} and #{contig2}" unless log.nil?
      
      # Remove those with insufficient numbers of links
      unless abutting_links.length >= options[:min_links]
        log.debug "Not putting in linkage between #{contig1} and #{contig2}, as there was insufficient links (#{max[1]} found)"
        next
      end
      
      # create a new edge between the appropriate sides.
      log.debug "Assigning the link between #{contig1} and #{contig2}, as max #{max[1]}"
      
      distance_between_contigs = linkset.distance_between_contigs(
                                                                   [linkset.contig1_length, linkset.contig2_length],
                                                                   [abutting1, abutting2],
                                                                   abutting_links,
                                                                   mean_insert_size
                                                                   )
      log.debug "Predicting a distance of #{distance_between_contigs} between #{contig1} and #{contig2}"
      
      assign_name = lambda do |contig_name, ending|
        if ending == ContigLink::START_OF_CONTIG
          from_node = get_start_node_name.call(contig_name)
        elsif ending == ContigLink::END_OF_CONTIG
          from_node = get_end_node_name.call(contig_name)
        else
          raise "programming error"
        end
      end
      from_node = assign_name.call(contig1, max[0][0])
      to_node = assign_name.call(contig2, max[0][1])
  
      graphviz.add_edges(from_node, to_node, :label => "#{max[1]}links_dist#{distance_between_contigs}")
    end
    
    #return the entire graph
    return graphviz
  end
end





if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = SCRIPT_NAME.gsub('.rb','')
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :max_distance => 4000,
    :min_links => 3,
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -l <linkage_file>
      
      Take an input linkage file, the 'all_links.csv' file, and output:
        
        1. a graphviz file <linkage_file>.gv containing linkages between the contigs\n\n"
      
    opts.on("-l", "--linkage-file PATH", "linkage file from ionPairer.pl [required]") do |arg|
      options[:linkage_file] = arg
    end
    opts.on("-f", "--reference-fasta-file PATH", "fasta file of reference sequences mapped to the sam file [required]") do |arg|
      options[:reference_fasta_file] = arg
    end
    opts.on("-u", "--mean MEAN_INSERT_LENGTH", "Average insert size, to estimate amount of sequence in contig joins [required]") do |arg|
      options[:mean_insert_size] = arg.to_i
      raise unless options[:mean_insert_size] >= 0
    end
    opts.on("-m", "--max-distance ARG", "I'm looking for pairs that map near the ends of contigs. How far is too far from the end, in base pairs? [default #{options[:max_distance]}]") do |arg|
      options[:max_distance] = arg.to_i
      raise unless options[:max_distance] > 0
    end
    opts.on("-L", "--min-linkages ARG", "How many linkages are sufficient to call a join [default #{options[:min_links]}]") do |arg|
      options[:min_links] = arg.to_i
      raise unless options[:min_links] >= 0
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 0 or options[:linkage_file].nil? or options[:reference_fasta_file].nil? or options[:mean_insert_size].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  
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
  
  
  # read in the linkage file
  contig_linkset = ContigLinkSet.new
  ignored_as_the_same_contig = 0
  csv.open(options[:linkage_file], :col_sep => "\t", :headers => true).each do |row|
    # contig1_name  position1  mapping_quality1  direction1  contig2_name  position2  mapping_quality2  direction2
    #70	25937	149	1	70	24818	17	0
    #54	32081	20	1	54	29539	142	0
    #86	8990	9	1	86	7389	153	0
    #88	14033	20	1	88	12003	24	0
    raise "Unexpected number of columns in all_links file" if row.length != 8
    
    if row[0]==row[4]
      ignored_as_the_same_contig += 1
    else
      # convert some columns to integers
      [1, 5].each do |i|
        row[i] = row[i].to_i
      end
      
      # convert direction into '+' or '-'
      [3, 7].each do |i|
        if row[i]=='0'
          row[i] = '+'
        elsif row[i]=='1'
          row[i] = '-'
        else
          raise "unexpected direction"
        end
      end
  
      contig_linkset.add_contig_set(
        row[0], row[1], row[3],
        row[4], row[5], row[7]
      )
    end
  end
  eg_key = contig_linkset.keys[0]
  eg = contig_linkset[eg_key]
  log.info "Read in linkages between #{contig_linkset.length} pairs of contigs. For instance between #{eg_key[0]} and #{eg_key[1]} there are #{eg.length} links"
  
  contig_linkset.contig_lengths = contig_lengths
  graphviz = contig_linkset.generate_graphviz(options[:max_distance], {:min_links => options[:min_links], :mean_insert_size => options[:mean_insert_size]})
  
  # print outputs
  graphviz.output :dot => "#{options[:linkage_file]}.dot"
  graphviz.output :png => "#{options[:linkage_file]}.png", :use => :neato
  graphviz.output :svg => "#{options[:linkage_file]}.svg", :use => :neato
  
end #end if running as a script





