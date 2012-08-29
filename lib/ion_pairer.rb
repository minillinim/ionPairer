require 'graphviz'
require 'bio-logger'
require 'yaml'

require 'ion_pairer/bl2seq_runner'


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


module IonPairer
  # Assemblies are made of scaffolds
  class Assembly < Array
  end
  
  # Scaffolds are made of components (contigs or unresolveds)
  class Scaffold < Array
    # Return a string YAML representation of this scaffold in a format s
    def scaffolder_yaml
      all = []
      each do |component|
        if component.kind_of?(IonPairer::Contig)
          component.name.force_encoding 'UTF-8' #otherwise the YAML outputter outputs them as ASCII-8BIT, which is encoded in binary, no good.
          name_bit = {'source' => component.name}
          base = {
            'sequence' => name_bit
          }
          if component.reverse
            base['sequence']['reverse'] = true
          end
          all.push base
          
        elsif component.kind_of?(IonPairer::Unresolved)
          all.push 'unresolved' => {'length' => 25}
          
        else
          raise "Unexpected component type in assembly: #{component.class}"
        end
      end
      
      return all.to_yaml
    end
  end
  
  class Component
  end
  
  class Contig < Component
    attr_accessor :name, :reverse
    
    def initialize
      @reverse = false
    end
    
    def to_yaml
      p @name
      yaml = {:name => @name}
      if reverse
        yaml[:reverse => 'true']
      end
    end
  end
  
  class Unresolved < Component
    attr_accessor :length
    
    def initialize
      @length = 25
    end
  end
  
  class GraphvizReader
    class GraphvizContig
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
    
    LOG_NAME = 'ionpairer'
    def log
      if Bio::Log::LoggerPlus[LOG_NAME].nil?
        Bio::Log::LoggerPlus.new(LOG_NAME)
      end
      Bio::Log::LoggerPlus[LOG_NAME]
    end
    
    
    # Return an Assembly object, filled with the assembly as defined in the graphviz file
    def read_assembly(graphviz_file)
      graph = GraphViz.parse(graphviz_file)
      log.info "Parsed in graph file containing #{graph.node_count} nodes (i.e. #{graph.node_count/2} contigs) and #{graph.edge_count} edges (i.e. #{graph.edge_count-graph.node_count/2} links between contigs)"
      
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
      
      # Find nodes that have no edges. These are the ones that will begin or end each scaffold
      stoppers = node_edges.select{|node, edges| edges.length==1}.keys
      raise unless stoppers.length.modulo(2) == 0 #is this demonstrably true by this point regardles of input? Maybe, still being paranoid.
      log.debug "Ends of scaffolds: #{stoppers.join(", ")}"
      
      assembly = Assembly.new
      
      
      add_contig = lambda do |start, stop|
        contig = Contig.new
        
        gv_contig = GraphvizContig.new
        gv_contig.node1_name = start
        gv_contig.node2_name = stop
        contig.name = gv_contig.contig_name
        if gv_contig.orientation == '-'
          contig.reverse = true
        end
        contig
      end
      
      
      while !stoppers.empty?
        scaffold = Scaffold.new
        assembly.push scaffold
        
        start = stoppers[0]
        stoppers.delete start
        edges = node_edges[start]
        raise unless edges.length == 1
        edge = edges[0]
        stop = edge.other_node(start)
        node_edges[stop].delete edge
        
        scaffold.push add_contig.call start, stop

        # For each new contig until we've run out
        while !stoppers.include?(stop)
          scaffold.push Unresolved.new
          
          start = node_edges[stop][0].other_node(stop) #node that connects to the trailing end of the last contig
          last_trailing_node = stop
          edges = node_edges[start].reject{|edge| edge.node_one == last_trailing_node or edge.node_two == last_trailing_node}
          
          raise unless edges.length == 1
          edge = edges[0]  
          stop = edge.other_node(start)
          node_edges[stop].delete edge
          
          scaffold.push add_contig.call start, stop
          #log.debug "Finished with #{stop}, which now has #{node_edges[stop].length} edges"
        end
        stoppers.delete stop
      end
      
      return assembly
    end
  end
end
