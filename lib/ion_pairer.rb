require 'graphviz'
require 'bio-logger'

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
  class Assembly < Array
    
  end
  
  class Scaffold < Array
    
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
    def self.read_assembly(graphviz_file)
      
    end
  end
end
