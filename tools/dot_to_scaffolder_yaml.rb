#!/usr/bin/env ruby

require 'optparse'
require 'bio-logger'
require 'graphviz'
require 'bio'

$:.unshift File.join(File.dirname(__FILE__),'..','lib')
require 'ion_pairer'

if __FILE__ == $0 #needs to be removed if this script is distributed as part of a rubygem
  SCRIPT_NAME = File.basename(__FILE__); LOG_NAME = 'ionpairer'
  
  # Parse command line options into the options hash
  options = {
    :logger => 'stderr',
    :default_distance => 25,
    :output_directory => '.',
  }
  o = OptionParser.new do |opts|
    opts.banner = "
      Usage: #{SCRIPT_NAME} -g <dot_file> [-d <output_directory>]
      
      Takes a dot file originally output from ionPairer.pl/contig_linker.rb, after it has been ironed out (all scaffolds are linear).
      
      Outputs YAML suitable for input to the scaffolder rubygem (http://next.gs) which can then be used to build contigs (potentially after further manual adjustment of the scaffolding process e.g. resolving overlapping contig ends.)
      \n\n"
      
    opts.on("-g", "--graphviz-dot-file PATH", "dot file that defines the scaffolding [required]") do |arg|
      options[:dot_file] = arg
    end
    opts.on("-d", "--output-directory DIRECTORY", "Scaffolds are output one scaffold per file. Where do you want these files to go? [default: #{options[:output_directory]}]") do |arg|
      options[:output_directory] = arg
    end
    
    # logger options
    opts.on("-q", "--quiet", "Run quietly, set logging to ERROR level [default INFO]") {Bio::Log::CLI.trace('error')}
    opts.on("--logger filename",String,"Log to file [default #{options[:logger]}]") { |name| options[:logger] = name}
    opts.on("--trace options",String,"Set log level [default INFO]. e.g. '--trace debug' to set logging level to DEBUG"){|s| Bio::Log::CLI.trace(s)}
  end
  o.parse!
  if ARGV.length != 0 or options[:dot_file].nil?
    $stderr.puts o
    exit 1
  end
  # Setup logging. bio-logger defaults to STDERR not STDOUT, I disagree
  Bio::Log::CLI.logger(options[:logger]); log = Bio::Log::LoggerPlus.new(LOG_NAME); Bio::Log::CLI.configure(LOG_NAME)
  
  log.warn "Note that this script currently cannot handle fully circular inputs (i.e. cyclic graphs). To handle these, break it up and then modify the output file manually"

  assembly = IonPairer::GraphvizReader.new.read_assembly(options[:dot_file])
  
  Dir.chdir(options[:output_directory]) do
    assembly.each_with_index do |scaffold, scaffold_counter|   
      File.open("scaffold#{scaffold_counter}.yml",'w') do |f|
        f.print scaffold.scaffolder_yaml
      end
    end
  end
  
  
end #end if running as a script