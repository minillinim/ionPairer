require 'rspec'
require 'pp'
require 'open3'
require 'bio-logger'
require 'tmpdir'
require 'graphviz'

# To run this test:
# $ rspec /path/to/test_script_being_tested.rb

# Assumes that the name of the file being tested is ../something.rb, and the name of this script is test_something.rb
$:.unshift File.join(File.dirname(__FILE__),'..')
script_under_test = File.basename(__FILE__).gsub(/^test_/,'')
require script_under_test
def assert_equal(e,o); o.should eq(e); end
path_to_script = File.join(File.dirname(__FILE__),'..',script_under_test)

# setup logging
log_name = File.basename(script_under_test,'.rb')
Bio::Log::LoggerPlus.new(log_name)
Bio::Log::CLI.logger('stderr')
#Bio::Log::CLI.trace('warn')
#Bio::Log::CLI.configure(log_name) # when commented out no debug is printed out

describe 'ContigLink' do
  it 'should predict_link_direction on longer contigs' do
    # long contigs
    #predict_link_direction(max_distance, contig_length, position, direction)
    #START_OF_CONTIG = 'contig_start'
    #END_OF_CONTIG = 'contig end'
    #INCONSISTENT_POSITION_AND_DISTANCE = 'inconsistent position and distance on contig ends'
    #NOT_NEAR_EITHER_END = 'not near either end of the contig'
    ContigLink.predict_link_direction(10, 30, 5, '-').should eq(ContigLink::START_OF_CONTIG)
    ContigLink.predict_link_direction(2, 30, 5, '-').should eq(ContigLink::NOT_NEAR_EITHER_END)
    ContigLink.predict_link_direction(10, 30, 1, '+').should eq(ContigLink::INCONSISTENT_POSITION_AND_DISTANCE)
    ContigLink.predict_link_direction(10, 30, 25, '+').should eq(ContigLink::END_OF_CONTIG)
    ContigLink.predict_link_direction(10, 30, 25, '-').should eq(ContigLink::INCONSISTENT_POSITION_AND_DISTANCE)
    ContigLink.predict_link_direction(10, 30, 15, '-').should eq(ContigLink::NOT_NEAR_EITHER_END)
  end
  
  it 'should predict_link_direction on longer contigs' do
    # shorter contigs
    #predict_link_direction(max_distance, contig_length, position, direction)
    ContigLink.predict_link_direction(100, 30, 5, '-').should eq(ContigLink::START_OF_CONTIG)
    ContigLink.predict_link_direction(100, 30, 1, '+').should eq(ContigLink::END_OF_CONTIG)
  end
  
  it 'ContigLink should work when there is multiple pairs' do
    #    attr_accessor :position1, :position1,
    #:direction1, :direction2
    link = ContigLink.new
    link.position1 = 10
    link.position2 = 400
    link.direction1 = '-'
    link.direction2 = '+'
    link.predict_link_ends(200, 100, 500).should eq([ContigLink::START_OF_CONTIG, ContigLink::END_OF_CONTIG])
    link.predict_link_ends(200, 100, 5000).should eq([ContigLink::START_OF_CONTIG, ContigLink::NOT_NEAR_EITHER_END])
  end
  
  it 'should estimate insert size reasonably' do
    link = ContigLink.new
    link.position1 = 10
    link.position2 = 400
    link.direction1 = '-'
    link.direction2 = '+'
    #def estimate_distance_between_contigs(contig_lengths, contig_ends, mean_insert_size)
    link.estimate_distance_between_contigs([500,600], [ContigLink::START_OF_CONTIG, ContigLink::START_OF_CONTIG], 2000).should eq(2000-400-10)
    link.estimate_distance_between_contigs([500,600], [ContigLink::START_OF_CONTIG, ContigLink::END_OF_CONTIG], 2000).should eq(2000-10-(600-400))
  end
end

describe 'ContigLinkageSet' do
  it 'should classify right' do
    set = ContigLinkageSet.new
    set.contig1_length = 30
    set.contig2_length = 100

    link1 = ContigLink.new
    link1.position1 = 30
    link1.position2 = 80
    link1.direction1 = '-'
    link1.direction2 = '+'
    link2 = ContigLink.new
    link2.position1 = 20
    link2.position2 = 80
    link2.direction1 = '-'
    link2.direction2 = '+'
    set.links = [link1, link2]
    set.classified_abuttings_hash(50).should eq(
      [ContigLink::START_OF_CONTIG, ContigLink::END_OF_CONTIG] =>
      [link1, link2]
    )
  end
  
  it 'should predict_best_abutting' do
    set = ContigLinkageSet.new
    set.contig1_length = 30
    set.contig2_length = 100

    link1 = ContigLink.new
    link1.position1 = 30
    link1.position2 = 80
    link1.direction1 = '-'
    link1.direction2 = '+'
    link2 = ContigLink.new
    link2.position1 = 20
    link2.position2 = 80
    link2.direction1 = '-'
    link2.direction2 = '+'
    link3 = ContigLink.new
    link3.position1 = 20
    link3.position2 = 80
    link3.direction1 = '+' #this is what is different to link2
    link3.direction2 = '+'
    set.links = [link1, link2, link3]
    abutting1, abutting2, abutting_links = set.predict_best_abutting(50, 1)
    abutting1.should eq(ContigLink::START_OF_CONTIG)
    abutting2.should eq(ContigLink::END_OF_CONTIG)
    abutting_links.should eq([link1, link2])
    
    abutting1, abutting2, abutting_links = set.predict_best_abutting(50, 3)
    abutting1.should eq(nil)
  end
  
  it 'predict_best_abutting should fail when there there is an equal number of links in 2 directions' do
    set = ContigLinkageSet.new
    set.contig1_length = 30
    set.contig2_length = 100

    link1 = ContigLink.new
    link1.position1 = 30
    link1.position2 = 80
    link1.direction1 = '-'
    link1.direction2 = '+'
    link2 = ContigLink.new
    link3 = ContigLink.new
    link3.position1 = 20
    link3.position2 = 80
    link3.direction1 = '+' #this is what is different to link2
    link3.direction2 = '+'
    set.links = [link1, link3]
    abutting1, abutting2, abutting_links = set.predict_best_abutting(50, 1)
    abutting1.should eq(nil)
  end
  
  it 'should predict correct insert size' do
    # self.distance_between_contigs(contig_lengths, abuttings, abutting_links, mean_insert_size)
    #set.contig1_length = 30
    #set.contig2_length = 100
    #
    link1 = ContigLink.new
    link1.position1 = 30 #=> 30
    link1.position2 = 80 #=> 20
    link1.direction1 = '-'
    link1.direction2 = '+'
    link2 = ContigLink.new
    link2.position1 = 20 #=> 20
    link2.position2 = 80 #=> 20
    link2.direction1 = '-'
    link2.direction2 = '+'
    #set.links = [link1, link2]
    #set.classified_abuttings_hash(50).should eq(
    #  [ContigLink::START_OF_CONTIG, ContigLink::END_OF_CONTIG] =>
    #  [link1, link2]
    #)
    ContigLinkageSet.distance_between_contigs(
      [30, 100],
      [ContigLink::START_OF_CONTIG, ContigLink::END_OF_CONTIG],
      [link1, link2],
      1000
    ).should eq(1000-(50+40)/2)
    
    ContigLinkageSet.distance_between_contigs(
      [30, 100],
      [ContigLink::START_OF_CONTIG, ContigLink::END_OF_CONTIG],
      [link1, link2],
      100
    ).should eq(100-(50+40)/2)
  end
end

describe 'ContigLinkSet' do
  it 'should add contig links ok' do
    #add_contig_set(contig1_name, pair1_position, pair1_direction,
    #                   contig2_name, pair2_position, pair2_direction)
    set = ContigLinkSet.new
    set.length.should eq(0)
    set.add_contig_set('contig2', 200, '-', 'contig1', 100, '+')
    set.values.length.should eq(1)
    set.keys[0].should eq(['contig1', 'contig2'])
    set.values[0].length.should eq(1)
    set.values[0][0].position1.should eq(100)
  end
  
  it 'should graphviz on a single link. Easy right?' do
    set = ContigLinkSet.new
    set.length.should eq(0)
    set.add_contig_set('contig2', 2500, '-', 'contig1', 100, '+')
    set.add_contig_set('contig2', 2501, '-', 'contig1', 10, '+')
    #set.add_contig_set('contig2', 2501, '+', 'contig1', 10, '+')
    set.contig_lengths = {'contig1' => 1000, 'contig2' => 3000}
    graph = set.generate_graphviz(54000, :min_links => 1, :mean_insert_size => 2000)
    graph.node_count.should eq(4)
    graph.edge_count.should eq(3)
    expecteds = {
      'contig1START' => 'contig1END',
      'contig2START' => 'contig2END',
      'contig1END' => 'contig2START'
    }
    edge_names = []
    graph.each_edge do|edge|
      edge_names.push [edge.node_one, edge.node_two]
      expecteds[edge.node_one].should eq(edge.node_two), "from #{edge.node_one}, expected #{expecteds[edge.node_one]}"
    end
    edge_names.uniq.length.should eq(3)
  end
  
  it 'should graphviz with 2 sets of contigs' do
    set = ContigLinkSet.new
    set.length.should eq(0)
    set.add_contig_set('contig2', 2500, '-', 'contig1', 100, '+')
    set.add_contig_set('contig2', 2501, '-', 'contig1', 10, '+')
    set.add_contig_set('contig3', 288, '+', 'contig4', 10, '-')
    set.contig_lengths = {'contig1' => 1000, 'contig2' => 3000, 'contig3' => 300, 'contig4' => 400}
    graph = set.generate_graphviz(54000, :min_links => 1, :mean_insert_size => 2000)
    
    graph.node_count.should eq(8)
    graph.edge_count.should eq(6)
    expecteds = {
      'contig1START' => 'contig1END',
      'contig2START' => 'contig2END',
      'contig1END' => 'contig2START',
      'contig3START' => 'contig3END',
      'contig4START' => 'contig4END',
      'contig3END' => 'contig4START',
    }
    edge_names = []
    graph.each_edge do|edge|
      edge_names.push [edge.node_one, edge.node_two]
      expecteds[edge.node_one].should eq(edge.node_two), "from #{edge.node_one}, expected #{expecteds[edge.node_one]}"
    end
    edge_names.uniq.length.should eq(6)
  end
  
  it 'should gracefully handle no links' do
    set = ContigLinkSet.new
    set.length.should eq(0)
    set.add_contig_set('contig2', 2500, '-', 'contig1', 100, '+')
    set.contig_lengths = {'contig1' => 1000, 'contig2' => 3000}
    graph = set.generate_graphviz(10, :min_links => 1, :mean_insert_size => 2000) #set the max distance to something impossible for testing purposes
    
    graph.node_count.should eq(4)
    graph.edge_count.should eq(2)
    expecteds = {
      'contig1START' => 'contig1END',
      'contig2START' => 'contig2END',
    }
    edge_names = []
    graph.each_edge do|edge|
      edge_names.push [edge.node_one, edge.node_two]
      expecteds[edge.node_one].should eq(edge.node_two), "from #{edge.node_one}, expected #{expecteds[edge.node_one]}"
    end
    edge_names.uniq.length.should eq(2)
  end
  
  it 'should respect min_links' do
    set = ContigLinkSet.new
    set.length.should eq(0)
    set.add_contig_set('contig2', 2500, '-', 'contig1', 100, '+')
    set.contig_lengths = {'contig1' => 1000, 'contig2' => 3000}
    graph = set.generate_graphviz(54000, :min_links => 1, :mean_insert_size => 2000).edge_count.should eq(3)
    graph = set.generate_graphviz(54000, :min_links => 2, :mean_insert_size => 2000).edge_count.should eq(2)
    graph = set.generate_graphviz(54000, :min_links => 3, :mean_insert_size => 2000).edge_count.should eq(2)
  end
end




describe 'the script' do
  
  it 'should work somewhat' do
    #  Dir.chdir tmpdirectory #run everything in a temporary directory

    # Create a dummy input file
    input = [
      'h1 h2 h3 h4 h5 h6 h7 h8',
      "contig1	25937	149	0	contig2	23	17	1",
      "contig1	25941	20	0	contig2	54	142	1",
      "contig3	8990	9	1	contig4	7389	153	0",
      "contig5	14033	20	1	contig6	12003	24	0",
      ].collect{|s| s.split(/\s/).join("\t")}
    fasta = [
      '>contig1',
      'A'*26000,
      '>contig2',
      'A'*29500,
      '>contig3',
      'A'*100000,
      '>contig4',
      'A'*100000,
      '>contig5',
      'A'*100000,
      '>contig6',
      'A'*100000,
    ]
    
    
    err = nil
    links_basename = 'contig_linker_test.links'
    fasta_filename = 'contig_linker_test.fasta'
    File.open(links_basename,'w') do |links_file|
      File.open(fasta_filename,'w') do |fasta_file|
        links_file.puts input.join("\n")
        links_file.close
        
        fasta_file.puts fasta.join("\n")
        fasta_file.close
        
        command = [
          path_to_script,
          '-l',
          links_file.path,
          '-f',
          fasta_file.path,
          '-L',
          '0',
          '-q',
          '-u',
          '2000',
        ]
        res = `#{command.join(' ')} 2>error`
        res.should eq("")
        err = File.open('error').readlines
        err.length.should eq(0), "error was #{err}"
      end
    end
    
    dot_filename = "#{links_basename}.dot"
    png_filename = "#{links_basename}.png"
    svg_filename = "#{links_basename}.svg"
    File.exists?(dot_filename).should be true
    File.exists?(png_filename).should be true
    File.exists?(svg_filename).should be true
    
    # Read in the graphviz file, make sure it is all good
    GraphViz.parse(dot_filename) do |g|
      g.node_count.should eq(12)
      g.edge_count.should eq(6+1), "graph was #{g}, error output was #{err.join("")}"
    end
    
    #clean up
    [
      links_basename,
      fasta_filename,
      'error',
      dot_filename,
      png_filename,
      svg_filename,
    ].each do |f|
      File.delete f
    end
  end
end
