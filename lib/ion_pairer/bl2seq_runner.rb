require 'tempfile'
require 'bio'

# Given 2 Bio:seq objects, blast them against each other and return the result.
module Bio
  class BlastPlus
    class Hsp
      attr_accessor :qseqid, :sseqid, :pident, :length, :mismatch, :gapopen, :qstart, :qend, :sstart, :subject_end, :evalue, :bitscore
    end
    
    class Runner
      # Run a Bio::Seq object against another. Assumes bl2seq is working correctly
      def bl2seq(seq1, seq2)
        Tempfile.open('rubybl2seq') { |t1|  
          t1.puts '>one'
          t1.puts seq1
          t1.close
          
          Tempfile.open('rubybl2seq') { |t2|  
            t2.puts '>two'
            t2.puts seq2
            t2.close

            # Run the bl2seq. Assume nucleotide blast for the moment
            cmd = [
              'blastn',
              '-task',
              'blastn',
              '-query',
              t1.path,
              '-subject',
              t2.path,
              '-outfmt',
              '6',
            ]
            
            hsps = []
            Bio::Command.call_command_open3(cmd) do |stdin, stdout, stderr|
              # Create the report from the output file
              err = stderr.read
              
              raise err unless err.nil? or err == ''
              str = stdout.read
              
              str.each_line do |line|
                hsp = Bio::BlastPlus::Hsp.new
                
                splits = line.strip.split("\t")
                fields = [
                  :qseqid=, :sseqid=, :pident=, :length=, :mismatch=, :gapopen=, :qstart=, :qend=, :sstart=, :subject_end=, :evalue=, :bitscore=
                ]
                raise unless splits.length == fields.length
                 
                fields.each_with_index do |sym, i|
                  hsp.send sym, splits[i]
                end
                hsps.push hsp
              end
            end
            return hsps
          }
        }
      end
    end
  end
end