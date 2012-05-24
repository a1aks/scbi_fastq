$:.unshift(File.dirname(__FILE__)) unless
  $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require 'scbi_fastq/fastq_file'
module ScbiFastq
   VERSION = '0.0.15'
end
