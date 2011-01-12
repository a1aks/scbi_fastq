require File.dirname(__FILE__) + '/test_helper.rb'

class TestScbiFastq < Test::Unit::TestCase

  def setup
   	@test_file='/tmp/sanger.fastq';
  	
  	@seq_fasta='ACTG'
		@seq_qual=[25]
	  @seq_name='SEQ'
	  
  end
  
  
  def fill_file(n,offset=33)
  	f=File.new(@test_file,'w')
  	
  	n.times do |c|
  	  i = c+1
  	  name = "@#{@seq_name+i.to_s} comments"
  	  
  		f.puts(name)
  		f.puts(@seq_fasta*i)
  		f.puts('+'+name)
  		f.puts((@seq_qual*i*@seq_fasta.length).map{|e| (e+offset).chr}.join)
  	end
  	
  	f.close
  end

  def test_each

    	 # make new file and fill with data
		  fill_file(100) 	
		

			fqr=FastqFile.new(@test_file)

			i=1
			
			fqr.each do |n,s,q|
			  
			  assert_equal(@seq_name+i.to_s,n)
			  assert_equal(@seq_fasta*i,s)
			  assert_equal((@seq_qual*i*@seq_fasta.length).join(' '),q)

				i+=1
			end
			 
		  fqr.close			
  end

  def test_each_comments

    	 # make new file and fill with data
		  fill_file(100) 	
		

			fqr=FastqFile.new(@test_file)

			i=1
			
			fqr.each do |n,s,q,c|
			  
			  assert_equal(@seq_name+i.to_s,n)
			  assert_equal(@seq_fasta*i,s)
			  assert_equal((@seq_qual*i*@seq_fasta.length).join(' '),q)
			  assert_equal('comments',c)

				i+=1
			end
			 
		  fqr.close			
  end

  
  # def test_open_file
  #   fill_file(100)
  #   fq=FastqFile.new('test/sanger.fastq')
  #   
  #   fq.each do |n,f,q|
  #      puts n,f,q
  #      puts fq.num_seqs
  #   end
  #   
  #   fq.close
  #   
  # end
end
