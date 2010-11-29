
# add ord method to ruby 1.8
if !String.instance_methods.include?('ord')
   class String 
      
     def ord
       return self[0]
     end
     
   end
end
      


class FastqFile

  attr_accessor :num_seqs

  #------------------------------------
  # Initialize instance
  #------------------------------------
  def initialize(fasta_file_name, fastq_type = :sanger, qual_to_array=false, qual_to_phred=true)

  	
  	if !File.exist?(fasta_file_name)
    	raise "File #{fasta_file_name} doesn't exists" 
  	end
    @file_fasta = File.open(fasta_file_name)

    @num_seqs = 0
    @fastq_type=fastq_type
    @qual_to_array = qual_to_array
    
    @qual_to_phred = qual_to_phred
    
  end
  
  def close
      @file_fasta.close
  end
 
  
  #------------------------------------
  # Scans all sequences
  #------------------------------------
  def each
        
    rewind

	  n,f,q=next_seq
	  
	  while (!n.nil?)
			yield(n,f,q)
			n,f,q=next_seq
		end

  	rewind
  	
  end

  def rewind
     
     @num_seqs = 0 ;      
     @file_fasta.pos=0
    
  end

  #------------------------------------
  # Scans a file, firing events to process content
  #------------------------------------
  def next_seq
    #init variables
    res = read_fastq
  	return res
  end
  
  
  private 
  
  #------------------------------------
  #  Read one sequence fasta
  #------------------------------------
  # @GEM-108-D02
  # AAAAGCTGG
  # +
  # :::::::::
  
  def read_fastq

    seq_name=nil
    seq_fasta=nil
    seq_qual=nil
    
    reading = :fasta

    # mientras hay lineas en el fichero
    while (!@file_fasta.eof)
      revert_pos = @file_fasta.pos
      # lee una linea
      line_fasta = @file_fasta.readline; line_fasta.chomp!

      # si la linea es una nueva secuencia
      if ((line_fasta =~ /^@/))
        #get only name
        
        # si hay algo que devolver, romper bucle
        if !seq_fasta.nil?
          @file_fasta.pos=revert_pos
          break
        end
        
        # remove starting
        line_fasta.gsub!(/^@\s*/,'')
        
        line_fasta =~ /(^[^\s]+)\s*(.*)$/
        # remove comments
        seq_name = $1
        comments=$2
        
        seq_fasta=''
        seq_qual=''
        
        reading = :fasta
        
        
      elsif ((line_fasta =~ /^\+/))
        reading = :qual
        
      else # no es una linea de nombre, añadir al fasta o al qual
        
        line_fasta.strip! if !line_fasta.empty?
        
        if reading == :fasta
          #add line to fasta of seq
          seq_fasta+=line_fasta;
          seq_fasta.strip! if !seq_fasta.empty?
	      
        else # reading qual 
          seq_qual+=line_fasta;
          seq_qual.strip! if !seq_qual.empty?
        end  
      end
    end
    
    if !seq_name.nil? && !seq_qual.empty?

       @num_seqs += 1
       
       if @qual_to_phred
         seq_qual=seq_qual.each_char.map{|e| (e.ord-33)}

         if !@qual_to_array
             seq_qual=seq_qual.join(' ')
         end
       end
       
    end
    
    return [seq_name,seq_fasta,seq_qual]
  end

  
end
