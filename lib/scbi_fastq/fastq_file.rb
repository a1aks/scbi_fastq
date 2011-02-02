
# add ord method to ruby 1.8
if !String.instance_methods.include?(:ord)
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
  def initialize(fasta_file_name, mode='r', fastq_type = :sanger, qual_to_array=true, qual_to_phred=true)

  	
    if mode.upcase.index('W')
      @fastq_file = File.open(fasta_file_name,'w')
    elsif mode.upcase.index('A')
      if !File.exist?(fasta_file_name)
      	raise "File #{fasta_file_name} doesn't exists" 
      end
    	
      @fastq_file = File.open(fasta_file_name,'a')
    else #read only
      if !File.exist?(fasta_file_name)
      	raise "File #{fasta_file_name} doesn't exists" 
      end
    	
      if fasta_file_name.is_a?(IO)
    	  @fastq_file = fasta_file_name
      else
        @fastq_file = File.open(fasta_file_name,'r')
      end
    end
    
    @mode = mode
    @num_seqs = 0
    @fastq_type=fastq_type
    
    #  S - Sanger        Phred+33,  raw reads typically (0, 40)
    #  X - Solexa        Solexa+64, raw reads typically (-5, 40)
    #  I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)
    #  J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)
    # > >>> def solexa_quality_from_phred(phred_quality) :
    # > ...     return 10*log(10**(phred_quality/10.0) - 1, 10)
    # > ...
    # > >>> solexa_quality_from_phred(90)
    # > 89.999999995657035
    # > >>> solexa_quality_from_phred(50)
    # > 49.99995657033466
    # > >>> solexa_quality_from_phred(10)
    # > 9.5424250943932485
    # > >>> solexa_quality_from_phred(1)
    # > -5.8682532438011537
    # > >>> solexa_quality_from_phred(0.1)
    # > -16.32774717238372
    # > 
    # > >>> def phred_quality_from_solexa(solexa_quality) :
    # > ...     return 10*log(10**(solexa_quality/10.0) + 1, 10)
    # > ...
    # > >>> phred_quality_from_solexa(90)
    # > 90.000000004342922
    # > >>> phred_quality_from_solexa(10)
    # > 10.41392685158225
    # > >>> phred_quality_from_solexa(0)
    # > 3.0102999566398116
    # > >>> phred_quality_from_solexa(-20)
    # > 0.043213737826425784
    
    
    #sanger by default
    @to_phred = lambda{|q| q - 33}
    @from_phred = lambda{|q| (q+33).chr}
    
    if @fastq_type == :ilumina
        @to_phred = lambda{|q| q - 64}
        # @from_phred = lambda{|q| (q+64).chr}
        
    elsif @fastq_type == :solexa
       # 
       # solexa to phred quals
       
       @to_phred = lambda{|q| (10*Math.log(10**(q/10.0)+1,10)).round}
       # @from_phred = lambda{|q| (10*Math.log(10**(q/10.0)-1,10)).round.chr}
       
       #phred to solexa quals
       
    end
    
    @qual_to_array = qual_to_array
    
    @qual_to_phred = qual_to_phred
    
  end
  
  def close
      @fastq_file.close
  end
 
  
  #------------------------------------
  # Scans all sequences
  #------------------------------------
  def each
        
    rewind

	  n,f,q,c=next_seq
	  
    while (!n.nil?)
			yield(n,f,q,c)
			n,f,q,c=next_seq
    end

  	rewind
  	
  end

  def rewind
     
     @num_seqs = 0 ;
     @fastq_file.pos=0
    
  end

  #------------------------------------
  # Scans a file, firing events to process content
  #------------------------------------
  def next_seq
    #init variables
    res = read_fastq
  	return res
  end
  
  
  def write_seq(seq_name,seq_fasta,seq_qual,comments='')
	  name = ""
	  
		@fastq_file.puts("@#{seq_name} #{comments}")
		@fastq_file.puts(seq_fasta)
		@fastq_file.puts("+#{seq_name} #{comments}")
		
		if seq_qual.is_a?(Array)
		  @fastq_file.puts(seq_qual.map{|e| @from_phred.call(e)}.join)
	  else
		  @fastq_file.puts(seq_qual.split(/\s+/).map{|e| @from_phred.call(e.to_i)}.join)
		end
    
  end
  
  def with_qual?
    true
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

    seq_name = nil
    seq_fasta = nil
    seq_qual = nil
    comments = nil
    
    reading = :fasta
    
    if !@fastq_file.eof
      
      begin
        #read four lines
        name_line = @fastq_file.readline.chomp
        seq_fasta = @fastq_file.readline.chomp
        name2_line = @fastq_file.readline.chomp
        seq_qual = @fastq_file.readline.chomp

      
        # parse name
        if name_line =~ /^@\s*([^\s]+)\s*(.*)$/
          # remove comments
          seq_name = $1
          comments=$2
        else
          raise "Invalid sequence name in #{line}"
        end
      
        # parse fasta
        seq_fasta.strip! if !seq_fasta.empty?

        # parse qual_name
    
        if !seq_name.nil? && !seq_qual.empty?

           @num_seqs += 1
       
           if @qual_to_phred
             seq_qual=seq_qual.each_char.map{|e| (@to_phred.call(e.ord))}

             if !@qual_to_array
                 seq_qual=seq_qual.join(' ')
             end
           end
       
        end
      rescue EOFError
        raise "Bad format in FastQ file"
      end
    end
    
    return [seq_name,seq_fasta,seq_qual,comments]
  end

  
  def read_fastq_old3

    seq_name = nil
    seq_fasta = nil
    seq_qual = nil
    comments = nil
    
    reading = :fasta
    
    

    # mientras hay lineas en el fichero
    while (!@fastq_file.eof)
      revert_pos = @fastq_file.pos
      # lee una linea
      reading = :name
      
      line = @fastq_file.readline; line.chomp!

      # if reading == :seq_name
        
      if reading == :name
        # si la linea es una nueva secuencia comienza por @
        if ((line =~ /^@/))
          #get only name
        
          # si hay algo que devolver, romper bucle
          if !seq_fasta.nil?
            @fastq_file.pos=revert_pos
            break
          end
        
          # remove starting
          line.gsub!(/^@\s*/,'')
        
          line =~ /(^[^\s]+)\s*(.*)$/
          # remove comments
          seq_name = $1
          comments=$2
        
          seq_fasta=''
          seq_qual=''
        
          reading = :fasta
        elsif ((line =~ /^+/))
        
          reading = :qual
        end
        
      # elsif ((line =~ /^\+#{seq_name}/)) || ((line =~ /^\+\s*$/))
        
      elsif reading == :fasta
          #add line to fasta of seq
  
          line.strip! if !line.empty?
        
          seq_fasta+=line;
          seq_fasta.strip! if !seq_fasta.empty?
      elsif reading == :qual
  
          line.strip! if !line.empty?
        
          seq_qual+=line;
          seq_qual.strip! if !seq_qual.empty?
      end
    end
    
    if !seq_name.nil? && !seq_qual.empty?

       @num_seqs += 1
       
       if @qual_to_phred
         seq_qual=seq_qual.each_char.map{|e| (@to_phred.call(e.ord))}

         if !@qual_to_array
             seq_qual=seq_qual.join(' ')
         end
       end
       
    end
    
    return [seq_name,seq_fasta,seq_qual,comments]
  end

  def read_fastq_old

    seq_name = nil
    seq_fasta = nil
    seq_qual = nil
    comments = nil
    
    reading = :fasta

    # mientras hay lineas en el fichero
    while (!@fastq_file.eof)
      revert_pos = @fastq_file.pos
      # lee una linea
      line = @fastq_file.readline; line.chomp!

      # si la linea es una nueva secuencia
      if ((line =~ /^@/))
        #get only name
        
        # si hay algo que devolver, romper bucle
        if !seq_fasta.nil?
          @fastq_file.pos=revert_pos
          break
        end
        
        # remove starting
        line.gsub!(/^@\s*/,'')
        
        line =~ /(^[^\s]+)\s*(.*)$/
        # remove comments
        seq_name = $1
        comments=$2
        
        seq_fasta=''
        seq_qual=''
        
        reading = :fasta
        
        
      elsif ((line =~ /^\+#{seq_name}/)) || ((line =~ /^\+\s*$/))
        reading = :qual
        
      else # no es una linea de nombre, añadir al fasta o al qual
        
        line.strip! if !line.empty?
        
        if reading == :fasta
          #add line to fasta of seq
          seq_fasta+=line;
          seq_fasta.strip! if !seq_fasta.empty?
	      
        else # reading qual 
          seq_qual+=line;
          seq_qual.strip! if !seq_qual.empty?
        end  
      end
    end
    
    if !seq_name.nil? && !seq_qual.empty?

       @num_seqs += 1
       
       if @qual_to_phred
         seq_qual=seq_qual.each_char.map{|e| (@to_phred.call(e.ord))}

         if !@qual_to_array
             seq_qual=seq_qual.join(' ')
         end
       end
       
    end
    
    return [seq_name,seq_fasta,seq_qual,comments]
  end
  
  
end
