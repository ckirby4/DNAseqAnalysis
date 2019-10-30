#DNAseq Analysis
#1. check for records in a file (count_records)
#2. calculate the length of each DNA sequence (check_length)
#3. identify ORF in each DNA sequence (orf_identifier)
#4. identify repeated motif in sequence (repeats_identifier)
#Inorder to use, please assign the path and file name(Must be FASTA) to the class function indicated below 

class dna_tool_sets ():
  
    def __init__(self, file_name):
        
        self.file_name = file_name
        self.dict = {} 
        f_reader = open (self.file_name)
        for line in f_reader:
            line = line.strip("\n") 
            if ">" in line:
                header = line 
                self.dict[header] = "" 
            else:
                self.dict[header] += line 
        f_reader.close()
    
    def count_records (self):
        """ this function count number of records in the file"""
        number_of_records = len(self.dict) 
        print "How many records are in the multi-FASTA file: %d \n"\
               %number_of_records
        
        
    def check_length(self):
        """
        this function will check the  length of each record and the according header 
        of the record
        """
        length_dict = {} 
        for key, value in self.dict.items():
            length_dict[key] = len(value)
            
        lengths = length_dict.values() 
        
        max_length = max(lengths)
        min_length = min(lengths)
        
        record_max_length = [item for item in length_dict if length_dict[item] == max_length]
               
        record_min_length = [item for item in length_dict if length_dict[item] == min_length]
        
        # length of the longest sequence in the file
        print "The length of the longest sequence: %d \n"%max_length, \
              "The number of longest sequence: %d \n"%len(record_max_length)
        # length of the shortest sequence in the file
        print "The length of the shortest sequence: %d \n"%min_length, \
              "The number of shortest sequence: %d \n"%len(record_min_length)
        
    
               
    def find_pos(self, dna):
        """ 
        help function for orf_identifier: 
        find start position and length for each read frame for a given sequence
        dna: sequence, string
        """
        start_code = "ATG"
        stop_codes = ["TAA", "TAG", "TGA"]

        pos_dict = {} 
        
        for i in range(3): # 3 different frames
            
            pos = []
            # generate frames
            if i == 0:
                frame = [dna[j:j+3] for j in range(i, len(dna), 3)]
            else:
                frame = [dna[:i]] + [dna[j:j+3] for j in range(i, len(dna), 3)]
            # all possible start postions and stop positions in seq
            start_pos = []
            stop_pos = []
            try:
                index_start_pos = [m for m, y in enumerate(frame) if \
                                  y == start_code]
                start_pos += index_start_pos # possible positions for start code "ATG"
            except ValueError:
                pos.append((-1, 0)) # if no ATG return -1 as start position
                continue
 
            for stop_code in stop_codes:
                try:
                    
                    index_stop_code = [n for n, x in enumerate(frame) if \
                                       x == stop_code and n > min(start_pos)]
                    stop_pos += index_stop_code
                except ValueError:
                    continue
            if len(stop_pos) == 0: # if no stop code find add -1 
                 pos.append((-1, 0))
            else:
                #closest paired code
                 while len(start_pos) != 0:
                     start = min(start_pos)
                     try:
                         end = min([stop for stop in stop_pos if stop > start])
                     except ValueError:
                         break     
                 # start position and length
                     s_pos = len("".join(frame[:start])) + 1
                     pos.append((s_pos, (end - start + 1)*3))
                     start_pos.remove(start) 
            pos_dict["frame%d"%(i+1)] = pos 
            
        return pos_dict
        
    def revs_complement(dna):
        """
        help function for orf_identifier:
        to transform a sequence to reverse complementary sequence
        """
        pairs = {"A": "T", "C": "G", "G": "C", "T": "A"} 
        c_dna = [pairs[s] for s in dna] 
        return "".join(c_dna)[::-1] 
        
    
    def orf_identifier (self):
        """
        This function return all the orf information with start posiotion and
        length of orf in read frame 1, 2 and 3
        the values for frame 1, 2, and are represented as pairs of tuple in a list
        for example: {"header1": {"frame1":[(0, 100)], "frame2":[(20, 400)], "frame3":[(-1, 0)]},...} 
        represents for header1:
        frame1- start position: 0, length of orf: 100
        frame2- start position: 20, length of orf: 400
        frame3- No ORF detected
        """
        orf = {}
        for header, dna_seq in self.dict.items(): # generate orf for the whole file
            pos = self.find_pos(dna_seq)
            orf[header] = pos
        # find header for question 7
        id_key = [key for key in orf if "gi|142022655|gb|EQ086233.1|16" in key]
        idx = id_key[0]  
        # generate list of frames for questions 4 to 7
        frame1, frame2, frame3, all_frames, id_frames = [], [], [], [], []
        for key, dict_value in orf.items():
            frame1 += dict_value["frame1"]
            frame2 += dict_value["frame2"]
            frame3 += dict_value["frame3"]
            frames = dict_value["frame1"] + dict_value["frame2"] + dict_value["frame3"]
            all_frames += frames
            if key == idx:
                id_frames = dict_value["frame1"] + dict_value["frame2"] + dict_value["frame3"]
            
        
        #frame 2 of any of the sequences?
        frame2_max_length = max(frame2, key = lambda x: x[1])
        print "The length of longest ORF in frame2: %d\n"%frame2_max_length[1]
        
       
        
        frame1_max_length_pos = max(frame3, key = lambda x: x[1])
        print "The start position of longest ORF in frame3: %d\n"%frame1_max_length_pos[0]
        
        max_length = max(all_frames, key = lambda x: x[1])
        print "The longest ORF of all frames and sequences: %d\n"%max_length[1]
        

        max_length_id = max(id_frames, key = lambda x: x[1])
        print "The length of longest ORF for ", idx, "is: %d \n" %max_length_id[1]
        
        
            
    def find_repeats(self, dna, n):
        """
        This help function for repeats_identifier find and count repeats for 
        each dna sequence
        dna: sequence, string
        n: number of repeats, int
        """
        repeats = {}
        for i in range(0, len(dna)):
            repeat = dna[i:i+n] # generate possible repeats
            if len(repeat) == n:
                if repeat not in repeats:
                    repeats [repeat] = 1 # initiate record
                else:
                    # count repeated repeats
                    repeats[repeat] = repeats.get(repeat) + 1
        return repeats
    
    def repeats_identifier(self, n):
        """
        This function generates repeats with counts for each record 
        (repeats_set) and for the whole file (combined_repeats)
        n: number of repeats, int
        """
        # record the repeats with counts for each record 
        repeats_set = {}
        for header, dna_seq in self.dict.items():
            repeats = self.find_repeats(dna_seq, n)
            repeats_set[header] = repeats 
        # record the repeats with counts for the whole file
        combined_repeats = {}
        for dict_value in repeats_set.values():
            for key in dict_value:
                if key not in combined_repeats:
                    combined_repeats[key] = dict_value[key]
                else:
                    combined_repeats[key] = combined_repeats.get(key) \
                                            + dict_value[key]
        
        if n == 6:
            most_freq_7 = max (combined_repeats.values())
            print "The most frequently repeats occur: %d times \n"%most_freq_7
       
        if n == 7:
            most_freq_7_seq = [key for key in combined_repeats if \
                       combined_repeats[key] == max(combined_repeats.values())]
            print "The following repeats occured most frequently: \n", most_freq_7_seq
        
        if n == 12:
            
            count_most_freq_10 = len([value for value in combined_repeats.values()\
                             if value == max(combined_repeats.values())])
        
        
            print "The number of different 10-base sequences occur max times: %d \n"\
                  %count_most_freq_10
        
        # uncomment the return code to return the repeats with counts for each record 
        # (repeats_set) and for the whole file (combined_repeats)
        #return repeats_set, combined_repeats
        
if __name__ == "__main__":
    
    file_name = 'Example.fasta'    #INSERT FILE NAME HERE 
    dna_tools = dna_tool_sets (file_name)
    # How many records are in the multi-FASTA file?
    dna_tools.count_records()
    
    # What is the length of the longest sequence in the file?
    # What is the length of the shortest sequence in the file?
    dna_tools.check_length()
    
    # What is the length of the longest ORF appearing in reading
    # frame 2 of any of the sequences?
    # What is the starting position of the longest ORF in reading frame 1 
    # any of the sequences? 
    # What is the length of the longest ORF appearing in any sequence and 
    #in any forward reading frame?
   
    dna_tools.orf_identifier()
    
    # Find the most frequently occurring repeat of length 7 in all 
    # sequences. How many times does it occur in all?
    dna_tools.repeats_identifier(6)
    
    dna_tools.repeats_identifier(12)
    
    
    dna_tools.repeats_identifier(7)
    
    
