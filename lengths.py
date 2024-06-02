# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 18:43:45 2024

@author: Sebastian
"""

import os, re, csv, time, multiprocessing, subprocess

def main():
    print()
    
    print("Number of cpu : ", multiprocessing.cpu_count())
    performance_param = True
    
    
    
    #The following lines must NOT be commented out
    lengths = length_class(performance_param)
    lengths.parse_params()
    lengths.set_outputs()
    lengths.check_outputs()
    lengths.set_executables()
    lengths.check_executables()
    
    lengths.set_samplelist()
    print("Get primers")
    lengths.set_primers()
    
    
    #All parts of the pipeline, individual lines can be commented out depending what parts of the pipeline you want to run.
    start = time.time()
    #errormessage, number, number_filled, number_empty, total_size = lengths.runTrimomatic()
    #errormessage, number, number_filled, number_empty, total_size = lengths.runUsearchMergePairs()
    print("Start")
    #number, number_filled, number_empty, total_size, primer_number, surviving_sequences, abandondend_sequences = lengths.runDemultiplex()
    #print(total_size)
    #print(surviving_sequences)
    #print(abandondend_sequences)
    #lowest_length, filename_lowest, largest_length, filename_largest, number, number_filled, number_empty, total_size = lengths.getLengthStatistics()
    lengths.GenotypeLength()
    end = time.time()
    
    
    
    
    
    print(end - start)
    
    
    
    
    
    return 0


class length_class():
    def __init__(self, performance_param):
        self.workspace = os.getcwd()
        
        self.performance_param = performance_param

        

        self.params = {
            "OutputFolder" : os.path.normpath(self.workspace + "/Default"),
            "Bin" : "None",
            "RawData" : "None",
            "PrimerFile" : "None",
            "SampleFile" : "None",
            "Adapters" : "TrueSeqAdaptersInUsage.fa",
            "RExecutable" : "R",
            "AlleleList" : "Missing",
            "Operatingsystem" : "windows",
            "Ploidy" : "diploid",
            "Indexcomboposition" : 1,
            "Uniqueidentifier" : "Default",
            "MinCount" : 20,
            "MinLength" : 290,
            "MaxMismatch" : 2,
            "ConsensusThreshold" : 0.7,
            "LengthWindow" : ['310', '600'],
            "Executables" : {
                "Trimmdir" :  "None",
                "Adaptersdir" : "None",
                "trimmexec" : "None",
                "usearchexec" : "None",
                "scriptRDiploid" : "None",
                "scriptRHaploid" : "None"
            },
            "SamplesOldNew" : {},
            "Primers" : {
                "PrimerInformation" : {},
                "Boundaries" : {}
            },
            "Samples" : [],
            "DemultiplexParams" : {
                "Number_Surviving_Sequences" : 0,
                "Number_Abandonded_Sequences" : 0,
            },
            "LengthStatistics" : {
                "LowestLength" : 0,
                "FileNameLowest" : "",
                "LargestLength" : 0,
                "FileNameLargest" : ""
            },
            "Multiprocessing" : {
                "Cores" :  multiprocessing.cpu_count()
            }
        }
        
        self.list_error = []
        #self.primers_used = []
        
        self.lengths_counter = 0
        
        self.counter_boundaries = 1

        

    def parse_params(self):
        params = os.path.normpath(self.workspace + '/parameters.txt')
        if not(os.path.exists(params)):
            raise Exception("There is no parameter file with the name parameters.txt")
        
        with open(params, 'r') as paramfile:
            for line in paramfile:
                line = line.rstrip("\r\n")
                if (not line) or line.startswith("#"):
                    continue
                self.set_param(line)
        
        #print(self.params)
        
    def set_param(self, line):
        elements = line.split("=")
        key = elements[0].rstrip("\r\n").strip()
        value = elements[1].rstrip("\r\n").strip(".").strip() #.strip(".")
        if (value and not(value == ".")):
            self.params[key] = value
            

    def set_outputs(self):
        output = self.params["OutputFolder"]
        try:
            os.mkdir(output)
        except OSError:
            print ("Creation of the directory %s failed because it's already there." % output) # if it fail it produces this error message
        else:
            print ("Successfully created the directory %s " % output)
        
        for folder in ['QC', 'SeparatOut', 'MergedOut', 'MarkerStatistics', 'AlleleLenghtCounts', 'MarkerPlots']:
            try:
                os.mkdir(os.path.normpath(output + '/' + folder))
            except OSError:
                print ("Creation of the directory %s failed because it's already there." % folder) # if it fail it produces this error message
            else:
                print ("Successfully created the directory %s " % folder)
    
    def set_executables(self):
        binpath = self.params["Bin"]
        for subdir, dirs, files in os.walk(binpath):
            if (os.path.basename(subdir).split('-')[0] == 'Trimmomatic'):
                trimmdir = os.path.normpath(subdir)
                self.params["Executables"]["Trimmdir"] = os.path.normpath(trimmdir)
            elif (os.path.basename(subdir) == 'adapters'):
                adaptersdir = os.path.normpath(subdir)
                self.params["Executables"]["Adaptersdir"] = os.path.normpath(adaptersdir)
                  
            for file in files:
                if (file.split("-")[0] == "trimmomatic" and ".jar" in file):
                    self.params["Executables"]["trimmexec"] = os.path.normpath(trimmdir + '/' + file)
                elif file.startswith('usearch'):
                    if (self.params["Operatingsystem"] == "windows" and "win" in file):
                        self.params["Executables"]["usearchexec"] = os.path.normpath(binpath + '/' + file)
                    elif (self.params["Operatingsystem"] == "linux" and not("win" in file)):
                        self.params["Executables"]["usearchexec"] = os.path.normpath(binpath + '/' + file)
                elif file.startswith('Rscript_Markerlength_STUTTER_Color_BaryzentricMinsize20_notVerbose'):
                    self.params["Executables"]["scriptRDiploid"] = os.path.normpath(binpath + '/' + file)
                elif file.startswith('Rscript_Markerlength_Haploid'):
                    self.params["Executables"]["scriptRHaploid"] = os.path.normpath(binpath + '/' + file)
                
    
    def check_outputs(self):
        outputs = ["Bin", "RawData", "PrimerFile", "SampleFile"]
        for output in outputs:
            if (self.params[output] == "None" or not(os.path.exists(os.path.normpath(self.params[output])))):
                raise Exception(f"{output} is not properly set or does not exist!")
    
    def check_executables(self):
        executables = ["Trimmdir", "Adaptersdir", "trimmexec", "usearchexec", "scriptRDiploid", "scriptRHaploid"]
        for executable in executables:
            if (self.params["Executables"][executable] == "None" or not(os.path.exists(os.path.normpath(self.params["Executables"][executable])))):
                raise Exception(f"{executable} is not properly set or does not exist!")
                
    def return_params(self):
        return self.params["OutputFolder"], self.params["SampleFile"], self.params["PrimerFile"], self.params["RawData"]
    
    
    
    
    
    
    def set_samplelist(self):
        rows_new = []
        samplesheet = self.params["SampleFile"]
        
        with open(samplesheet) as file:
            for line in file:
                row = []
                self.params["Samples"].append(re.split(',|\t|;', line)[0])
                row.append(re.split(',|\t|;', line)[0])
                row.append(re.split(',|\t|;', line)[1])
                row.append(re.split(',|\t|;', line)[1].replace("_", "")) #.replace("-", "")
                rows_new.append(row)
        
        #print(rows_new)
        with open(samplesheet, 'w', newline='') as new_samplesheet:
            writer = csv.writer(new_samplesheet, delimiter =',', quoting=csv.QUOTE_MINIMAL)
            for row in rows_new:
                row = [element.replace("\n", "") for element in row]
                writer.writerow(row)
                
    
    def checkInputDir(self):
        FilesData = {}
        for subdir, dirs, files in os.walk(self.params["RawData"]): #instead of: for fastq in os.listdir(self.raw_data_dir):
            for fastq in files: 
                sample = fastq.split('_')[int(self.params["Indexcomboposition"])-1]
                if sample in self.params["Samples"]:
                    print(sample)
                    if (sample not in FilesData.keys()):
                        FilesData[sample] = []
                    FilesData[sample].append(os.path.normpath(subdir + "/" + fastq))
        
        return FilesData
    
    
    
    # def runTrimomatic(self):
    #     FileData = self.checkInputDir()
        
    #     binpath = self.params["Bin"]
    #     operatingsystem = self.params["Operatingsystem"]
    #     trimdir = os.path.basename(os.path.normpath(self.params["Executables"]["Trimmdir"]))
    #     trimexec = os.path.basename(os.path.normpath(self.params["Executables"]["trimmexec"]))
    #     adaptersdir = os.path.basename(os.path.normpath(self.params["Executables"]["Adaptersdir"]))
    #     adapterfile = self.params["Adapters"]
        
    #     list_error = []
    #     for sample, files in FileData.items():
    #         files.sort() # get R1 and R2 list and sort them to make sure that R1 always shows first
    #         R1 = files[0] # get path of R1 input file
    #         R2 = files[1]
            
    #         if (operatingsystem == "windows"):
    #             outFileR1 = R1.split('\\')[-1].split('.')[0]
    #             outFileR2 = R2.split('\\')[-1].split('.')[0]
    #         elif (operatingsystem == "linux"):
    #             outFileR1 = R1.split('/')[-1].split('.')[0]
    #             outFileR2 = R2.split('/')[-1].split('.')[0]
            
    #         R1paired = os.path.normpath(self.params["OutputFolder"] + "/QC/" + outFileR1 + "_QTpaired.fastq")
    #         R1unpaired = os.path.normpath(self.params["OutputFolder"] + "/QC/" + outFileR1 + "_QTunpaired.fastq")
    #         R2paired = os.path.normpath(self.params["OutputFolder"] + "/QC/" + outFileR2 + "_QTpaired.fastq")
    #         R2unpaired = os.path.normpath(self.params["OutputFolder"] + "/QC/" + outFileR2 + "_QTunpaired.fastq")
            
            
    #         adapter = "bin/" + trimdir + "/" + adaptersdir + "/" + adapterfile
    #         exe =  os.path.normpath(binpath + "/" + trimdir + "/" + trimexec)
            
    #         code = 'java -jar ' + exe + ' PE ' + R1 + ' ' + R2 + ' ' + R1paired + ' ' + R1unpaired + ' ' + R2paired + ' ' + R2unpaired + ' ILLUMINACLIP:' + adapter + ':2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20'
    #         print('running: ' + code) # print code
    #         os_code = os.system(code)
    #         list_error.append(os_code)
    #         if os_code != 0:
    #             print(f"The command exited with error code: {os_code}")
    #             #list_error.append(os_code)
    #         #else:
    #             #return os_code
        
    #     number, number_filled, number_empty = self.check_results(os.path.normpath(self.params["OutputFolder"] + "/QC/"))
        
    #     for error in list_error:
    #         if (error != 0):
    #             return error, number, number_filled, number_empty
            
        
    #     return 0, number, number_filled, number_empty
                
    
    
    
    
    def runTrimomatic(self):
        FileData = self.checkInputDir()
        self.list_error = []
        
        if (self.performance_param):
            with multiprocessing.Pool(processes = multiprocessing.cpu_count()-1) as pool:
                processes = []
                args_list = [(sample, files) for sample, files in FileData.items()]
                #chunk_size = max(1, len(args_list) // (multiprocessing.cpu_count() * 4))
                results = pool.starmap_async(self.run_trimmo_per_file, args_list) #chunksize=chunk_size
                processes = results.get()
        else:
            for sample, files in FileData.items():
                self.run_trimmo_per_file(sample, files)
        
        # if (self.performance_param):
        #     with multiprocessing.Pool(processes = multiprocessing.cpu_count()) as pool:
        #         processes = []
        #         args_list = [(sample, files) for sample, files in FileData.items()]
        #         for sample, files in FileData.items():
        #             result = pool.apply_async(self.run_trimmo_per_file, args = (sample, files))
        #             processes.append(result)
        #             print("hh")
        #         output = [result.get() for result in processes]
                
        # else:
        #     for sample, files in FileData.items():
        #         self.run_trimmo_per_file(sample, files)
        
        # if (self.performance_param):
        #     processes = []
        #     for sample, files in FileData.items():
        #         p = multiprocessing.Process(target=self.run_trimmo_per_file, args=(sample, files))
        #         processes.append(p)
        #         p.start()
        #     for p in processes:
        #         p.join()   
        # else:
        #     for sample, files in FileData.items():
        #         self.run_trimmo_per_file(sample, files)
        
        number, number_filled, number_empty, total_size = self.check_results(os.path.normpath(self.params["OutputFolder"] + "/QC/"))
        
        for error in self.list_error:
            if (error != 0):
                return error, number, number_filled, number_empty
            
        return 0, number, number_filled, number_empty, total_size
    
    
    def run_trimmo_per_file(self, sample, files):
        files.sort() # get R1 and R2 list and sort them to make sure that R1 always shows first
        R1 = files[0] # get path of R1 input file
        R2 = files[1]
        
        if (self.params["Operatingsystem"] == "windows"):
            outFileR1 = R1.split('\\')[-1].split('.')[0]
            outFileR2 = R2.split('\\')[-1].split('.')[0]
        elif (self.params["Operatingsystem"] == "linux"):
            outFileR1 = R1.split('/')[-1].split('.')[0]
            outFileR2 = R2.split('/')[-1].split('.')[0]
        
        R1paired = os.path.normpath(self.params["OutputFolder"] + "/QC/" + outFileR1 + "_QTpaired.fastq")
        R1unpaired = os.path.normpath(self.params["OutputFolder"] + "/QC/" + outFileR1 + "_QTunpaired.fastq")
        R2paired = os.path.normpath(self.params["OutputFolder"] + "/QC/" + outFileR2 + "_QTpaired.fastq")
        R2unpaired = os.path.normpath(self.params["OutputFolder"] + "/QC/" + outFileR2 + "_QTunpaired.fastq")
        
        adapter = "bin/" + os.path.basename(os.path.normpath(self.params["Executables"]["Trimmdir"])) + "/" + os.path.basename(os.path.normpath(self.params["Executables"]["Adaptersdir"])) + "/" + self.params["Adapters"]
        exe =  os.path.normpath(self.params["Bin"] + "/" + os.path.basename(os.path.normpath(self.params["Executables"]["Trimmdir"])) + "/" + os.path.basename(os.path.normpath(self.params["Executables"]["trimmexec"])))
        
        code = 'java -jar ' + exe + ' PE ' + R1 + ' ' + R2 + ' ' + R1paired + ' ' + R1unpaired + ' ' + R2paired + ' ' + R2unpaired + ' ILLUMINACLIP:' + adapter + ':2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20'
        print('running: ' + code) # print code
        if (self.params["Operatingsystem"] == "windows"):
            os_code = subprocess.run(code)
            #os_code = os.system(code)
            self.list_error.append(os_code)
            if os_code != 0:
                print(f"The command exited with error code: {os_code}")
        else:
            subprocess.call(code, shell=True)
            
            
    def runUsearchMergePairs(self):
        self.list_error = []
        input_directory = os.path.normpath(self.params["OutputFolder"] + '/QC') #The input for usearch are the trimmed sequences from the QC folder
        output_folder = os.path.normpath(self.params["OutputFolder"] + '/MergedOut')
        
        if (self.performance_param):
            with multiprocessing.Pool(processes = multiprocessing.cpu_count()-1) as pool:
                processes = []
                args_list = [(file, input_directory, output_folder) for file in os.listdir(input_directory)]
                results = pool.starmap_async(self.runUsearchMergePairs_per_file, args_list)
                processes = results.get()
        else:
            for file in os.listdir(input_directory):
                self.runUsearchMergePairs_per_file(file, input_directory, output_folder)
        
        # if (self.performance_param):
        #     processes = []
        #     for file in os.listdir(input_directory):
        #         p = multiprocessing.Process(target=self.runUsearchMergePairs_per_file, args=(file, input_directory, output_folder))
        #         processes.append(p)
        #         p.start()
        #     for p in processes:
        #         p.join()
        # else:
        #     for file in os.listdir(input_directory):
        #         self.runUsearchMergePairs_per_file(file, input_directory, output_folder)
        
        number, number_filled, number_empty, total_size = self.check_results(output_folder)
            
        for error in self.list_error:
            if (error != 0):
                return error, number, number_filled, number_empty, total_size
            
        return 0, number, number_filled, number_empty, total_size
    
    
    def runUsearchMergePairs_per_file(self, file, input_directory, output_folder):
        print(file)
        if file.endswith('QTpaired.fastq') and "R1" in file:
            R1 = os.path.normpath(input_directory + '/' + file) # get path of R1 input file
            R2 = R1.replace("_R1_", "_R2_") # get path of R2 input file
            out = os.path.normpath(output_folder + '/' + file.split('_')[int(self.params["Indexcomboposition"])-1] + '_joined.fastq') # get path of file to save merge reads
            exe = os.path.normpath(self.params["Bin"] + '/' + os.path.basename(os.path.normpath(self.params["Executables"]["usearchexec"]))) #self.usearch_exec_path
            code = exe + ' -fastq_mergepairs ' + R1 + ' -reverse ' + R2 + ' -fastqout ' + out + ' -fastq_pctid 80 -fastq_maxdiffs 40'
            print('running: ' + code) # print code
            print()
            if (self.params["Operatingsystem"] == "windows"):
                os_code = subprocess.run(code)
                #os_code = os.system(code)
                #subprocess.call(code, shell=True) # run program
                self.list_error.append(os_code)
                if os_code != 0:
                    print(f"The command exited with error code: {os_code}")
            else:
                subprocess.call(code, shell=True)
    
    
    
    def set_primers(self):
        primerFile = self.params["PrimerFile"]
        primer_dict={}
        primer_dict["PrimerInformation"] = {}
        primer_dict["Boundaries"] = {}
        with open(primerFile) as primer_file:
            for line in primer_file:
                line_list = re.split('\t|,|;', line.rstrip('\r\n')) # make list based seperating elemntes based on tabs
                if (len(line_list) > 3):
                    primers = line_list[1:3]
                    boundaries = line_list[3:]
                    print(line_list[3:])
                    boundaries = [(int(element.lstrip('\r\n').rstrip('\r\n').split()[0]), int(element.lstrip().rstrip().split()[1])) for element in line_list[3:]]
                    print(boundaries)
                    primer_dict["Boundaries"][line_list[0]] = boundaries
                else:
                    primers = line_list[1:] # primer sequences (contains forwards and reverse primer)
                primer_dict["PrimerInformation"][line_list[0]]=primers # dict: {locus_name: [Forward, Reverse]}
        
        self.params["Primers"]["PrimerInformation"] = {k : v for k, v in primer_dict["PrimerInformation"].items() if k != ""}
        self.params["Primers"]["Boundaries"] = primer_dict["Boundaries"]
        print(self.params["Primers"]["Boundaries"])
        #print(self.params["Primers"])
        return self.params["Primers"]["PrimerInformation"]
        
    
    
    def runDemultiplex(self): 
        input_directory = os.path.normpath(self.params["OutputFolder"] + '/MergedOut')
        demultiplexFile = self.params["PrimerFile"] #demultiplexing based on primers
        SampleFile = self.params["SampleFile"]
        outDir = os.path.normpath(self.params["OutputFolder"] + '/SeparatOut')
        
        
        #self.params["DemultiplexParams"]["Number_Surviving_Sequences"] = 0
        #self.params["DemultiplexParams"]["Number_Abandonded_Sequences"] = 0
        # save dictionary with primer information
        # primer_dict={}
        # with open(demultiplexFile) as primer_file:
        #     for line in primer_file:
        #         line_list = re.split('\t|,|;', line.rstrip('\r\n')) # make list based seperating elemntes based on tabs
        #         primers = line_list[1:] # primer sequences (contains forwards and reverse primer)
        #         primer_dict[line_list[0]]=primers # dict: {locus_name: [Forward, Reverse]}
        
        # self.params["Primers"] = {k : v for k, v in primer_dict.items() if k != ""}
        self.set_primers()
        primer_number = len(self.params["Primers"]["PrimerInformation"])
        #print(self.params["Primers"]["PrimerInformation"])
        
        # save dictionary with sample information
        sample_dir={}
        with open(SampleFile) as sample_file:
            for line in sample_file:
                line_list = re.split('\t|,|;', line.rstrip('\r\n')) # make list based seperating elemntes based on commas
                sampleNew = line_list[2] # new sample name
                sample_dir[line_list[0]]=sampleNew # save as dictionary: {old_sample_name: new_sample_name}
                #We make this dict in order to create a link between old and new sample_name
        self.params["SamplesOldNew"] = sample_dir
        
        self.params["DemultiplexParams"]["Number_Surviving_Sequences"] = 0  # 'i' for integer
        self.params["DemultiplexParams"]["Number_Abandonded_Sequences"] = 0

        
        if (self.performance_param):
            with multiprocessing.Pool(processes = multiprocessing.cpu_count()-1) as pool:
                processes = []
                args_list = [(fastq, input_directory, outDir) for fastq in os.listdir(input_directory)]
                results = pool.starmap_async(self.runDemultiplex_per_file, args_list)
                processes = results.get()
        else:
            for fastq in os.listdir(input_directory):
                self.runDemultiplex_per_file(fastq, input_directory, outDir)
        
        # if (self.performance_param):
        #     processes = []
        #     for fastq in os.listdir(input_directory):
        #         p = multiprocessing.Process(target=self.runDemultiplex_per_file, args=(fastq, input_directory, outDir))
        #         processes.append(p)
        #         p.start()
        #     for p in processes:
        #         p.join()
        # else:
        #     for fastq in os.listdir(input_directory):
        #         self.runDemultiplex_per_file(fastq, input_directory, outDir)
        
        number, number_filled, number_empty, total_size = self.check_results(outDir)
        
        print("Finished \n")
        #print(str(self.params["DemultiplexParams"]["Number_Surviving_Sequences"]))
            
        return number, number_filled, number_empty, total_size, primer_number, self.params["DemultiplexParams"]["Number_Surviving_Sequences"], self.params["DemultiplexParams"]["Number_Abandonded_Sequences"]
                
      
    def runDemultiplex_per_file(self, fastq, input_directory, outDir):
        newName = self.params["SamplesOldNew"][fastq.split('_')[0]] # get old sample name from file and retrieve new name
        print('processing sample ' + newName + ' (old name: ' + str(fastq.split('_')[0]) + ') ...')
        
        print(str(self.params["DemultiplexParams"]["Number_Surviving_Sequences"]))
        #self.sequences_used = []
        for locus, primers in self.params["Primers"]["PrimerInformation"].items():
            primerF = primers[0] # get sequence of primer forward
            primerR = self.rev_comp(primers[1]) # get sequence of primer reverse and reverse complement it
            print('processing primer '+ locus)
            out = os.path.normpath(outDir + '/' + newName + '_' + locus + '.fasta') #  get name of output file.
            fastq_path = os.path.normpath(input_directory + '/' + fastq) #path of the joined fastqfile from MergedOut folder
            self.extract_primer_combOv(fastq_path, primerF, primerR, out, locus)
    
    
    
    
    def parseFastqOv(self, fastq):
        with open(fastq) as file:
            result = {}
            #print()
            cn = None
            for i,line in enumerate(file):
                if i % 4 == 0:
                    cn = line.rstrip('\r\n')
                elif i % 4 in [1,2,3]:
                    result[cn] = result.get(cn,'')+line
        return result
    
    
    def rev_comp(self, seq):
        # dictionary with Key for complementing nucleotide
        code = {'A':'T','C':'G','T':'A','G':'C', 'Y':'R', 'R':'Y', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'B':'V', 'D':'H', 'H':'D', 'V':'B', 'N':'N'}
        result_sequence = '' # resulting string will be saved here
        for i in seq: # iterate through nucleotide in input string
            complement = code.get(i.upper(), i) # get complement. if it does not exist get the input nucleotide
            result_sequence = complement + result_sequence #add complemented nucleotide in the beggining of the output string. This way it is reversed
        return result_sequence
    
    def mismatch(self, seq_a, seq_b):
        # dictionary with Key of nucleotides that can match to each other. For example, an A in one string can matck with 'A','N','R','W','M' in the other string
        IUPAC = {'A':['A', 'N','R','W','M'], 'C':['C', 'N','Y','S','M'], 'G':['G', 'N','R','S','K'], 'T':['T', 'N','Y','W','K'], '-':['-'],
        'R':['A','G', 'N'], 'Y':['C','T', 'N'], 'S':['G','C', 'N'], 'W':['A','T', 'N'], 'K':['G','T', 'N'], 'M':['A','C', 'N'], 'B':['C','G','T', 'N'],
        'D':['A','G','T', 'N'], 'H':['A','C','T', 'N'], 'V':['A','C','G', 'N'], 'N':['A','T','G','C', 'N'], 'I':['A','T','G','C', 'N']}
        len1 = len(seq_a)
        len2 = len(seq_b)
        mismatches = 0 # mismatched is assumed to be 0. Only if a mismatch is found +1 is added
        for pos in range (0,min(len1,len2)): # in case sequences with different lengths are added, iterate only through indexxes of the small one.
            base_b = seq_b[pos] # get nucleotide of seq_b
            if base_b in IUPAC:
                base_b_IUPAC = IUPAC[base_b] # get possible matches of seq_b nucleotide in a list format
                # add mismatch only if corresponding nucleotode in seq_a is not found in the list of possible matches of seq_b
                if not seq_a[pos] in base_b_IUPAC:
                    mismatches+=1
            else:
                mismatches+=1
        return mismatches
    
    def extract_primer_combOv(self, fastq, primerF, primerR, out, locus):
        with open(out, 'w') as output_file:
            for rec, data in self.parseFastqOv(fastq).items(): #self.parseFastqOv(fastq) returns a dictionary with the keys the headers of the merged fastq files and the values being 3 strings (sequence, +, quality)
                read=data.split('\n')[0] # get sequence
                #qual=data.split('\n')[2] #  get quality
                motifF=str(read[:len(primerF)]) #  get region of the sequence coresponding to forward primer
                motifR=str(read[-len(primerR):]) #  get region of the sequence coresponding to reverse primer
                if int(self.mismatch(motifF, primerF)) <= int(self.params["MaxMismatch"]) and int(self.mismatch(motifR, primerR)) <= int(self.params["MaxMismatch"]):
                    if (locus in self.params["Primers"]["Boundaries"]):
                        for boundary in self.params["Primers"]["Boundaries"][locus]:
                            if (len(read) >= boundary[0] and len(read) <= boundary[1]):
                                output_file.write('>' + rec[1:] + '\n' + read + '\n')
                                self.params["DemultiplexParams"]["Number_Surviving_Sequences"] += 1
                                break
                    elif (len(read) >= int(self.params["MinLength"])):
                        #if (read in self.sequences_used):
                            #continue
                        #self.sequences_used.append(read)
                        self.params["DemultiplexParams"]["Number_Surviving_Sequences"] += 1
                        #print(str(self.params["DemultiplexParams"]["Number_Surviving_Sequences"]))
                        output_file.write('>' + rec[1:] + '\n' + read + '\n') # saves sequence in fasta format if mismatch with each primer is below maximum and length is above the defined length
                    else:
                        self.params["DemultiplexParams"]["Number_Abandonded_Sequences"] += 1
                else:
                    self.params["DemultiplexParams"]["Number_Abandonded_Sequences"] += 1
        #The demultiplexing process is not extracting a subset of the merged reads which fulfill the conditions of their beginning and end parts having mismatches with the primers below a certain self.max_mm and the reads themselves need to have lengths larger than at least self.min_length
        #This is done for each merged fastq file per primer triple (locus, forward_primer, reverse primer)
        #Example: 2 merged fastq files, 3 primer triples: 6 output files
    
    
    
    
    
    
    def getLengthStatistics(self):
        input_folder = os.path.normpath(self.params["OutputFolder"] + '/SeparatOut') #we want to get the full path of the folder, not just the name
        output_folder = os.path.normpath(self.params["OutputFolder"] + '/MarkerStatistics')
        output_folder2 = os.path.normpath(self.params["OutputFolder"] + '/AlleleLenghtCounts')
        self.lengths_counter = 0
        list_original_markers = []
        if (len(self.params["Primers"]["Boundaries"]) > 0):
            duplicate_folder = os.path.normpath(self.params["OutputFolder"] + '/MarkerStatisticsBoundaries')
            try:
                os.mkdir(duplicate_folder)
            except OSError:
                print ("Creation of the directory %s failed because it's already there." % duplicate_folder) # if it fail it produces this error message
            else:
                print ("Successfully created the directory %s " % duplicate_folder)
        for file in os.listdir(input_folder): #filenames in SeparatOut are of the form sample_locus.fasta
            print('get length profile from file ' + file + "\n")
            #samplename = file.rstrip('.fasta').split('_') #1st element of this list is samplename and then after a _ comes the locus name
            samplename = os.path.splitext(file)[0].split('_') #1st element of this list is samplename and then after a _ comes the locus name
            print(samplename)
            marker = '_'.join(samplename[1:]) # get marker name
            
            output_file_path = os.path.normpath(output_folder + '/' + marker + '_' + samplename[0] + '_.Statistics') # output name: marker_sampleName_.Statistics
            #output_file_path2 = os.path.normpath(output_folder2 + '/' + marker + '_' + samplename[0] + '_.txt') # output name: marker_sampleName_.Statistics
            input_file_path = os.path.normpath(input_folder + '/' + file)
            self.getLengthProfilePerSample(input_file_path, output_file_path, marker, samplename[0]) # save statistic file
            
            
            if (marker in self.params["Primers"]["Boundaries"] and ("dupl" in marker)):
                pattern = r'dupl[1-9]|dupl10'
                marker_original = re.sub(pattern, '', marker)
                file_original = re.sub(pattern, '', file)
                print("Duplicates")
                print(marker_original)
                print(marker)
                print(file_original)
                #if (marker_original not in list_original_markers):
                print(file)
                print("Hit for " + marker_original)
                output_file_path_dupl = os.path.normpath(duplicate_folder + '/' + marker_original + '_' + samplename[0] + '_.Statistics')
                input_file_path = os.path.normpath(input_folder + '/' + file_original)
                self.getLengthProfilePerSample(input_file_path, output_file_path_dupl, marker_original, samplename[0])
                list_original_markers.append(marker_original)
                
                output_file_path_dupl = os.path.normpath(duplicate_folder + '/' + marker + '_' + samplename[0] + '_.Statistics')
                input_file_path = os.path.normpath(input_folder + '/' + file)
                self.getLengthProfilePerSample(input_file_path, output_file_path_dupl, marker, samplename[0])
            
            # i = 1
            # if ("dupl" in marker):
            #     output_file_path_dupl = os.path.normpath(duplicate_folder + '/' + marker + '_' + samplename[0] + '_.Statistics')
            #     input_file_path = os.path.normpath(input_folder + '/' + file)
            #     self.getLengthProfilePerSample(input_file_path, output_file_path_dupl, marker, samplename[0])
                
            #     marker = marker.replace("dupl" + str(i), "")
            #     output_file_path_dupl = os.path.normpath(duplicate_folder + '/' + marker + '_' + samplename[0] + '_.Statistics')
            #     input_file_path = os.path.normpath(input_folder + '/' + file)
            #     self.getLengthProfilePerSample(input_file_path, output_file_path_dupl, marker, samplename[0])
                
            
            # output_file_path = os.path.normpath(output_folder + '/' + marker + '_' + samplename[0] + '_.Statistics') # output name: marker_sampleName_.Statistics
            # #output_file_path2 = os.path.normpath(output_folder2 + '/' + marker + '_' + samplename[0] + '_.txt') # output name: marker_sampleName_.Statistics
            # input_file_path = os.path.normpath(input_folder + '/' + file)
            # self.getLengthProfilePerSample(input_file_path, output_file_path, marker, samplename[0]) # save statistic file
            # #self.getLengthProfilePerSample(input_file_path, output_file_path2, marker, samplename[0])
        
        print(self.params["LengthStatistics"]["LowestLength"])
        print()
        print("Belong to filename " + self.params["LengthStatistics"]["FileNameLowest"] + "\n")
        print()
        print(self.params["LengthStatistics"]["LargestLength"])
        print()
        print("Belong to filename " + self.params["LengthStatistics"]["FileNameLargest"] + "\n")
        print()
        
        number, number_filled, number_empty, total_size = self.check_results(output_folder)
        
        return self.params["LengthStatistics"]["LowestLength"], self.params["LengthStatistics"]["FileNameLowest"], self.params["LengthStatistics"]["LargestLength"], self.params["LengthStatistics"]["FileNameLargest"], number, number_filled, number_empty, total_size
        
    def getLengthProfilePerSample(self, fasta_path, out_path, marker, sample):
        # get dictionary {length: nr of reads}
        result = {}
        count_total = 0
        with open(fasta_path) as fasta_file:
            fasta_file_base = os.path.basename(os.path.normpath(fasta_path))
            with open(out_path, 'w') as output_file:
                for i, line in enumerate(fasta_file, start=1):
                    if i % 2 == 0:
                        seq_length = len(line.rstrip('\r\n'))
                        # if (marker in self.params["Primers"]["Boundaries"] and ((seq_length < int(self.params["Primers"]["Boundaries"][marker][0])) or seq_length > int(self.params["Primers"]["Boundaries"][marker][1]))):
                        #     print("Ignoring length " + str(seq_length) + " for marker " + marker + " and sample " + str(sample) + "\n")
                        #     continue
                        if (seq_length not in result.keys()):
                            result[seq_length] = 0
                        result[seq_length] += 1
                
                #list of counts and sort them; get dictionary of lenghs per count
                resultCount = {} # lengths per count
                counts = [] # count list
                for length, count in result.items():
                    counts.append(count)
                    if (count not in resultCount.keys()):
                        resultCount[count] = []
                    resultCount[count].append(length)
                    #resultCount[c] = resultCount.get(c, []) + [l]
                counts = list(set(counts))
                counts.sort()
                
                if result:
                    length_with_max_counts = max(result, key=result.get)
                    length_with_min_counts = min(result, key=result.get)
        
                # iterate through count list
                for count in counts:
                    lengths = resultCount[count] # get all lengths
                    lengths.sort() # sort lengths
                    for length in lengths:
                        if (self.lengths_counter == 0):
                            self.params["LengthStatistics"]["LowestLength"] = length
                            self.params["LengthStatistics"]["LargestLength"] = length
                            #self.lowest_length = length
                            #self.largest_length = length
                        elif (self.lengths_counter > 0):
                            if (length < self.params["LengthStatistics"]["LowestLength"]):
                                self.params["LengthStatistics"]["LowestLength"] = length
                                self.params["LengthStatistics"]["FileNameLowest"] = out_path
                            elif (length > self.params["LengthStatistics"]["LargestLength"]):
                                self.params["LengthStatistics"]["LargestLength"] = length
                                self.params["LengthStatistics"]["FileNameLargest"] = out_path
                        self.lengths_counter += 1
                        if (marker in self.params["Primers"]["Boundaries"]):
                            line = str(length) + ' ' + str(count)
                            for boundary in self.params["Primers"]["Boundaries"][marker]:
                                line += ' ' + str(boundary[0]) + ' ' + str(boundary[1])
                            output_file.write(line + ' \n')
                            print("boundaries: \n")
                            print(line)
                        else:
                            output_file.write(str(length) + ' ' + str(count) + '\n') #  save: length count
                            print("No boundaries: \n")
                            print(str(length) + '' + str(count))
                        
            # fasta_file_name = fasta_file_base
            # self.lengths_per_separatfile[fasta_file_name] = {}
            # for length, count in result.items():
            #     count_total += count
            #     self.lengths_per_separatfile[fasta_file_name][length] = count
            
            # if result:
            #     count_frequency_of_max = float(result[length_with_max_counts]/count_total) * 100
            #     count_frequency_of_min = float(result[length_with_min_counts]/count_total) * 100
            #     if (count_frequency_of_max > self.highest_countfrequency):
            #         self.length_highest_countfrequency = length_with_max_counts
            #         self.highest_countfrequency = count_frequency_of_max
            #         self.filename_highestfrequence = out_path
            #     if (count_frequency_of_min > self.lowest_countfrequency):
            #         self.length_lowest_countfrequency = length_with_min_counts
            #         self.lowest_countfrequency = count_frequency_of_min
            #         self.filename_lowestfrequence = out_path
        #print("Both dict and reversedict: ")
        print()
        #print(result) #form: length_of_sequence : abundance/count
        print()
        #print(resultCount) #form: abundance/count : length_of_sequence
        print()
    
    
    
    
    
    
    def GenotypeLength(self):
        if (self.params["Operatingsystem"].lower() == "windows"):
            R = os.path.normpath(self.params["RExecutable"])
            
        else:
            R = "Rscript"
        
        list_error = []
        #GenotypeLength(RExDip, OutDir + 'MarkerStatistics/', minNrReads, '0.7', OutDir + 'MarkerPlots/', LengthWindow)
        input_folder = os.path.normpath(self.params["OutputFolder"] + '/MarkerStatistics')
        output_folder = os.path.normpath(self.params["OutputFolder"] + '/MarkerPlots')
        Rscript = None
        if (self.params["Ploidy"] == "diploid"):
            #Rscript = self.scriptRDiploid_path
            Rscript = self.params["Executables"]["scriptRDiploid"]
        elif (self.params["Ploidy"] == "haploid"):
            #Rscript = self.scriptRHaploid_path
            Rscript = self.params["Executables"]["scriptRHaploid"]
            
        print()
        print(Rscript)
        # make code
        markerplots_file_path = os.path.normpath(output_folder + '/' + 'markerplots.pdf ')
        markermatrix_file_path = os.path.normpath(output_folder + '/' + 'markermatrix.csv ')
        markerlist_file_path = os.path.normpath(output_folder + '/' + 'markerlist.csv ')
        #code = R + ' --vanilla ' + Rscript + ' ' + input_folder +  ' ' + str(self.minCount) + ' ' + markerplots_file_path + str(self.alpha) + ' ' + markermatrix_file_path + markerlist_file_path
        #code = R + ' --vanilla ' + Rscript + ' ' + input_folder + '/' +  ' ' + str(self.minCount) + ' ' + output_folder + '/' + 'markerplots.pdf ' + str(self.alpha) + ' ' + output_folder + '/' + 'markermatrix.csv ' + output_folder + '/' + 'markerlist.csv ' + ' '.join(str(self.lengthWindow))
        #code = R + ' --vanilla ' + Rscript + ' ' + input_folder +  ' ' + str(self.minCount) + ' ' + output_folder + 'markerplots.pdf ' + str(self.alpha) + ' ' + output_folder + 'markermatrix.csv ' + output_folder + 'markerlist.csv' + ' '.join(str(self.lengthWindow))
        #code = R + ' --vanilla ' + Rscript + ' ' + input_folder +  ' ' + str(self.minCount) + ' ' + markerplots_file_path + str(self.alpha) + ' ' + markermatrix_file_path + markerlist_file_path + ' '.join(self.lengthWindow)
        print(self.params["LengthWindow"])
        print(len(self.params["LengthWindow"]))
        liste = self.params["LengthWindow"].split(",")
        print(len(liste))
        
        
        #liste = []
        #for element in self.params["LengthWindow"]:
            #liste.append(str(element))
        #self.params["LengthWindow"] = liste
        if (self.params["Operatingsystem"].lower() == "windows"):
            # duplicate_folder = os.path.normpath(self.params["OutputFolder"] + '/MarkerStatisticsBoundaries')
            # try:
            #     os.mkdir(duplicate_folder)
            # except OSError:
            #     print ("Creation of the directory %s failed because it's already there." % duplicate_folder) # if it fail it produces this error message
            # else:
            #     print ("Successfully created the directory %s " % duplicate_folder)
            code = r'''{Rpath} --vanilla {Rscript} {in_folder} {minCount} {plots_path} {alpha} {matrix_path} {marker_path} {lengthWindow}'''.format(   #{marker_path} between matrix_path and LengthWindow
                Rpath=R,
                Rscript=Rscript,
                in_folder=input_folder,
                minCount=self.params["MinCount"],
                plots_path=markerplots_file_path,
                alpha=self.params["ConsensusThreshold"],
                matrix_path=markermatrix_file_path,
                marker_path=markerlist_file_path,
                lengthWindow= ' '.join(liste) #' '.join(self.params["LengthWindow"])
            )
            os_code = subprocess.run(code)
            if os_code != 0:
                print(f"The command exited with error code: {os_code}")
                #list_error.append(os_code)
            #else:
                #return os_code
        
            for error in list_error:
                if (error != 0):
                    return error
                
        else:
            code = R + ' --vanilla ' + Rscript + ' ' + input_folder +  ' ' + str(self.params["MinCount"]) + ' ' + markerplots_file_path + str(self.params["ConsensusThreshold"]) + ' ' + markermatrix_file_path + markerlist_file_path + ' '.join(liste)
            subprocess.call(code, shell=True)
            
        
        if (len(self.params["Primers"]["Boundaries"]) > 0):
            input_folder = os.path.normpath(self.params["OutputFolder"] + '/MarkerStatisticsBoundaries')
            markerplots_file_path = os.path.normpath(output_folder + '/' + 'markerplotsDuplicates.pdf ')
            markermatrix_file_path = os.path.normpath(output_folder + '/' + 'markermatrixBoundaries.csv ')
            markerlist_file_path = os.path.normpath(output_folder + '/' + 'markerlistBoundaries.csv ')
            if (self.params["Operatingsystem"].lower() == "windows"):
                code = r'''{Rpath} --vanilla {Rscript} {in_folder} {minCount} {plots_path} {alpha} {matrix_path} {marker_path} {lengthWindow}'''.format(   #{marker_path} between matrix_path and LengthWindow
                    Rpath=R,
                    Rscript=Rscript,
                    in_folder=input_folder,
                    minCount=self.params["MinCount"],
                    plots_path=markerplots_file_path,
                    alpha=self.params["ConsensusThreshold"],
                    matrix_path=markermatrix_file_path,
                    marker_path=markerlist_file_path,
                    lengthWindow= ' '.join(liste) #' '.join(self.params["LengthWindow"])
                )
                os_code = subprocess.run(code)
                if os_code != 0:
                    print(f"The command exited with error code: {os_code}")
                    #list_error.append(os_code)
                #else:
                    #return os_code
            
                for error in list_error:
                    if (error != 0):
                        return error
                    
            else:
                code = R + ' --vanilla ' + Rscript + ' ' + input_folder +  ' ' + str(self.params["MinCount"]) + ' ' + markerplots_file_path + str(self.params["ConsensusThreshold"]) + ' ' + markermatrix_file_path + markerlist_file_path + ' '.join(liste)
                subprocess.call(code, shell=True)
            
            print(code)
            
        #code = R + ' --vanilla ' + Rscript + ' ' + input_folder +  ' ' + str(self.params["MinCount"]) + ' ' + markerplots_file_path + str(self.params["ConsensusThreshold"]) + ' ' + markermatrix_file_path + markerlist_file_path + ' '.join(liste)
        #print(code) # print code
        #subprocess.call(code, shell=True) # call code
        #os_code = os.system(code)
        #subprocess.call(code, shell=True) # run program
        
        
        
        # os_code = subprocess.run(code)
        # if os_code != 0:
        #     print(f"The command exited with error code: {os_code}")
        #     #list_error.append(os_code)
        # #else:
        #     #return os_code
    
        # for error in list_error:
        #     if (error != 0):
        #         return error
        
        self.open_file(os.path.normpath(self.params["OutputFolder"] + "/Markerplots/markerplots.pdf"))
        
        return 0
    
    
    
    
    def check_results(self, folder_path):
        print()
        number = 0
        number_empty = 0
        number_filled = 0
        total_size = 0
        if os.path.exists(folder_path):
            for root, dirs, files in os.walk(folder_path):
                number += len(files)
                for file in files:
                    file_path = os.path.join(root, file)
                    file_size = os.path.getsize(file_path)  # Get the size of the file
                    total_size += file_size
                    if (not(os.path.getsize(file_path) == 0)):
                        number_filled += 1
                    else:
                        number_empty += 1
        return number, number_filled, number_empty, total_size

    

    def open_file(self, path):
        if (self.params["Operatingsystem"] == "windows"):
            os.startfile(path)
        elif (self.params["Operatingsystem"] == "linux"):
            subprocess.call(["xdg-open", path])
    


        
        



if __name__ == "__main__":
    main()