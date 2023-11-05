# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 16:59:42 2023

@author: Sebastian
"""

import os, subprocess, re, sys
from datetime import datetime
from collections import Counter
import csv

def main():
    
    paramter_file_path = os.path.normpath(os.getcwd() + '/Parameter.txt')
    if not(os.path.exists(paramter_file_path)):
        raise Exception("There is no paramter file with the name Parameter.txt")
    
    with open(paramter_file_path) as parameters:
        for line in parameters:
            if (not line) or line.startswith("#"):
                continue
            if (line.startswith("DataDir")):
                raw_data_dir = line.split("=")[1].rstrip("\r\n").strip()
                if (not raw_data_dir):
                    raise Exception("No name for data_dir in parameter file!")
                continue
            if (line.startswith("PrimerFile")):
                primerfile = line.split("=")[1].rstrip("\r\n").strip()
                if (not primerfile):
                    raise Exception("No name for primerfile in parameter file!")
                continue
            if (line.startswith("SampleFile")):
                sample_sheet = line.split("=")[1].rstrip("\r\n").strip()
                if (not sample_sheet):
                    raise Exception("No name for sample file in parameter file!")
                continue
            if (line.startswith("WorkingDir")):
                output_folder = line.split("=")[1].rstrip("\r\n").strip()
                if (not output_folder):
                    output_folder = "Output_default"
                continue
            if (line.startswith("Ploidy")):
                datatype = line.split("=")[1].rstrip("\r\n").strip()
                if (not datatype):
                    datatype = "diploid"
                continue
            
            
            if (line.startswith("numeberMM")):
                max_mismatches = line.split("=")[1].rstrip("\r\n").strip()
                if (not max_mismatches):
                    max_mismatches = 2
                else:
                    max_mismatches = int(max_mismatches)
                continue
            if (line.startswith("minCount")):
                minCount = line.split("=")[1].rstrip("\r\n").strip()
                if (not minCount):
                    minCount = 10
                else:
                    minCount = int(minCount)
                continue
            if (line.startswith("minLength")):
                min_length = line.split("=")[1].rstrip("\r\n").strip()
                if (not min_length):
                    min_length = 290
                else:
                    min_length = int(min_length)
                continue
            if (line.startswith("consThreashold")):
                alpha = line.split("=")[1].rstrip("\r\n").strip()
                if (not alpha):
                    alpha = 0.7
                else:
                    alpha = float(alpha)
                continue
            if (line.startswith("LengthWindow")):
                values = line.split("=")[1].rstrip("\r\n").strip().split(",")
                if (not values):
                    lengthWindow = ['310', '600']
                else:
                    lengthWindow = []
                    lengthWindow.append(values[0])
                    lengthWindow.append(values[1])
                continue
            if (line.startswith("AlleleList")):
                allele_list = line.split("=")[1].rstrip("\r\n").strip()
                continue
            
            if (line.startswith("operating_system")):
                par_os = line.split("=")[1].rstrip("\r\n").strip()
                if (not par_os):
                    par_os = "linux"
                continue
            if (line.startswith("unique_identifier")):
                unique_name = line.split("=")[1].rstrip("\r\n").strip()
                if (not unique_name):
                    unique_name = "Default"
                continue
            
            
    
    print(raw_data_dir)
    print(primerfile)
    print(sample_sheet)
    print(output_folder)
    print(datatype)
    
    print(max_mismatches)
    print(min_length)
    print(minCount)
    print(alpha)
    print(lengthWindow)
    
    print(allele_list)
    print(par_os)
    print(unique_name)
    
    print()
    print(datatype)
    print()

    print("AlleList: " + allele_list)
    
    project_folder = "Default1"
    
    
    Amplicon_inst = amplicon_class(project_folder, output_folder, primerfile, raw_data_dir, sample_sheet, max_mismatches, min_length, minCount, alpha, lengthWindow, datatype, par_os)
    
    print("Running!")
    Amplicon_inst.print_folderfiles_inputs()
    print()
    Amplicon_inst.set_executables()
    print()
    Amplicon_inst.set_samplelist()
    print()
    Amplicon_inst.checkInputDir()
    print()
    Amplicon_inst.runTrimomatic()
    print()
    Amplicon_inst.runUsearchMergePairs()
    print()
    Amplicon_inst.runDemultiplex()
    print()
    print("Now length statistics:")
    Amplicon_inst.getLengthStatistics()
    print()
    print()
    #For Windows:
    r_exec_path = "C:\\Programme\\R\\R-4.3.1\\bin\\Rscript.exe"
    print("Now getting genotypes based on RScript:")
    #Amplicon_inst.GenotypeLength(r_exec_path)
    
    
    
    
    
    
    return 0


class amplicon_class():
    def __init__(self, project_folder, output_folder, primerfile, raw_dir, sample_sheet, max_mm, min_length, minCount, alpha, lengthWindow, datatype, par_os):
        self.workspace_path = os.getcwd() #Setting the workspace to the directory where the script is located
        
        self.par_os = par_os
        
        
        self.project_path = os.path.normpath(self.workspace_path + "/" + project_folder)
        if (project_folder == ""):
            self.project_folder = "Default"
        else:
            self.project_folder = project_folder
            
        try:
            os.mkdir(self.project_path)
        except OSError:
            print ("Creation of the directory %s failed because it's already there." % self.project_folder) # if it fail it produces this error message
        else:
            print ("Successfully created the directory %s " % self.project_folder)
        
        
        
        self.output_folder_path = os.path.normpath(self.project_path +  "/" + output_folder)
        if (project_folder == ""):
            self.output_folder = "DefaultOutput"
        else:
            self.output_folder = output_folder
            
        try:
            os.mkdir(self.output_folder_path)
        except OSError:
            print ("Creation of the directory %s failed because it's already there." % self.output_folder) # if it fail it produces this error message
        else:
            print ("Successfully created the directory %s " % self.output_folder)
        
        for folder in ['QC', 'SeparatOut', 'MergedOut', 'MarkerStatistics', 'MarkerPlots']:
            try:
                os.mkdir(self.output_folder_path + '/' + folder)
            except OSError:
                print ("Creation of the directory %s failed because it's already there." % folder) # if it fail it produces this error message
            else:
                print ("Successfully created the directory %s " % folder)
        
        
        
        
        bin_path = os.path.normpath(self.workspace_path + '/bin')
        if os.path.exists(bin_path):
            self.bin_path = bin_path
            self.bin = "bin"
            #print("Bin path: " + self.bin + "\n")
        else:
            raise Exception("There is no bin! Pipeline cannot start!")
            
            
        raw_dir_path = os.path.normpath(self.workspace_path + '/' + raw_dir)
        if os.path.exists(raw_dir_path):
            self.raw_data_dir = raw_dir
            self.raw_data_dir_path = raw_dir_path
            #print("Raw data path: " + self.raw_data_dir + "\n")
        else:
            raise Exception("There is no directory containing raw data! Pipeline cannot start!")
            
            
        primerfile_path = os.path.normpath(self.workspace_path + '/' + primerfile)
        if os.path.exists(primerfile_path):
            self.primerfile_path = primerfile_path
            self.primerfile = primerfile
            #print("Primerfile: " + self.primerfile + "\n")
        else:
            raise Exception("There is no Primerfile!")
            
            
        sample_sheet_path = os.path.normpath(self.workspace_path + '/' + sample_sheet)
        if os.path.exists(sample_sheet_path):
            self.sample_sheet_path = sample_sheet_path
            self.sample_sheet = sample_sheet
            #print("Samplesheet (csv): " + self.sample_sheet + "\n")
        else:
            raise Exception("There is no Samplesheet!")
            
            
        
        self.trim_dir = ""
        self.trim_dir_path = ""
        self.adapters_dir = ""
        self.adapters_dir_path = ""
        
        self.trimmo_exec = ""
        self.trimmo_exec_path = ""
        self.usearch_exec = ""
        self.usearch_exec_path = ""
        self.scriptRDiploid_path = ""
        self.scriptRDiploid = ""
        self.scriptRHaploid_path = ""
        self.scriptRHaploid = ""
        
        self.sample_list = []
        
        self.max_mm  = max_mm
        self.min_length = min_length
        self.minCount = minCount
        self.alpha = alpha
        self.lengthWindow = lengthWindow
        self.datatype = datatype
        
        
        self.analysis_folder_path = os.path.normpath(self.output_folder_path +  "/Analysis")
        self.analysis_folder = "Analysis"
        try:
            os.mkdir(self.analysis_folder_path)
        except OSError:
            print ("Creation of the directory %s failed because it's already there." % self.analysis_folder) # if it fail it produces this error message
        else:
            print ("Successfully created the directory %s " % self.analysis_folder)
        
        self.analysis_file_path = os.path.normpath(self.analysis_folder_path + "/analysis.txt")
        
        
        self.lengths_per_separatfile =  {}
        
    def return_params(self):
        return self.project_folder, self.output_folder, self.sample_sheet, self.primerfile, self.raw_data_dir
    
    def return_r_path_diploid(self):
        return self.scriptRDiploid_path
    
    def print_folderfiles_inputs_paths(self): #printing out all the paths of nessecary folders and files:
        print("Workspace path: " + self.workspace_path + "\n")
        print("Bin path: " + self.bin_path + "\n")
        print("Raw data path: " + self.raw_data_dir_path + "\n")
        print("Primerfile path: " + self.primerfile_path + "\n")
        print("Samplesheet (csv) path: " + self.sample_sheet_path + "\n")
        
    def print_folderfiles_inputs(self): #printing out all the paths of nessecary folders and files:
        print("Workspace: " + self.workspace_path + "\n")
        print("Bin: " + self.bin + "\n")
        print("Raw data: " + self.raw_data_dir + "\n")
        print("Primerfile: " + self.primerfile + "\n")
        print("Samplesheet (csv): " + self.sample_sheet + "\n")
        
        
        
        
    
    def set_executables(self):
        for subdir, dirs, files in os.walk(self.bin):
            if (os.path.basename(subdir).split('-')[0] == 'Trimmomatic'):
                self.trim_dir = os.path.basename(subdir)
                self.trim_dir_path = os.path.normpath(self.bin + '/' + self.trim_dir)
            elif (os.path.basename(subdir) == 'adapters'):
                self.adapters_dir = os.path.basename(subdir)
                self.adapters_dir_path = os.path.normpath(self.bin + '/' + self.adapters_dir)
                  
            for file in files:
                if (file.split("-")[0] == "trimmomatic" and ".jar" in file):
                    self.trimmo_exec = file
                    self.trimmo_exec_path = os.path.normpath(self.trim_dir_path + '/' + file)
                elif file.startswith('usearch'):
                    if (self.par_os == "windows" and "win" in file):
                        self.usearch_exec = file
                        self.usearch_exec_path = os.path.normpath(self.bin + '/' + file)
                    elif (self.par_os == "linux" and not("win" in file)):
                        self.usearch_exec = file
                        self.usearch_exec_path = os.path.normpath(self.bin + '/' + file)
                elif file.startswith('Rscript_Markerlength_STUTTER_Color_BaryzentricMinsize20_notVerbose'):
                    self.scriptRDiploid = file
                    self.scriptRDiploid_path = os.path.normpath(self.bin + '/'  + file)
                elif file.startswith('Rscript_Markerlength_Haploid'):
                    self.scriptRHaploid = file
                    self.scriptRHaploid_path = os.path.normpath(self.bin + '/'  + file)
      
    def return_executables(self):
        return self.trimmo_exec, self.usearch_exec, self.scriptRDiploid, self.scriptRHaploid
    
    def run_script(self):
        print("Running!")
        self.print_folderfiles_inputs()
        print()
        self.set_executables()
        print()
        self.set_samplelist()
        print()
        self.checkInputDir()
        print()
        self.runTrimomatic()
        print()
        #self.runUsearchMergePairs()
    
    def print_executables(self):
        print("Trimmomatic exec: " + self.trimmo_exec + "\n")
        print("Usearch exec: " + self.usearch_exec + "\n")
        print("RDiploid exec: " + self.scriptRDiploid + "\n")
        print("RHaploid exec: " + self.scriptRHaploid + "\n")
        
    def set_samplelist(self):
        rows_new = []
        
        with open(self.sample_sheet_path) as file:
            for line in file:
                row = []
                self.sample_list.append(re.split(',|\t|;', line)[0])
                row.append(re.split(',|\t|;', line)[0])
                row.append(re.split(',|\t|;', line)[1])
                row.append(re.split(',|\t|;', line)[1].replace("_", "").replace("-", ""))
            
                rows_new.append(row)
        
        print(rows_new)
        with open(self.sample_sheet_path, 'w', newline='') as new_samplesheet:
            writer = csv.writer(new_samplesheet, delimiter =',', quoting=csv.QUOTE_MINIMAL)
            for row in rows_new:
                row = [element.replace("\n", "") for element in row]
                print(row)
                writer.writerow(row)
        
        
    def checkInputDir(self):
        FilesData = {}
        for subdir, dirs, files in os.walk(self.raw_data_dir): #instead of: for fastq in os.listdir(self.raw_data_dir):
            for fastq in files: 
                sample = fastq.split('_')[0]
                if sample in self.sample_list:
                    if (sample not in FilesData.keys()):
                        FilesData[sample] = []
                    FilesData[sample].append(os.path.normpath(subdir + "/" + fastq))
                            
        return FilesData
    
    
    def runTrimomatic(self):
        #Intput files can be either / or \\ if given the path but the output files always need to be with / slashes, even on windows
        FileData = self.checkInputDir()
        print(FileData)
        list_error = []
        for sample, files in FileData.items():
            files.sort() # get R1 and R2 list and sort them to make sure that R1 always shows first
            R1 = files[0] # get path of R1 input file
            R2 = files[1]
            
            if (self.par_os == "windows"):
                outFileR1 = R1.split('\\')[-1].split('.')[0]
                outFileR2 = R2.split('\\')[-1].split('.')[0]
            elif (self.par_os == "linux"):
                outFileR1 = R1.split('/')[-1].split('.')[0]
                outFileR2 = R2.split('/')[-1].split('.')[0]
            print(outFileR1)
            
            R1paired = self.project_folder + "/" + self.output_folder + "/QC/" + outFileR1 + "_QTpaired.fastq"
            R1unpaired = self.project_folder + "/" + self.output_folder + "/QC/" + outFileR1 + "_QTunpaired.fastq"
            R2paired = self.project_folder + "/" + self.output_folder + "/QC/" + outFileR2 + "_QTpaired.fastq"
            R2unpaired = self.project_folder + "/" + self.output_folder + "/QC/" + outFileR2 + "_QTunpaired.fastq"
            
            adapter = "bin/" + self.trim_dir + "/" + self.adapters_dir + "/TrueSeqAdaptersInUsage.fa"
            
            exe =  "bin/" + self.trim_dir + "/" + self.trimmo_exec
            
            code = 'java -jar ' + exe + ' PE ' + R1 + ' ' + R2 + ' ' + R1paired + ' ' + R1unpaired + ' ' + R2paired + ' ' + R2unpaired + ' ILLUMINACLIP:' + adapter + ':2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20'
            print('running: ' + code) # print code
            #subprocess.call(code, shell=True)
            os_code = os.system(code)
            list_error.append(os_code)
            if os_code != 0:
                print(f"The command exited with error code: {os_code}")
                #list_error.append(os_code)
            #else:
                #return os_code
        
        for error in list_error:
            if (error != 0):
                return error
        
        return 0
                
            
        
        
            
            
    def runUsearchMergePairs(self):
        # iterate through input fastq files
        list_error = []
        input_directory = os.path.normpath(self.output_folder_path + '/QC') #The input for usearch are the trimmed sequences from the QC folder
        output_folder = os.path.normpath(self.output_folder_path + '/MergedOut')
        #For example: If we have 2 forward/reverse read files (meaning 4 fastq files in total), we get 2  joined fastq files
        for file in os.listdir(input_directory):
            # only consider R1 files and only paired files, unpaired are ignored
            if file.endswith('R1_001_QTpaired.fastq'):
                R1 = os.path.normpath(input_directory + '/' + file) # get path of R1 input file
                R2 = R1.replace("_R1_", "_R2_") # get path of R2 input file
                out = os.path.normpath(output_folder + '/' + file.split('_')[0] + '_joined.fastq') # get path of file to save merge reads
                exe = os.path.normpath("bin/" + self.usearch_exec) #self.usearch_exec_path
                code = exe + ' -fastq_mergepairs ' + R1 + ' -reverse ' + R2 + ' -fastqout ' + out + ' -fastq_pctid 80 -fastq_maxdiffs 40'
                print('running: ' + code) # print code
                print()
                os_code = os.system(code)
                #subprocess.call(code, shell=True) # run program
                list_error.append(os_code)
                if os_code != 0:
                    print(f"The command exited with error code: {os_code}")
                    #list_error.append(os_code)
                #else:
                    #return os_code
            
        for error in list_error:
            if (error != 0):
                return error
            
        return 0
    
    
    
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
            base_b_IUPAC = IUPAC[base_b] # get possible matches of seq_b nucleotide in a list format
            # add mismatch only if corresponding nucleotode in seq_a is not found in the list of possible matches of seq_b
            if not seq_a[pos] in base_b_IUPAC:
                mismatches+=1
        return mismatches
    
    def extract_primer_combOv(self, fastq, primerF, primerR, out):
        with open(out, 'w') as output_file:
            for rec, data in self.parseFastqOv(fastq).items(): #self.parseFastqOv(fastq) returns a dictionary with the keys the headers of the merged fastq files and the values being 3 strings (sequence, +, quality)
                read=data.split('\n')[0] # get sequence
                qual=data.split('\n')[2] #  get quality
                motifF=str(read[:len(primerF)]) #  get region of the sequence coresponding to forward primer
                motifR=str(read[-len(primerR):]) #  get region of the sequence coresponding to reverse primer
                if int(self.mismatch(motifF, primerF)) <= int(self.max_mm) and int(self.mismatch(motifR, primerR)) <= int(self.max_mm):
                    if (len(read) >= int(self.min_length)):
                        output_file.write('>' + rec[1:] + '\n' + read + '\n') # saves sequence in fasta format if mismatch with each primer is below maximum and length is above the defined length
        #The demultiplexing process is not extracting a subset of the merged reads which fulfill the conditions of their beginning and end parts having mismatches with the primers below a certain self.max_mm and the reads themselves need to have lengths larger than at least self.min_length
        #This is done for each merged fastq file per primer triple (locus, forward_primer, reverse primer)
        #Example: 2 merged fastq files, 3 primer triples: 6 output files
    
    def runDemultiplex(self): 
        input_directory = os.path.normpath(self.output_folder_path + '/MergedOut')
        demultiplexFile = self.primerfile #demultiplexing based on primers
        SampleFile = self.sample_sheet
        
        outDir = os.path.normpath(self.output_folder_path + '/SeparatOut')
        
        count = 0

        # save dictionary with primer information
        primer_dict={}
        with open(demultiplexFile) as primer_file:
            for line in primer_file:
                line_list = re.split('\t|,|;', line.rstrip('\r\n')) # make list based seperating elemntes based on tabs
                primers = line_list[1:] # primer sequences (contains forwards and reverse primer)
                primer_dict[line_list[0]]=primers # dict: {locus_name: [Forward, Reverse]}
            

        # save dictionary with sample information
        sample_dir={}
        with open(SampleFile) as sample_file:
            for line in sample_file:
                line_list = re.split('\t|,|;', line.rstrip('\r\n')) # make list based seperating elemntes based on commas
                sampleNew = line_list[2] # new sample name
                sample_dir[line_list[0]]=sampleNew # save as dictionary: {old_sample_name: new_sample_name}
                #We make this dict in order to create a link between old and new sample_name


        for fastq in os.listdir(input_directory):
            print()
            newName = sample_dir[fastq.split('_')[0]] # get old sample name from file and retrieve new name
            print('processing sample ' + newName + ' (old name: ' + str(fastq.split('_')[0]) + ') ...')
            
            for locus, primers in primer_dict.items():
                primerF = primers[0] # get sequence of primer forward
                primerR = self.rev_comp(primers[1]) # get sequence of primer reverse and reverse complement it
                print('processing primer '+ locus)
                out = os.path.normpath(outDir + '/' + newName + '_' + locus + '.fasta') #  get name of output file.
                fastq_path = os.path.normpath(input_directory + '/' + fastq) #path of the joined fastqfile from MergedOut folder
                self.extract_primer_combOv(fastq_path, primerF, primerR, out)
            count += 1
        return count
                
    
    
    
    def getLengthStatistics(self):
        input_folder = os.path.normpath(self.output_folder_path + '/SeparatOut') #we want to get the full path of the folder, not just the name
        output_folder = os.path.normpath(self.output_folder_path + '/MarkerStatistics')
        for file in os.listdir(input_folder):
            print('get length profile from file ' + file + "\n")
            samplename = file.rstrip('.fasta').split('_') #before it was called nameED, why?
            #print(samplename)
            marker = '_'.join(samplename[1:]) # get marker name
            #print(marker)
            output_file_path = os.path.normpath(output_folder + '/' + marker + '_' + samplename[0] + '_.Statistics') # output name: marker_sampleName_.Statistics
            input_file_path = os.path.normpath(input_folder + '/' + file)
            self.getLengthProfilePerSample(input_file_path, output_file_path) # save statistic file
    
    def getLengthProfilePerSample(self, fasta_path, out_path):
        # get dictionary {length: nr of reads}
        result = {}
        with open(fasta_path) as fasta_file:
            fasta_file_base = os.path.basename(os.path.normpath(fasta_path))
            with open(out_path, 'w') as output_file:
                for i, line in enumerate(fasta_file, start=1):
                    if i % 2 == 0:
                        seq_length = len(line.rstrip('\r\n'))
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
                
        
                # iterate through count list
                for count in counts:
                    lengths = resultCount[count] # get all lengths
                    lengths.sort() # sort lengths
                    for length in lengths:
                        output_file.write(str(length) + ' ' + str(count) + '\n') #  save: length count
        
            fasta_file_name = fasta_file_base
            self.lengths_per_separatfile[fasta_file_name] = {}
            for length, count in result.items():
                self.lengths_per_separatfile[fasta_file_name][length] = count
            
        
        #print("Both dict and reversedict: ")
        #print()
        #print(result) #form: length_of_sequence : abundance/count
        #print()
        #print(resultCount) #form: abundance/count : length_of_sequence
        #print()
    
    def get_unique_lengths(self):
        return len(self.lengths_per_separatfile)
    
    def write_to_analysis_lengths(self):
        with open(self.analysis_file_path, 'w') as analysis_file:
            for fasta in self.lengths_per_separatfile.keys():
                analysis_file.write(fasta + ": \n")
                analysis_file.write("unique length, " + "abundance \n")
                for length, count in self.lengths_per_separatfile[fasta].items():
                    if (len(self.lengths_per_separatfile[fasta]) == 0):
                        analysis_file.write("File does not have sequences")
                    else:
                        analysis_file.write(str(length) + " , " + str(count) + "\n")
       
    
        
    def GenotypeLength(self, R_exec):
        if (self.par_os == "windows"):
            R = os.path.normpath(R_exec)
            
        else:
            R = "Rscript"
        
        list_error = []
        #GenotypeLength(RExDip, OutDir + 'MarkerStatistics/', minNrReads, '0.7', OutDir + 'MarkerPlots/', LengthWindow)
        input_folder = os.path.normpath(self.output_folder_path + '/MarkerStatistics')
        output_folder = os.path.normpath(self.output_folder_path + '/MarkerPlots')
        Rscript = None
        if (self.datatype == "diploid"):
            Rscript = self.scriptRDiploid_path
        elif (self.datatype == "haploid"):
            Rscript = self.scriptRHaploid_path
        # make code
        markerplots_file_path = os.path.normpath(output_folder + '/' + 'markerplots.pdf ')
        markermatrix_file_path = os.path.normpath(output_folder + '/' + 'markermatrix.csv ')
        markerlist_file_path = os.path.normpath(output_folder + '/' + 'markerlist.csv ')
        #code = R + ' --vanilla ' + Rscript + ' ' + input_folder +  ' ' + str(self.minCount) + ' ' + markerplots_file_path + str(self.alpha) + ' ' + markermatrix_file_path + markerlist_file_path
        #code = R + ' --vanilla ' + Rscript + ' ' + input_folder + '/' +  ' ' + str(self.minCount) + ' ' + output_folder + '/' + 'markerplots.pdf ' + str(self.alpha) + ' ' + output_folder + '/' + 'markermatrix.csv ' + output_folder + '/' + 'markerlist.csv ' + ' '.join(str(self.lengthWindow))
        #code = R + ' --vanilla ' + Rscript + ' ' + input_folder +  ' ' + str(self.minCount) + ' ' + output_folder + 'markerplots.pdf ' + str(self.alpha) + ' ' + output_folder + 'markermatrix.csv ' + output_folder + 'markerlist.csv' + ' '.join(str(self.lengthWindow))
        code = R + ' --vanilla ' + Rscript + ' ' + input_folder +  ' ' + str(self.minCount) + ' ' + markerplots_file_path + str(self.alpha) + ' ' + markermatrix_file_path + markerlist_file_path + ' '.join(self.lengthWindow)
        print(code) # print code
        #subprocess.call(code, shell=True) # call code
        os_code = os.system(code)
        #subprocess.call(code, shell=True) # run program
        if os_code != 0:
            print(f"The command exited with error code: {os_code}")
            #list_error.append(os_code)
        #else:
            #return os_code
    
        for error in list_error:
            if (error != 0):
                return error
        
        return 0


if __name__ == "__main__":
    main()