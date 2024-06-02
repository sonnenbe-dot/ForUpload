# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 20:06:33 2023

@author: Sebastian
"""

import os, subprocess, re, sys
from datetime import datetime
from collections import Counter
import csv


def main():
    
    paramter_file_path = os.path.normpath(os.getcwd() + '/parameters.txt')
    if not(os.path.exists(paramter_file_path)):
        raise Exception("There is no paramter file with the name parameters.txt")
    
    with open(paramter_file_path) as parameters:
        for line in parameters:
            if (not line) or line.startswith("#"):
                continue
            if (line.startswith("OutputFolder")):
                output_folder = line.split("=")[1].rstrip("\r\n").strip()
                if (not output_folder):
                    output_folder = "OutputDefault"
                continue
            if (line.startswith("Bin")):
                bin_dir = line.split("=")[1].rstrip("\r\n").strip()
                if (not bin_dir):
                    print("No name for data_dir in parameter file!")
                continue
            if (line.startswith("RawData")):
                raw_data_dir = line.split("=")[1].rstrip("\r\n").strip()
                if (not raw_data_dir):
                    print("No name for data_dir in parameter file!")
                continue
            
            if (line.startswith("PrimerFile")):
                primerfile = line.split("=")[1].rstrip("\r\n").strip()
                if (not primerfile):
                    print("No name for primerfile in parameter file!")
                #param_list_files.append(os.path.normpath(primerfile))
                continue
            if (line.startswith("SampleFile")):
                sample_sheet = line.split("=")[1].rstrip("\r\n").strip()
                if (not sample_sheet):
                    print("No name for sample file in parameter file!")
                #param_list_files.append(os.path.normpath(sample_sheet))
                continue
            if (line.startswith("RExecutable")):
                Rexecutable = line.split("=")[1].rstrip("\r\n").strip()
                if (not Rexecutable):
                    Rexecutable = "R"
                #param_list_files.append(os.path.normpath(Rexecutable))
                continue
            if (line.startswith("AlleleList")):
                allele_list = line.split("=")[1].rstrip("\r\n").strip()
                #param_list_files.append(os.path.normpath(allele_list))
                continue
        
            
            
            if (line.startswith("MaxMismatch")):
                max_mismatches = line.split("=")[1].rstrip("\r\n").strip()
                if (not max_mismatches):
                    max_mismatches = 2
                else:
                    max_mismatches = int(max_mismatches)
                #calc_param_list.append(max_mismatches)
                continue
            if (line.startswith("MinCount")):
                minCount = line.split("=")[1].rstrip("\r\n").strip()
                if (not minCount):
                    minCount = 10
                else:
                    minCount = int(minCount)
                #calc_param_list.append(minCount)
                continue
            if (line.startswith("MinLength")):
                min_length = line.split("=")[1].rstrip("\r\n").strip()
                if (not min_length):
                    min_length = 290
                else:
                    min_length = int(min_length)
                #calc_param_list.append(min_length)
                continue
            if (line.startswith("ConsensusThreshold")):
                alpha = line.split("=")[1].rstrip("\r\n").strip()
                if (not alpha):
                    alpha = 0.7
                else:
                    alpha = float(alpha)
                #calc_param_list.append(alpha)
                continue
            if (line.startswith("LengthWindow")):
                values = line.split("=")[1].rstrip("\r\n").strip().split(",")
                if (not values):
                    lengthWindow = ['310', '600']
                else:
                    lengthWindow = []
                    lengthWindow.append(values[0])
                    lengthWindow.append(values[1])
                #calc_param_list.append(lengthWindow)
                continue
            
            if (line.startswith("Ploidy")):
                ploidy = line.split("=")[1].rstrip("\r\n").strip()
                if (not ploidy):
                    ploidy = "diploid"
                #param_list_special.append(ploidy)
                continue
            if (line.startswith("Operatingsystem")):
                par_os = line.split("=")[1].rstrip("\r\n").strip()
                if (not par_os):
                    par_os = "windows"
                #param_list_special.append(par_os.lower())
                continue
            if (line.startswith("Uniqueidentifier")):
                unique_name = line.split("=")[1].rstrip("\r\n").strip()
                if (not unique_name):
                    unique_name = "Default"
                #param_list_special.append(unique_name)
                continue
            if (line.startswith("Indexcomboposition")):
                index_combo = line.split("=")[1].rstrip("\r\n").strip()
                if (not index_combo):
                    index_combo = 1
                else:
                    index_combo = int(index_combo)
                #param_list_special.append(index_combo)
                continue
            
            
    
    print(raw_data_dir)
    print(primerfile)
    print(sample_sheet)
    print(output_folder)
    print(ploidy)
    
    print(max_mismatches)
    print(min_length)
    print(minCount)
    print(alpha)
    print(lengthWindow)
    
    print(allele_list)
    print(par_os)
    print(unique_name)
    
    
    print()

    print("AlleList: " + allele_list)
    
    
    #allele_list = ""
    textbox_checkbox = True
    
    
    
    Sequence_inst = class_quality(output_folder, primerfile, sample_sheet, minCount, alpha, ploidy, unique_name, par_os, allele_list, textbox_checkbox)
    
    #Sequence_inst.print_input_paths()
    print()
    allele_length_dict, number, number_filled, number_empty, samples_without_lengths, samples, sample_duplicates = Sequence_inst.extract_AlleleLengths_haploid()
    print()
    
    number, number_filled, number_empty = Sequence_inst.RunConsensusAll()
    print()
    #loci_number, loci = Sequence_inst.extract_marker_diploid()
    number, number_filled, number_empty = Sequence_inst.joinSamplesSameMarker()
    print()
    number, number_filled, number_empty = Sequence_inst.correctAllSeqs()
    print()
    #sample_number, samples = Sequence_inst.extract_samples()
    additional_loci, markers_more_than2 = Sequence_inst.AlleleCall()
    
    return 0





class class_quality():
    def __init__(self, output_folder_path, primerfile, samplesheet, minCount, consThreashold, datatype, unique_name, par_os, allele_list_path, textbox_bool):
        self.workspace_path = os.getcwd()
        
        
        self.output_folder_path = ""
        output_folder_path = os.path.normpath(output_folder_path)
        if (output_folder_path == "" or not(os.path.exists(output_folder_path)) or output_folder_path == "."):
            self.output_folder_path = os.path.normpath(self.workspacepath + "/DefaultOutput")
        else:
            self.output_folder_path  = output_folder_path 
        
        
        new_workspace_folder_path = os.path.normpath(self.output_folder_path)
        print(new_workspace_folder_path)
        if os.path.exists(new_workspace_folder_path):
            self.new_workspace_folder_path = new_workspace_folder_path
            self.new_workspace_folder = os.path.basename(self.output_folder_path)
        else:
            raise Exception("There is no output folder " + new_workspace_folder_path + ". 2nd part cannot start!")
        
        separatout_path = os.path.normpath(self.new_workspace_folder_path + '/SeparatOut')
        if os.path.exists(separatout_path):
            self.separatout_path = separatout_path
        else:
            raise Exception("There is no output folder " + separatout_path + ". 2nd part cannot start!")
            
        
        matrix_file_path = os.path.normpath(self.new_workspace_folder_path + '/MarkerPlots/markermatrix.csv')
        if os.path.exists(matrix_file_path):
            self.matrix_file_path = matrix_file_path
            self.matrix_file = 'markermatrix.csv'
        else:
            raise Exception("There is no output matrix file (csv) " + self.matrix_file + ". 2nd part cannot start!")
        
        primerfile_path = os.path.normpath(primerfile)
        if os.path.exists(primerfile_path):
            self.primerfile_path = primerfile_path
            self.primerfile = os.path.basename(os.path.normpath(primerfile_path))
            #print("Primerfile: " + self.primerfile + "\n")
        else:
            raise Exception("There is no Primerfile" + self.primerfile + " in the workspace " + self.old_workspace + ". \n")
        
        samplesheet_path = os.path.normpath(samplesheet)
        if os.path.exists(samplesheet_path):
            self.samplesheet_path = samplesheet_path
            self.samplesheet = os.path.basename(os.path.normpath(samplesheet_path))
            #print("Primerfile: " + self.primerfile + "\n")
        else:
            raise Exception("There is no Primerfile" + self.samplesheet + " in the workspace " + self.old_workspace + ". \n")
        
        
        for folder in ['AllelesOut', 'ConsensusOut', 'ConsensusTogether', 'Corrected', 'AlleleCall', 'AdditionalInfo']:
            try:
                os.mkdir(self.new_workspace_folder_path + '/' + folder)
            except OSError:
                print ("Creation of the directory %s failed because it's already there." % folder) # if it fail it produces this error message
            else:
                print ("Successfully created the directory %s " % folder)
        
        
        self.unique_seq_file = 'unique_sequences.txt'
        self.unique_seq_file_path = os.path.normpath(self.new_workspace_folder_path + '/AdditionalInfo/' + self.unique_seq_file)
        
        self.denovo_bool = 0
        #allele_list_path = allele_list_path
        #allele_list_path = os.path.normpath(self.new_workspace + '/' + allele_list)
        if (not(allele_list_path == "") and os.path.exists(os.path.normpath(allele_list_path))):
            print("Hei")
            self.denovo_bool = 1
            self.allele_list_path = os.path.normpath(allele_list_path)
            self.allele_list = os.path.basename(os.path.normpath(allele_list_path))
        else:
            print("There is no Allelelist file path: \n" + os.path.normpath(allele_list_path) + " \n")
            print("Denovo Call!")
            self.allele_list = ""
            self.allele_list_path = ""
            
        
        self.minCount = minCount
        self.datatype = datatype
        self.consThreashold = consThreashold
        
        self.loci = []
        self.loci_number = 0
        self.samples = []
        self.sample_number = 0
        
        self.unique_name = unique_name
        
        self.analysis = {}
        
        self.textbox_bool = textbox_bool
        
        
        self.files_addedN = []
        self.additional_loci = []
    
    def run_script(self):
        print()
        self.print_input_paths()
        print()
        self.get_Marker_Sample_Alleles_length_diploid()
        print()
        print(self.sample_names)
        self.RunConsensusAll()
        print()
        self.joinSamplesSameMarker()
        print()
        self.correctAllSeqs()
        print()
        self.AlleleCall()
        print()
    
    def print_input_paths(self): #printing out all the nessecary folders and files:
        print()
        print("Workspace path: " + self.new_workspace_folder_path + "\n")
        print("Separatout path: " + self.separatout_path + "\n")
        print("Matrixfile (csv) path: " + self.matrix_file_path + "\n")
        print("Primerfile path: " + self.primerfile_path + "\n")
        
        
    #This is the method that is used to parse the file with the location given by file_path
    #The file has to be of fasta or fastq type
    def parsing_fastq_fasta(self, file_path): 
        with open(file_path) as input_file:
            data = {}
            filetype = os.path.basename(os.path.normpath(file_path)).split('.')[1]
            if (filetype == 'fasta'):
                for line in input_file:
                    if line.startswith('>'): # head starts with '>'
                        head = line.rstrip('\r\n')[1:] # exclude new line and starting character '>'
                        if head in data.keys():
                            head += '_2' # if head already existing in dictionary add _2 sufix so seuqneces with the same head are considered
                        data[head] = '' #The value for each head is an empty string, this string will be a sequence in the next step
                    else:
                        if (head in data.keys()):
                            data[head] = line.rstrip('\r\n').upper()
                        
            elif (filetype == 'fastq'):
                for i, line in enumerate(input_file):
                    if i % 4 == 0: # headers line numebers when devide by four should always have a remainder of 0
                        name = line.rstrip('\r\n')[1:] # exclude new line and starting character '@'
                    elif i % 4 == 1: # sequences line numeber when devide by four should always have a remainder of 1
                        data[name] = line.rstrip('\r\n')
            else:
                raise Exception("No fasta or fastq. Pipeline can't continue. \n")
        
        return data
    
    
    def getFileType(self, folder_path):
        folder_name = os.path.basename(os.path.normpath(folder_path))
        for file in os.listdir(folder_path):
            filetype = file.split('.')[-1] # get extension
            if (filetype == 'fasta' or filetype == 'fastq'):
                return filetype
        raise Exception("No fasta or fastq in " + folder_name + "directory. Pipeline can't continue. \n")
        
    
    #filepath is here each fasta/fastq file from the seperatoutfolder, we extract the sequences that correspond to length parameter and has no Ns
    def extract_sequences_per_length(self, file_path, length):
        result = {}
        for rec, seq in self.parsing_fastq_fasta(file_path).items():
            if length == len(seq) and seq.count("N") == 0:
                result[rec] = seq
        return result
    
    def writeFasta(self, SeqDict, out):
        with open(out,'w') as outFasta:
            for rec, seq in SeqDict.items():
                outFasta.write('>' + rec[1:] + '\n' + seq + '\n')
                
    
            
    def get_Marker_Sample_Alleles_length_haploid(self):
        folder_path = self.separatout_path
        separatout_filetype = self.getFileType(folder_path) #Get the filetype of separatout files (should be fastq or fasta)
    
    
    
    
    
    def remove_duplicates(self, x):
        uniques = set()
        result = []
        for element in x:
            if element not in uniques:
                uniques.add(element)
                result.append(element)
        return result
    
    
    def extract_marker_diploid(self):
        with open(self.matrix_file_path) as allele_matrix:
            for line in allele_matrix:
                line_splitted = re.split('\t|,|;', line.rstrip('\r\n').replace('"', '')) # parse matrix lines
                sample = line_splitted[0].strip('"')
                if (sample == 'samplename'): #we are in the first line which is the header
                    loci = line_splitted[1:]
                    loci = ['_'.join(locus.split('_')[0:2]) for locus in loci] #trimming away the _1 and _2 parts of each loci name
                    self.loci = self.remove_duplicates(loci) #remove duplicates while keeping the same order!
                    self.loci_number = len(self.loci)
                    self.analysis["Number of loci"] = str(len(self.loci))
                    break
                
        return self.loci_number, self.loci
    
    
    def extract_marker_haploid(self):
        with open(self.matrix_file_path) as allele_matrix:
            for line in allele_matrix:
                line_splitted = re.split('\t|,|;', line.rstrip('\r\n').replace('"', '')) # parse matrix lines
                sample = line_splitted[0].strip('"')
                if (sample == 'samplename'): #we are in the first line which is the header
                    loci = line_splitted[1:]
                    loci = ['_'.join(locus.split('_')[0:2]) for locus in loci] #trimming away the _1 and _2 parts of each loci name
                    #self.loci = self.remove_duplicates(loci) #remove duplicates while keeping the same order!
                    self.loci_number = len(self.loci)
                    self.analysis["Number of loci"] = str(len(self.loci))
                    break
                
        return self.loci_number, self.loci
    
    
    def extract_samples(self):
        with open(self.matrix_file_path) as allele_matrix:
            for line in allele_matrix:
                line_splitted = re.split('\t|,|;', line.rstrip('\r\n').replace('"', '')) # parse matrix lines
                sample = line_splitted[0].strip('"')
                bool_sample_without_length = True
                if (not(sample == 'samplename')): #we are in the first line which is the header
                    #print(sample)
                    self.samples.append(sample)
        
        self.samples = self.remove_duplicates(self.samples)
        
        return len(self.samples), self.samples
    
    def extract_loci(self):
        with open(self.matrix_file_path) as allele_matrix:
            for line in allele_matrix:
                line_splitted = re.split('\t|,|;', line.rstrip('\r\n').replace('"', ''))
                sample = line_splitted[0].strip('"')
                bool_sample_without_length = True
                if (sample == 'samplename'): #we are in the first line which is the header
                    loci = line_splitted[1:]
                    loci = ['_'.join(locus.split('_')[0:2]) for locus in loci] #trimming away the _1 and _2 parts of each loci name
                    self.loci = self.remove_duplicates(loci) #remove duplicates while keeping the same order!
                    self.loci_number = len(self.loci)
                    self.analysis["Number of loci"] = str(len(self.loci))
                    break
        
        
    
    def extract_AlleleLengths_diploid(self):
        folder_path = self.separatout_path
        separatout_filetype = self.getFileType(folder_path) #Get the filetype of separatout files (should be fastq or fasta)
        
        self.loci = []
        self.samples = []
        allele_length_dict = {}
        allele_lengths_per_locus = {}
        allele_lengths_per_sample = {}
        
        sample_duplicates = []
        
        samples_without_lengths = []
        
        with open(self.matrix_file_path) as allele_matrix:
            for line in allele_matrix:
                line_splitted = re.split('\t|,|;', line.rstrip('\r\n').replace('"', '')) # parse matrix lines
                sample = line_splitted[0].strip('"')
                bool_sample_without_length = True
                if (sample == 'samplename'): #we are in the first line which is the header
                    loci = line_splitted[1:]
                    loci = ['_'.join(locus.split('_')[0:2]) for locus in loci] #trimming away the _1 and _2 parts of each loci name
                    self.loci = self.remove_duplicates(loci) #remove duplicates while keeping the same order!
                    self.loci_number = len(self.loci)
                    self.analysis["Number of loci"] = str(len(self.loci))
                    #print(loci_names)
                else:
                    if (sample in self.samples):
                        sample_duplicates.append(sample)
                    self.samples.append(sample)
                    allele_lengths = line_splitted[1:]
                    allele_lengths = [allele_length.strip('"').replace('empty', '0').replace('too little reads', '0').replace('NA', '0').split('_')[0] for allele_length in allele_lengths] #Using split('_')[0] to also eliminate the _man check string
                    #allele_lengths contains all the lengths (each length corresponds to a loci (1st row in markermatrix)) and per sample (1st column in markermatrix)
                     
                    if not(len(self.loci) == 0):
                        i = 0
                        for locus in self.loci:
                            allele_length_dict[str(sample)+ "_" + str(locus)] = {}
                            #locus = [i]
                            al1 = allele_lengths[i] # get allele 1 (most reads)
                            al2 = allele_lengths[i + 1] # get allele 2 (less reads)
                            print()
                            print('Processing sample ' + sample + ' for locus ' + locus + '\n')
                            print("Lengths: \n")
                            print(al1)
                            print(al2)
                            print()
                            if al1 != '0' and al2 != '0': # missing data (0) is not considered
                                print('......genotype: ' + al1 + ',' + al2)
                                if separatout_filetype == 'fasta':
                                    separatout_file_path = os.path.normpath(self.separatout_path + '/' + sample + '_' + locus + '.fasta')
                                elif separatout_filetype == 'fastq':
                                    separatout_file_path = os.path.normpath(self.separatout_path + '/' + sample + '_' + locus + '.fastq')
                                else:
                                    raise Exception("No fasta or fastq in SeparatOut directory. Pipeline can't continue. \n")
                                
                                #For each loci/sample combination we extract the sequences for the given allele length
                                for al in list(set([al1, al2])):
                                    outFasta = os.path.normpath(self.new_workspace_folder_path + '/AllelesOut/' + locus + '_' + sample +  '_Al_' + al + '.fasta')
                                    extracted_sequences = self.extract_sequences_per_length(separatout_file_path, int(al))
                                    #extracted sequences is a dictionary with the headers (the lines starting with '>' and keys being the sequences)
                                    
                                    if (len(extracted_sequences) == 0):
                                        print('No sequence found for allele length ' + al)
                                        allele_length_dict[str(sample)+ "_" + str(locus)][al] = 0
                                    else:
                                        bool_sample_without_length = False
                                        print(str(len(extracted_sequences)) + " (yes) sequences found for Sample " + sample + ", locus " + locus + " allele length " + al)
                                        allele_length_dict[str(sample)+ "_" + str(locus)][al] = len(extracted_sequences)
                                        self.writeFasta(extracted_sequences, outFasta) #  write to file
                            else:
                                print("....genotype missing.")
                            
                            i += 2
                    
                    else:
                        raise Exception("No Loci! \n")
                    
                    if (bool_sample_without_length):
                        samples_without_lengths.append(sample)
        
        number, number_filled, number_empty = self.check_results(os.path.normpath(self.new_workspace_folder_path + '/AllelesOut'))
         
        return allele_length_dict, number, number_filled, number_empty, list(set(samples_without_lengths)), self.samples, list(set(sample_duplicates))
    
    
    def extract_AlleleLengths_haploid(self):
        folder_path = self.separatout_path
        separatout_filetype = self.getFileType(folder_path) #Get the filetype of separatout files (should be fastq or fasta)
        
        self.loci = []
        self.samples = []
        allele_length_dict = {}
        allele_lengths_per_locus = {}
        allele_lengths_per_sample = {}
        
        sample_duplicates = []
        
        samples_without_lengths = []
        
        with open(self.matrix_file_path) as allele_matrix:
            for line in allele_matrix:
                line_splitted = re.split('\t|,|;', line.rstrip('\r\n').replace('"', '')) # parse matrix lines
                sample = line_splitted[0].strip('"')
                bool_sample_without_length = True
                if (sample == 'samplename'): #we are in the first line which is the header
                    loci = line_splitted[1:]
                    loci = ['_'.join(locus.split('_')[0:2]) for locus in loci] #trimming away the _1 and _2 parts of each loci name
                    self.loci = self.remove_duplicates(loci) #remove duplicates while keeping the same order!
                    self.loci_number = len(self.loci)
                    self.analysis["Number of loci"] = str(len(self.loci))
                    #print(loci_names)
                else:
                    if (sample in self.samples):
                        sample_duplicates.append(sample)
                    self.samples.append(sample)
                    allele_lengths = line_splitted[1:]
                    allele_lengths = [allele_length.strip('"').replace('empty', '0').replace('too little reads', '0').replace('NA', '0').split('_')[0] for allele_length in allele_lengths] #Using split('_')[0] to also eliminate the _man check string
                    #allele_lengths contains all the lengths (each length corresponds to a loci (1st row in markermatrix)) and per sample (1st column in markermatrix)
                     
                    if not(len(self.loci) == 0):
                        i = 0
                        for locus in self.loci:
                            allele_length_dict[str(sample)+ "_" + str(locus)] = {}
                            #locus = [i]
                            al = allele_lengths[i] # get allele 1 (most reads)
                            print()
                            print('Processing sample ' + sample + ' for locus ' + locus + '\n')
                            print("Lengths: \n")
                            print(al)

                            print()
                            if (al != '0'): # missing data (0) is not considered
                                print('......genotype: ' + al + ',')
                                if separatout_filetype == 'fasta':
                                    separatout_file_path = os.path.normpath(self.separatout_path + '/' + sample + '_' + locus + '.fasta')
                                elif separatout_filetype == 'fastq':
                                    separatout_file_path = os.path.normpath(self.separatout_path + '/' + sample + '_' + locus + '.fastq')
                                else:
                                    raise Exception("No fasta or fastq in SeparatOut directory. Pipeline can't continue. \n")
                                
                                #For each loci/sample combination we extract the sequences for the given allele length
                                
                                outFasta = os.path.normpath(self.new_workspace_folder_path + '/AllelesOut/' + locus + '_' + sample +  '_Al_' + al + '.fasta')
                                extracted_sequences = self.extract_sequences_per_length(separatout_file_path, int(al))
                                #extracted sequences is a dictionary with the headers (the lines starting with '>' and keys being the sequences)
                                    
                                if (len(extracted_sequences) == 0):
                                    print('No sequence found for allele length ' + al)
                                    allele_length_dict[str(sample)+ "_" + str(locus)][al] = 0
                                else:
                                    bool_sample_without_length = False
                                    print(str(len(extracted_sequences)) + " (yes) sequences found for Sample " + sample + ", locus " + locus + " allele length " + al)
                                    allele_length_dict[str(sample)+ "_" + str(locus)][al] = len(extracted_sequences)
                                    self.writeFasta(extracted_sequences, outFasta) #  write to file
                            else:
                                print("....genotype missing.")
                            
                            i += 1
                    
                    else:
                        raise Exception("No Loci! \n")
                    
                    if (bool_sample_without_length):
                        samples_without_lengths.append(sample)
        
        number, number_filled, number_empty = self.check_results(os.path.normpath(self.new_workspace_folder_path + '/AllelesOut'))
         
        return allele_length_dict, number, number_filled, number_empty, list(set(samples_without_lengths)), self.samples, list(set(sample_duplicates))
    
    
    def check_results(self, folder_path):
        print()
        number = 0
        number_empty = 0
        number_filled = 0
        if os.path.exists(folder_path):
            for root, dirs, files in os.walk(folder_path):
                number += len(files)
                for file in files:
                    file_path = os.path.join(root, file)
                    if (not(os.path.getsize(file_path) == 0)):
                        number_filled += 1
                    else:
                        number_empty += 1
        return number, number_filled, number_empty
    
    #The file contains all sequences for a sample-locus-allele_length combination
    def MakeConsensusPerFile(self, file_path): #the function works on one file from the allelesOut folder
        outDir = os.path.normpath(self.new_workspace_folder_path + '/ConsensusOut')
        
        fasta_dict = self.parsing_fastq_fasta(file_path) #this dictionary is from one fasta/fastq file from AllelesOut folder
        #keys are headers and values are sequences
        nr_seqs = len(fasta_dict) # number of sequences per file
        file_name = os.path.basename(os.path.normpath(file_path))
        
        # get sequence length. This is necessary to iterate through seuqence positions
        SeqLength_fixed = 0
        i = 0
        #We check if each sequence really has the same length!
        for header, sequence in fasta_dict.items():
            if (i == 0):
                SeqLength_fixed = len(sequence)
            if (i > 0):
                SeqLength_new = len(sequence)
                if (SeqLength_new != SeqLength_fixed):
                    raise Exception("Different Sequence Lengths in the file " + file_name + " \n")
            i += 1
            
        # creates a dictionary with nucleotide composition per alignemt collumn {position : [nucleotides for all sequences]}
        data = {}
        for i in range(0, SeqLength_fixed): # iterate through positions
            nuc = [] #list of all nucleotides for each sequence in the file at a certain position i
            for header, sequence in fasta_dict.items():
                nuc.append(sequence[i]) # save nucleotides for all sequences into a list
            data[i] = nuc
            
        outName = file_name.split('.')[0] + '_C_' + str(nr_seqs)
        con_sequence_file_path = os.path.normpath(outDir + '/' + outName + '_' + str(int(float(self.consThreashold)*100)) + '.fasta')
        frequency_file_path = os.path.normpath(outDir + '/' + outName + '.txt')
        
        with open (con_sequence_file_path, 'w') as cons_file:
            with open(frequency_file_path, 'w') as freq_file:
                conSeq = '>' + outName + '\n' # consensus sequence head
                freq_file.write('Position' + '\t' + 'Freq(A)' + '\t' + 'Freq(C)' + '\t' + 'Freq(G)' + '\t' + 'Freq(T)' + '\t' + '\n')
                for position, nuc in data.items(): #nuc is a list of nucleotides for the position rec for all sequences in the file
                    temp = str(position + 1) # position in sequence
                    # it assumes N as the defoult nucleotide and only replaces it if a nucleotyde has a frequency above the threashold
                    base = 'N' # define nucleotide as N
                    for n in ['A', 'C', 'G', 'T']:
                        freq = float(nuc.count(n))/len(nuc) # get frequency from each nucleotide #nr_seq should be the same value as len(nuc)
                        temp += '\t' + str(freq) # save fequency in freq file
                        #check of nucleotide frequncy is above the threashold
                        if freq >= float(self.consThreashold):
                            base = n 
                        #else:
                            #self.files_addedN.append(os.path.basename(os.path.normpath(con_sequence_file_path)))
                    freq_file.write(temp + '\n')
                    conSeq += base # save nucleotide in consensus sequence
                cons_file.write(conSeq + '\n')
                
    
    
    def RunConsensusAll(self):
        # just run the consensus function for all files
        allelesout_folder_path = os.path.normpath(self.new_workspace_folder_path + '/AllelesOut') #the input folder
        
        #self.files_addedN = []
        
        for file in os.listdir(allelesout_folder_path):
            print('making consensus of ' + file)
            filepath = os.path.normpath(allelesout_folder_path + '/' + file)
            self.MakeConsensusPerFile(filepath)
        
        number, number_filled, number_empty = self.check_results(os.path.normpath(self.new_workspace_folder_path + '/ConsensusOut'))
         
        return number, number_filled, number_empty
            
    def joinSamplesSameMarker(self):
        input_path = os.path.normpath(self.new_workspace_folder_path + '/ConsensusOut')
        output_path = os.path.normpath(self.new_workspace_folder_path + '/ConsensusTogether')
        
        for locus in self.loci:
            print("Locus:")
            print(locus)
            print()
            output_fasta_path = os.path.normpath(output_path + '/' + locus + '_together.fasta')
            with open(output_fasta_path, 'w') as output_fasta:
                for file in os.listdir(input_path):
                    if file.startswith(locus+"_") and file.endswith('.fasta'):
                        input_fasta_path = os.path.normpath(input_path + '/' + file)
                        with open(input_fasta_path) as input_fasta:
                            for line in input_fasta:
                                output_fasta.write(line)
        
        
        number, number_filled, number_empty = self.check_results(os.path.normpath(self.new_workspace_folder_path + '/ConsensusTogether'))
                                
        return number, number_filled, number_empty
    
    
    def correctAllSeqs(self):
        consensustogether_path = os.path.normpath(self.new_workspace_folder_path + '/ConsensusTogether')
        allelesout_path = os.path.normpath(self.new_workspace_folder_path + '/AllelesOut')
        OutDir = os.path.normpath(self.new_workspace_folder_path + '/Corrected')
        
        for file in os.listdir(consensustogether_path):
            print('correcting file ' + file)
            inFasta_path = os.path.normpath(consensustogether_path + '/' + file)
            outFasta_path = os.path.normpath(OutDir + '/' + file.split('.')[0] + '_Corr.fasta')
            if (self.datatype.upper() == "DIPLOID"):
                self.correctSequences(inFasta_path, allelesout_path, outFasta_path)
            elif (self.datatype.upper() == "HAPLOID"):
                self.correctSequences_Hap(inFasta_path, allelesout_path, outFasta_path)
        
        number, number_filled, number_empty = self.check_results(os.path.normpath(self.new_workspace_folder_path + '/Corrected'))
                                
        return number, number_filled, number_empty


    def get_Ns_records(self, sequence):
        result = []
        c = 0 # get position
        #iterate through sequence
        for nucleotide in sequence:
            if nucleotide.upper() == 'N':
                result.append(c) # add position with N to list
            c += 1
        result.sort()
        return result
    
    
    def get_seq_freq(self, fasta, pos_list): #fasta is a fasta file from AllelesOut folder
        haps = [] # list to save all linked nucleotides from the Ns positions
        for name, seq in self.parsing_fastq_fasta(fasta).items():
            # make linked nucleotide sequence (haplotype) and add it to list
            hap = ''
            for i in pos_list:
                hap += seq[i]
            haps.append(hap)
            # use Counter() function to make a dictionary of linked nucleotides frequency: {linked nucleotides: frequncy}
            #print("My haps:")
            #print(haps)
        return dict(Counter(haps))
    
    #Calling this function for each dict {nucleotide : number of occurence} and nrElements is either 1 or 2
    def get_to_el_from_dict(self, my_dict, nrElments): 
        # save the haplotype counts into a list
        print("NOW")
        print(my_dict)
        print()
        temp_count = []
        for hap, c in my_dict.items():
            temp_count.append(c)
        # sort count list
        temp_count.sort()
        print("NOW2")
        print(temp_count)
        print()
        # get most frequent counts
        temp_count_new = temp_count[-nrElments:]
        # get corresponding haplotypes and their relative frequency compared to the others
        result = [] # list to save included haplotypes
        incHapCount = 0 # count summ of the included haplotypes
        excHapCount = 0 # count summ of the excluded haplotypes
        for hap, c in my_dict.items():
            # get haplotypes corresponding to counts and save them in list
            if c in temp_count_new:
                result.append(hap)
                incHapCount += c
            else:
                excHapCount += c
        print('included: ' + str(incHapCount) + ' excluded: ' + str(excHapCount) + ' nr. haps: ' + str(len(result)))
        # Internal control #print('included: ' + str(incHapCount) + ' excluded: ' + str(excHapCount) + ' nr. haps: ' + str(len(result)))
        # decide if haplotypes should be considered. Include them if number of haplotypes are equal to number of elements and if their frequency is above than the remaining
        # sequences without Ns will retrieve empty result list and will go to the else condiction
        if len(result) == nrElments and incHapCount > excHapCount:
            return result
        else:
            return ['N' * len(result[0])] # save Ns has haplotypes in case the conditions are not matached
            # for sequences without Ns. since len(result[0]) == 0, the returned result will be ['']
    
    
    # def correctSequences_Hap(self, fasta, AlleleDir, outFasta):
    #     print("Not ready yet.")
        
    
    def correctSequences_Hap(self, fasta, AlleleDir, outFasta):
        # get positions with Ns and how many haplotypes should be consider
        positions = {} # dictionary to save positions per sequence
        # names_info = {} # dictionary to save number of alleles to be considered per sample
        for name, seq in self.parsing_fastq_fasta(fasta).items():
            if seq != '':
                # key is a tuple with the name of the sample and the sequence, while the item the list with the postions conataing Ns
                positions[(name, seq)] = self.get_Ns_records(seq)
                # # Check if a sample name was already process it by comparing it wiht the names_info dictionary.
                # sample = name.split('_')[2]
                # if sample in names_info:
                #     # If it alrady exists it can only have a maximum of one haplotype because it means it is heterozygote for length
                #     names_info[sample] = 1
                # else:
                #     # Others are homozygote for length and can have more than one allele
                #     names_info[sample] = 2
    
        # File to save output
        out = open(outFasta, 'w')
        # List of samples marked to be excluded
        Samples_to_exclude = []
        # dictionary to save sequences to be saved
        DataOut = {}
    
        # Correct and save sequences
        for name_seq, position in positions.items():
            # get inpute data:
            name = name_seq[0] # sequence name
            seq = name_seq[1] # sequence
            nr_seqs = int(name.split('_')[-1]) # nr of sequences used in consensus
            sample = name.split('_')[2] # sample name
            fileName = os.path.normpath(AlleleDir + '/' + '_'.join(name.split('_')[:5]) + '.fasta') # length extracted sequences file
            # het_info = names_info[sample] # number of haplotypes that can be recovered
            # get frequences of linked nucleotides (haplotypes)
            counter_dict = self.get_seq_freq(fileName, position)
            # get haplotypes to be used
            haps = self.get_to_el_from_dict(counter_dict, 1)
    
    
            # for internal control #print(counter_dict)
    
            # for internal control #print(name + ' nrHapsAllowed: ' + str(het_info) + ' N positions: ' + str(position) + ' allelesFile: ' + fileName)
    
            # Check if haplotypes conatin Ns. If they do, marke the sample to be exluded (add to Samples_to_exclude)
            if haps[0].count('N') > 0:
    
                Samples_to_exclude.append(sample)
                print(sample + ' excluded for ambiguous SNPs')
    
            # replace the sequences with haplotype information and save them in data out. In case of two haplotypes save two sequences
            out_temp = [] # List to save corrected sequences
            # only consider sequences contructed with a read number above the defined threashold
            if nr_seqs > int(self.minCount):
                # iterate through haplotypes
                # for sequences without Ns haps == [''] and thus they are iterable.
                for hap in haps:
                    temp = '' # variable to save the corect sequence
                    # p is position and i is the original sequence nucleotide
                    for p, i in enumerate(seq, 0):
                        # if sequence has an N repalce with corresponding nucleotide in the haplotype
                        if i == 'N' and p in position:
                            temp += hap[position.index(p)]
                        # otherwise just add original nucleotide
                        else:
                            temp += i
                    # add sequence to list
                    out_temp.append(temp)
                # save list of corrected sequences in dictionary
                DataOut[name] = out_temp
            # if read number used for consenus is too low mark sample to be exluded
            else:
                Samples_to_exclude.append(sample)
                print(sample + ' excluded for too little reads')
    
        # Save all corrected sequences for samples not marked to be excluded
        for name, outSeq in DataOut.items():
            sampleName = name.split('_')[2]
            if sampleName not in Samples_to_exclude:
                # All sequences will be save with '_0' at the end of the name. In case of heterozygotes SNPs the second sequence will be saved with '_1'
                for c, seq in enumerate(outSeq, 0):
                    out.write('>' + name + '_' + str(c) + '\n' + seq + '\n')
        out.close()
    
    
    def correctSequences(self, fasta, AlleleDir, outFasta): #the function is called for each fasta file in the consensus_together folder
        positions = {} # dictionary to save positions per sequence
        names_info = {} # dictionary to save number of alleles to be considered per sample
        for name, seq in self.parsing_fastq_fasta(fasta).items(): #Each fasta file in the consensus_together folder has only 1 header and 1 sequence, because they are all put together.
            if seq != '':
                print("header + sequence in consensustogether file: \n")
                print(name)
                print(seq)
                print("\n")
                #name is the complete header (starts with '>'), seq is the sequence underneath the header
                positions[(name, seq)] = self.get_Ns_records(seq) # key is a tuple with the name of the sample and the sequence, while the value is a list with the postions containg Ns
                
                sample = name.split('_')[2] # Check if a sample name was already process it by comparing it wiht the names_info dictionary.
                if sample in names_info:
                    # If it alrady exists it can only have a maximum of one haplotype because it means it is heterozygote for length
                    names_info[sample] = 1
                else:
                    # Others are homozygote for length and can have more than one allele
                    names_info[sample] = 2
        
        
        print("")
        print(positions)
        
        with open(outFasta, 'w') as output:
            Samples_to_exclude = []
            DataOut = {}
            print("Reached1")
            for name_seq, position in positions.items():
                name = name_seq[0] # sequence name, something like HH2_AC_MEL3_Al_468_C_14 (example)
                seq = name_seq[1] # sequence
                print("Reached here")
                #positions is a list with all Ns for this current sequence
                nr_seqs = int(name.split('_')[-1]) # nr of sequences used in consensus, is last part of the header
                sample = name.split('_')[2] # sample name
                fileName = os.path.normpath(AlleleDir + '/' + '_'.join(name.split('_')[:5]) + '.fasta')
                #filename is the name of the fasta file in the AllelesOut dir
                het_info = names_info[sample] #is either 1 or 2 (because only 1 or 2 alleles)
                
                counter_dict = self.get_seq_freq(fileName, position) #position is list with the elements being all positions of Ns of the fasta file filename (1 fasta file from allelesout folder)
                #print("My counter dict:")
                #print(counter_dict)
                haps = self.get_to_el_from_dict(counter_dict, het_info)
                
                #print("The haps:")
                #print(haps)
                
                if haps[0].count('N') > 0:
                    Samples_to_exclude.append(sample)
                    print(sample + ' excluded for ambiguous SNPs')
                    
                #(Samples_to_exclude)
                
                # replace the sequences with haplotype information and save them in data out. In case of two haplotypes save two sequences
                out_temp = [] # List to save corrected sequences
                # only consider sequences contructed with a read number above the defined threashold
                print("Comparing lengths:")
                print(nr_seqs)
                
                if nr_seqs > int(self.minCount):
                    # iterate through haplotypes
                    # for sequences without Ns haps == [''] and thus they are iterable.
                    for hap in haps:
                        temp = '' # varable to save the corect sequence
                        # p is position and i is the original sequence nucleotide
                        for p, i in enumerate(seq, 0):
                            # if sequence has an N repalce with corresponding nucleotide in the haplotype
                            if i == 'N' and p in position:
                                temp += hap[position.index(p)]
                            # otherwise just add original nucleotide
                            else:
                                temp += i
                        # add sequence to list
                        out_temp.append(temp)
                    # save list of corrected sequences in dictionary
                    DataOut[name] = out_temp
                # if read number used for consenus is too low mark sample to be exluded
                else:
                    Samples_to_exclude.append(sample)
                    print(sample + ' excluded for too little reads')
            
            # Save all corrected sequences for samples not marked to be excluded
            for name, outSeq in DataOut.items():
                print("Filtering")
                print(name)
                sampleName = name.split('_')[2]
                print(sampleName)
                if sampleName not in Samples_to_exclude:
                    # All sequences will be save with '_0' at the end of the name. In case of heterozygotes SNPs the second sequence will be saved with '_1'
                    for c, seq in enumerate(outSeq, 0):
                        output.write('>' + name + '_' + str(c) + '\n' + seq + '\n')    


    def AlleleCall(self):
        corrected_folder_path = os.path.normpath(self.new_workspace_folder_path + '/Corrected')
        allelecall_output_folder_path = os.path.normpath(self.new_workspace_folder_path + '/AlleleCall')
        
        #seq_file = 'unique_sequences.txt'
        #unique_seq_folder_path = os.path.normpath(self.new_workspace_folder_path + '/AdditionalInfo/' + seq_file)
        self.additional_loci = []
        markers_more_than2 = []
        if (self.denovo_bool == 0):
            print("Empty or no true path to allelelist")
            self.AelleleCall_denovo(corrected_folder_path, allelecall_output_folder_path)
        elif (self.denovo_bool == 1):
            print("Call listbased!")
            sample_loci_matrix, markers_more_than2 = self.AlleleCall_listbased(corrected_folder_path, allelecall_output_folder_path)
        else:
            print("What is happening?")
        
        return self.additional_loci, list(set(markers_more_than2))
            
    
    def write_genalex(self, sample_allele_dict):
        #dict_upload
        print("Writing Genalex Output")
        #print()
        #print(sample_allele_dict)
        print()
        #sample_allele_dict["samples"] = sample_allele_dict["samples"].split("\t")
        #print(sample_allele_dict["samples"].split("\t"))
        print()
        with open('genalex_auto', 'w', newline='') as file_csv:
            writer = csv.writer(file_csv, delimiter =',', quoting=csv.QUOTE_MINIMAL)
            writer.writerow("")
            writer.writerow("")
            first_row = ['samples'] + sample_allele_dict["samples"].split("\t")[1:]
            #print(first_row)
            writer.writerow(first_row)
            #saves the matrix the genotypes
            sample_dict = self.get_correct_samples()
            for sample in self.samples:
                next_row = [sample_dict[sample]] + sample_allele_dict[sample_dict[sample]].split("\t")[1:]
                #print("Row:")
                #print(next_row)
                writer.writerow(next_row)
                
            
            
            
    def AelleleCall_denovo(self, corrected_folder_path, allelecall_output_folder_path):
        print("Calling Denovo")
        locusAlleleList = {}
        loci = []
        for file in os.listdir(corrected_folder_path): #iterating through all fasta files from the corrected folder, each file has the name format: Locus_together_Corr
            locus = file.split('_')[0] + '_' + file.split('_')[1] # get locus
            print(locus)
            loci.append(locus)
            file_path = os.path.normpath(corrected_folder_path + '/' + file) #get the path to that file
            parsed_fasta = self.parsing_fastq_fasta(file_path) #parse that file
            #parsed_fasta is a dict with keys the heads (Loci_Sample_Length_C_ ... etc.)
            
            
            print()
            print("Printing the locus: \n")
            print(locus + "\n")
            print()
            
            alleles_dict = self.get_allele_dict(parsed_fasta)
            locusAlleleList[locus] = alleles_dict 
            #Each parsed fasta represents all sequences for 1 locus, thus the allele_dict is the value for 1 corresponding locus
        
        sample_loci_matrix = self.CallAlleles(locusAlleleList, corrected_folder_path, allelecall_output_folder_path)
        
        print()
        #print(locusAlleleList)
        print()
        #print(sample_loci_matrix)
        print()
    
    def AlleleCall_listbased(self, corrected_folder_path, allelecall_output_folder_path):
        print("Calling listbased")
        
        #if not(os.path.exists(self.allele_list_path)):
            #raise Exception("There is no allelelist for the listbased call!")
        
        sample_loci_matrix={}
        markers_more_than2 = []
        
        primer_dict = self.parse_primerfile()
        print(primer_dict)
        print("Hei1")
        AllelesInfo = self.parseAlleleList(primer_dict) #this dict has the form: locus : ( {allele_index : [sequence]}, list of sequences )
        #print(AllelesInfo)
        print("Hei2")
        NewLocusAlleleList = self.Complete_AlleleList_All(corrected_folder_path, AllelesInfo)
        print("Hei3")
        sample_loci_matrix, markers_more_than2 = self.CallAlleles(NewLocusAlleleList, corrected_folder_path, allelecall_output_folder_path)
        
        print()
        #print(NewLocusAlleleList)
        print()
        #print(sample_loci_matrix)
        print("Finished")
        #print(NewLocusAlleleList["TI12_TAC"])
        return sample_loci_matrix, markers_more_than2
        
        
    def Complete_AlleleList_All(self, corrected_folder_path, AllelesInfo):
        locusAlleleList = {}
        for file in os.listdir(corrected_folder_path):
            locus = file.split('_')[0] + '_' + file.split('_')[1] # get locus name
            file_path = os.path.normpath(corrected_folder_path + '/' + file) #get the path to that file
            parsed_fasta = self.parsing_fastq_fasta(file_path)
            if locus not in AllelesInfo: #check if key locus exists in AllelesInfo
                print('locus ' + locus + ' not present in allele list') # report non existing loci
                alleles_dict = self.get_allele_dict(parsed_fasta) # if it does not exist alleles are defined denovo
                locusAlleleList[locus] = alleles_dict
            else:
                alleles_info_locus = AllelesInfo[locus] # get infomation for a partivular locus
                alleles_dict_locus = alleles_info_locus[0]
                #print("alleles_dict_locus")
                #print(locus + "\n")
                #print(alleles_dict_locus) # get alleles dictionary
                seq_allele_locus = alleles_info_locus[1] # get list of already existing seuqences for that locus
                new_alleles_dict_locus = self.Complete_AlleleList_PerLocus(alleles_dict_locus, seq_allele_locus, parsed_fasta) # complete allele information
                locusAlleleList[locus] = new_alleles_dict_locus
            
        for locus in AllelesInfo:
            if locus not in locusAlleleList:
                locusAlleleList[locus] = AllelesInfo[locus][0]
        
        print("Access2:")
        print("Now")
        print(locusAlleleList["TI53_ATAG"])
        print("End")
        
        return locusAlleleList

    
    def Complete_AlleleList_PerLocus(self, alleles_dict_per_locus, seq_allele_List_locus, fasta_parsed):
        # define new lists and dictionraies to save the data with preexisting data
        alleles_dict_per_locus_New = alleles_dict_per_locus
        seq_allele_List_locus_New = seq_allele_List_locus
        # compare new sequences with existing sequence list
        for head, seq in fasta_parsed.items():
            allele_Number = len(alleles_dict_per_locus) # get number of alleles already existing. This is necessary to define new alleles, which will be allele_Number + 1
            #print("Number: " + str(allele_Number) + " \n")
            seqDeGap = seq.replace('-', '') # degap sequence in case of alignement
            # check if sequence exists
            if not seqDeGap in seq_allele_List_locus_New:
                allele_Number += 1 # make new allele name
                alleles_dict_per_locus_New[str(allele_Number)] = [seqDeGap] # add new allele and corresponding sequence to dictionary
                seq_allele_List_locus_New.append(seqDeGap) # add new sequence to sequence list
        return alleles_dict_per_locus_New
        
        
    def parseAlleleList(self, primer_dict):
        result = {}
        loci_in_primers = True
        
        primers_rev = []
        
        with open(self.allele_list_path) as allele_list:
            print("Opening")
            for line in allele_list:
                #print(line)
                led = line.rstrip('\r\n')
                if (len(led.split('_')) == 2): # This identify the marker name. All markers names have an underscore, thus if you slit it it will return a list of two ellements
                    locus = '_'.join(led.split('_')[:2])
                    #print()
                    #print("locus is " + locus + "\n")
                    alleles = {}
                    seq_allele_List = []
                    if (locus in primer_dict):
                        primers = primer_dict[locus]
                        primerF = primers[0]
                        primerR = self.rev_comp(primers[1])
                        primers_rev.append(primerR)
                        loci_in_primers = True
                    else:
                        self.additional_loci.append(locus)
                        loci_in_primers = False
                # get end of locus and save data
                elif led == '\\' or led == '' or led == "/": # the end of all loci is defined by a line with \
                    result[locus] = (alleles, seq_allele_List)
                # get allele and sequence information
                else:
                    led = led.split('\t')
                    allele = led[0].rstrip(':') # allele is saved with the format "name:"
                    seq = led[1].replace('-', '') # in case loci were made from aligned data we strip gaps out
                    # do not consider alleles marker for manual control. Tis means that they have an N and match to multiple alleles. These are marked with MC. This only aplies for the first versions of the script where
                    if allele != 'MC':
                        if (loci_in_primers):
                            #new_seq = primerF + seq[len(primerF):-len(primerR)] + primerR
                            new_seq = seq
                        else:
                            print("Reached")
                            #print(locus)
                            #print(str(allele))
                            new_seq = seq
                        alleles[allele] = alleles.get(allele, []) + [new_seq]
                        seq_allele_List.append(new_seq)
        
        #result["TI46_TAT"]
        #for locus in self.additional_loci:
            #result.pop(locus)
        print("Access:")
        print("Now")
        print(str(len(result)))
        #print(result["TI12_TAC"])
        print()
        #print(result["TI53_ATAG"])
        #print(result["TI12_TAC"])
        #print(result["TI32_AAAAT"])
        #print(result["TI46_TAT"])
        #print(result["TI45_ATATA"])
        print("Finished parsing AlleleList")
        return result
    
    def parse_primerfile(self):
        primer_dict={}
        with open(self.primerfile_path) as primerfile:
            for line in primerfile:
                line_list = re.split('\t|,|;', line.rstrip('\r\n')) # make list based seperating elemntes based on tabs
                #if list(set(line_list)) != [""]:
                if list(set(line_list)) != [""]:
                    primers = line_list[1:3] # primer sequences
                    primer_dict[line_list[0]]=primers # save {locus_name: [Forward, Reverse]}
        return primer_dict
    
    def get_correct_samples(self):
        sample_dict = {}
        with open(self.samplesheet_path) as samplesheet:
            for line in samplesheet:
                sample_script = re.split(',|\t|;', line.rstrip('\r\n'))[2]
                sample_original = re.split(',|\t|;', line.rstrip('\r\n'))[1]
                sample_dict[sample_script] = sample_original
        
        #print(sample_dict)    
        
        return sample_dict
    
    def CallAlleles(self, locusAlleleList, corrected_folder_path, allelecall_output_folder_path):
        print("Call Alleles")
        self.samples = self.remove_duplicates(self.samples)
        print(self.samples)
        
        sample_dict = self.get_correct_samples()
        
        markers_more_than2 = []
        
        seq_list_unique_entire = []
        out_dict = {'samples' : '\t'}
        for sample in self.samples:
            new_sample = sample_dict[sample]
            out_dict[new_sample] = ''
        
        locus_already_written = []
        #self.loci
        
        file_list = []
        
        for file in os.listdir(corrected_folder_path):
            file_list.append(file)
        
        sorted_filenames = sorted(file_list, key=lambda file: self.loci.index(file.split('_')[0] + '_' + file.split('_')[1]))
        
        output_allelelist_path = os.path.normpath(allelecall_output_folder_path + '/allelle_list.txt')
        output_allelelist_path_fasta = os.path.normpath(allelecall_output_folder_path + '/allelle_list.fasta')
        with open(output_allelelist_path, 'w') as allele_list:
            with open(output_allelelist_path_fasta, 'w') as allele_list_fasta:
                for file in sorted_filenames:
                #for file in os.listdir(corrected_folder_path): #each file represents the information for 1 unique locus
                    #Iterate through all files from the corrected folder
                    locus = file.split('_')[0] + '_' + file.split('_')[1] # get locus
                    locus_already_written.append(locus) ####
                    out_dict['samples'] += (locus + '\t')*2 #saves the current locus in the dictionary
                    print()
                    print("Get Alleles Dict for Locus " + locus +  "\n")
                    alleles_dict = locusAlleleList[locus]
                    
                    #alleles_dict = dict(sorted(alleles_dict.items(), key=int))
                    alleles_dict = dict(sorted(alleles_dict.items(), key=lambda item: int(item[0])))
                    
                    
                    file_path = os.path.normpath(corrected_folder_path + '/' + file)
                    parsed_fasta = self.parsing_fastq_fasta(file_path) #Get all sequences from one fasta file of the corrected folder
                    names_alleles, list_numbers = self.Define_genotypes(alleles_dict, parsed_fasta)
                    #if (locus == "TI12_TAC"):
                        #print("Here is numbers: \n")
                        #print(list_numbers)
                    #names_alleles is a dict with sample names as keys and and a list of values all corresponding alleles (saved through the indices as a list)
                    allele_list.write(locus + '\n')
                    for allele, allele_seq in alleles_dict.items(): #alleles_dict in the form index : sequence, one index can have 2 sequences
                        for seq in allele_seq:
                            seq_list_unique_entire.append(seq)
                            allele_list.write(allele + ':\t' + seq + '\n')
                            allele_list_fasta.write('>' + locus + ':' + allele + '\n')
                            allele_list_fasta.write(seq + '\n')
                    allele_list.write('\\\n')
                    
                    for sample in self.samples:
                        # call only if sample exists. Otherwise the genotype will be 0/t0 for missing data
                        if sample in names_alleles:
                            alleles = names_alleles[sample] #get all allele indices for this sample name
                            alleles.sort() #sorting according to keys which are the sample names
                            # if two alleles saves it as heterozygote
                            if len(alleles) == 2:
                                genotype = '\t'.join(alleles)
                            # if more report the genotype as ambiguous and mark genotype for manual control. It does not aply anymore
                            elif len(alleles) > 2:
                                markers_more_than2.append(locus)
                                print('locus ' + locus + ' for sample ' +  sample + '\thas ' + str(len(alleles)) + ' alleles: ' + str(alleles))
                                genotype = 'MC\tMC'
                            # if one saves as homozygote
                            else:
                                genotype = '\t'.join(alleles*2)
                        else:
                            genotype = '0\t0'
                        new_sample = sample_dict[sample]
                        out_dict[new_sample] += '\t' + genotype #it saves the genotypes in the dictionary
                    
                
                for locus in locusAlleleList:
                    if locus not in locus_already_written:
                        alleles_dict = locusAlleleList[locus]
                        alleles_dict = dict(sorted(alleles_dict.items(), key=lambda item: int(item[0])))
                        allele_list.write(locus + '\n')
                        for allele, allele_seq in alleles_dict.items(): #alleles_dict in the form index : sequence, one index can have 2 sequences
                            for seq in allele_seq:
                                seq_list_unique_entire.append(seq)
                                allele_list.write(allele + ':\t' + seq + '\n')
                                allele_list_fasta.write('>' + locus + ':' + allele + '\n')
                                allele_list_fasta.write(seq + '\n')
                        allele_list.write('\\\n')
                    
        matrix_path = os.path.normpath(allelecall_output_folder_path + '/matrix.txt')
        with open(matrix_path, 'w') as matrix_out:
            matrix_out.write('samples' + out_dict['samples'] + '\n')
            #saves the matrix the genotypes
            for sample in self.samples:
                new_sample = sample_dict[sample]
                matrix_out.write(new_sample + out_dict[new_sample] + '\n')
        
        unique_num = 0 #Here we could call the database and check the number of unique records
        #unique_seq_folder_path
        seq_list_unique_entire = list(set(seq_list_unique_entire))
        
        new_list = []
        
        with open(self.unique_seq_file_path, 'w') as output_seq:
            output_seq.write('Loci' + '\t' + 'index' + '\t' + 'Identifier' + '\t' + 'Sequence' + '\t' + '\n')
            for locus, allele_dict in locusAlleleList.items():
                for num, seq in allele_dict.items():
                    if (not(seq in new_list)):
                        unique_num += 1
                        unique_id = self.unique_name + str(unique_num)
                        output_seq.write(locus + '\t' + num + '\t' + unique_id + '\t' + str(seq[0]) + '\t' + '\n')
                        new_list.append(seq)
                        unique_id = ""
        
        #print("Not written to Genalex")
        #if (self.textbox_bool): #self.textbox_bool
            #self.write_genalex(out_dict)
            #print("Writing to Genalex")
        #else:
            #print("Not written to Genalex")
            
        
        return out_dict, markers_more_than2
        
    
    def get_allele_dict(self, parsed_fasta):
        alleles_dict = {}
        seq_list = [seq for name, seq in parsed_fasta.items()] #get all sequences from that parsed fasta, name is just the header
        # get unque sequences and sort them
        seq_list = list(set(seq_list))
        seq_list.sort()
        alleles_dict = {str(allele): [seq] for (allele, seq) in enumerate(seq_list, 1)} #allele is here a number, starts on 1 and gives each sequence a number starting from 1
        return alleles_dict #
    
    
    def Define_genotypes(self, alleles_dict_per_locus, fasta_parsed):
        result = {}
        list_numbers = []
        for head, seq in fasta_parsed.items():
            number = 0
            name = head.split('_')[2] # get sample name
            #print(name)
            seqDeGap = seq.replace('-', '') # degap sequences in case of aligned files
            for allele, allele_seq in alleles_dict_per_locus.items():
                # compare with sequences
                if (seqDeGap in allele_seq):
                    number += 1
                    #print("Number: \n")
                    #print(str(allele))
                    print()
                    # save information per sample
                    alleles_to_dict = result.get(name, []) + [allele]
                    result[name] = alleles_to_dict
                    #break
            list_numbers.append(number)
            #print(number)
        
            #if (name == "XN23"):
                #print("Get here::: \n")
                #print(result[name])
                #if (len(result[name]) == 3):
                    #print("CR>YYYYYY!!!")
       
        return result, list_numbers
    
    def rev_comp(self, seq):
        # dictionary with Key for complementing nucleotide
    	code = {'A':'T','C':'G','T':'A','G':'C', 'Y':'R', 'R':'Y', 'S':'S', 'W':'W', 'K':'M', 'M':'K', 'B':'V', 'D':'H', 'H':'D', 'V':'B', 'N':'N'}
    	result_sequence = '' # resulting string will be saved here
    	for i in seq: # iterate through nucleotide in input string
    		complement = code.get(i.upper(), i) # get complement. if it does not exist get the input nucleotide
    		result_sequence = complement + result_sequence # add complemented nucleotide in the beggining of the output string. This way it is reversed
    	return result_sequence
    
    
    
    def get_reference_dict(self, sample_loci_matrix):
        
        dict_reference = {}
        os.chdir(self.new_workspace_folder_path)
        
        self.reference_file_path = os.path.normpath("C:/Users/Sebastian/Desktop/Bio_work/local_scripts/Christina_workspace_1/CustomTkinter/test_ws/Presentation/reference.csv")
        
        list_all = []
        list_project = []
        list_organism = []
        list_country = []
        list_locality = []
        
        if os.path.exists(self.reference_file_path):
            with open(self.reference_file_path, "r") as reference:
                csv_reader = csv.reader(reference, delimiter=';')
                for i, line in enumerate(csv_reader):
                    #print(line)
                    sample_ref = line[2]
                    for sample in sample_loci_matrix.keys():
                        if (sample_ref == sample):
                            dict_reference[sample_ref] = {}
                            dict_reference[sample_ref]["Project"] = line[0] #Project
                            list_project.append(line[0])
                            dict_reference[sample_ref]["Organism"] = line[1] #Organism
                            list_organism.append(line[1])
                            dict_reference[sample_ref]["Country"] = line[5] #Country
                            list_country.append(line[5])
                            dict_reference[sample_ref]["Locality"] = line[6] #Locality
                            list_locality.append(line[6])
        
        list_all.append(list_project)
        list_all.append(list_organism)
        list_all.append(list_country)
        list_all.append(list_locality)
        
            
        
        return dict_reference, list_all
    
if __name__ == "__main__":
    main()    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
                
                


    
    
    
    