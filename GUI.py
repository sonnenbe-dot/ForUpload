# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 11:42:37 2023

@author: Sebastian
"""

import tkinter as tk
from tkinter import ttk
import tkinter.messagebox
import customtkinter as ctk
import os 
import re
import csv
import subprocess

from tkinter import filedialog 

from length_class import amplicon_class
from quality_class import class_quality
from database_class import SQdatabase

import threading


ctk.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
ctk.set_default_color_theme("blue")

def main():
    
    print("Testing")
    
    
    root = root_window()
    #frame = ssr_window(root)
    root.mainloop()
    
    
    
    return 0


class root_window(ctk.CTk):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.title("SSRGBAS GUI")
        self.geometry(f"{1000}x{500}")
        
        #Using grid_...configure only for column 0 and row 0 which means my root has only 1 row and 1 column
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        self.grid_propagate(True)
        
        self.main_frame = ssr_window(self)
        self.main_frame.grid(row=0, column=0, sticky="nsew")
        
        # Allow the container to expand in both directions
        self.main_frame.grid_rowconfigure(0, weight=1)
        self.main_frame.grid_columnconfigure(0, weight=1)
        
        
class ssr_window(ctk.CTkFrame):
    def __init__(self, parent):
        super().__init__(parent)
        
        
        self.columnconfigure(0, weight=1) 
        self.columnconfigure(1, weight=0)
        self.columnconfigure(2, weight=1)
        self.columnconfigure(3, weight=0)
        self.columnconfigure(4, weight=1)
        self.columnconfigure(5, weight=0)
        self.columnconfigure(6, weight=0)
        self.rowconfigure(0, weight=1)
        
        
        #self.grid(row=0 , column=0, sticky=tk.NSEW, padx=5, pady=5)
        #self.configure(background='light blue')
        
        
        
        self.parameters = ctk.CTkFrame(self, width=100, corner_radius=2)
        self.parameters.grid(row=0, column=0, rowspan=6, sticky="nsew")
        self.parameters.grid_columnconfigure(0, weight=1)
        #self.parameters.grid_columnconfigure(1, weight=1)
        self.parameters.grid_rowconfigure(5, weight=1)
        
        self.logo_parameters = ctk.CTkLabel(self.parameters, text="Parameters", font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_parameters.grid(row=0, column=0, padx=20, pady=(20, 20))
        #Buttons hell blau eine Spur
        
        self.general_param = ctk.CTkButton(self.parameters, border_width=1, border_color="black", text_color="black", text = "Pipeline Parameters", command = self.get_general_params)
        self.general_param.grid(row=1, column=0, padx=20, pady=(10, 10))
        
        self.genalex_param = ctk.CTkButton(self.parameters, border_width=1, border_color="black", text_color="black", text = "Genalex Parameter", command = self.get_genalex_params)
        self.genalex_param.grid(row=2, column=0, padx=20, pady=(10, 10))
        
        self.sample_import_param = ctk.CTkButton(self.parameters, border_width=1, border_color="black", text_color="black", text = "SampleImport Parameter", command = self.set_sample_import_params)
        self.sample_import_param.grid(row=3, column=0, padx=20, pady=(10, 10))
        
        self.checkit = ctk.CTkButton(self.parameters, border_width=1, border_color="black", text_color="black", text = "Workspace Status", command = self.check_status)
        self.checkit.grid(row=4, column=0, padx=20, pady=(10, 10))
        
        self.doc_param = ctk.CTkButton(self.parameters, border_width=1, border_color="black", text_color="black", text = "Instructions")
        self.doc_param.grid(row=6, column=0, padx=20, pady=(10, 5))
        
        self.appearance = ctk.CTkLabel(self.parameters, text="Appearance:", font=ctk.CTkFont(size=10, weight="bold"))
        self.appearance.grid(row=7, column=0, padx=20, pady=(5, 5))
        
        self.appearance_mode_optionemenu = ctk.CTkOptionMenu(self.parameters, values=["Dark", "Light"], command=self.change_appearance_mode_event)
        self.appearance_mode_optionemenu.grid(row=7, column=0, padx=20, pady=(5, 10))
        
        #self.appearance2 = ctk.CTkLabel(self.parameters, text="Appearance:", font=ctk.CTkFont(size=10, weight="bold"))
        #self.appearance2.grid(row=6, column=0, padx=(20, 20), pady=(5, 10))
        
        
        
        
        
        
        self.separator1 = tk.Canvas(self, width=2, bg='black')
        self.separator1.grid(row=0, column=1, sticky='ns')
        self.separator1.create_line(1, 0, 1, self.separator1.winfo_height())
        
        
        self.pipeline = ctk.CTkFrame(self, width=100, corner_radius=2)
        self.pipeline.grid(row=0, column=2, rowspan=4, sticky="nsew")
        self.pipeline.grid_columnconfigure(0, weight=1)
        self.pipeline.grid_rowconfigure(4, weight=1)
        
        self.logo_pipeline = ctk.CTkLabel(self.pipeline, text="Pipeline", font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_pipeline.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        self.project = ctk.CTkButton(self.pipeline, border_width=1, border_color="black", text_color="black", text = "Choose Project Name", command = self.choose_project)
        self.project.grid(row=1, column=0, padx=20, pady=(10, 40))
        
        self.first_script = ctk.CTkButton(self.pipeline, border_width=1, border_color="black", text_color="black", text = "Run Length Detection", command = self.run_first_script)
        self.first_script.grid(row=2, column=0, padx=20, pady=(10, 10))
        
        self.second_script = ctk.CTkButton(self.pipeline, border_width=1, border_color="black", text_color="black", text = "Run SNP Detection", command = self.run_second_script)
        self.second_script.grid(row=3, column=0, padx=20, pady=(10, 20))
        
        
        
        self.textbox = ctk.CTkTextbox(self.pipeline, width=250)
        self.textbox.grid(row=4, column=0, padx=(5, 5), pady=(5, 5), sticky="nsew")
        self.textbox.insert("0.0", "Pipeline:" + "\n\n" * 2)
        #self.textbox.insert("0.0", "Pipeline:" + "\n\n" * 2)
        
        
        
        self.analysis = ctk.CTkButton(self.pipeline, border_width=1, border_color="black", text_color="black", text = "Get Analysis", command = self.get_analysis)
        self.analysis.grid(row=5, column=0, padx=20, pady=(20, 5))
        
        self.genalex = ctk.CTkButton(self.pipeline, border_width=1, border_color="black", text_color="black", text = "Import as Genalex", command = self.import_genalex_latest)
        self.genalex.grid(row=6, column=0, padx=20, pady=(10, 5))
        
        
        
        
        
        
        
        
        
        
        self.separator2 = tk.Canvas(self, width=2, bg='black')
        self.separator2.grid(row=0, column=3, sticky='ns')
        self.separator2.create_line(1, 0, 1, self.separator1.winfo_height())
        
        
        
        self.database = ctk.CTkFrame(self, width=100, corner_radius=2)
        self.database.grid(row=0, column=4, rowspan=4, sticky="nsew")
        self.database.grid_columnconfigure(0, weight=1)
        self.database.grid_rowconfigure(5, weight=1)
        
        self.logo_database = ctk.CTkLabel(self.database, text="Database", font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_database.grid(row=0, column=0, padx=20, pady=(20, 20), sticky="nsew")
        
        
        self.all = ctk.CTkButton(self.database, border_width=1, border_color="black", text_color="black", text = "Add Dataset", command = self.add_data_to_local_database)
        self.all.grid(row=1, column=0, padx=20, pady=(10, 10))
        
        self.subset = ctk.CTkButton(self.database, border_width=1, border_color="black", text_color="black", text = "Extract Subset", command = self.extract_subset)
        self.subset.grid(row=2, column=0, padx=20, pady=(10, 15))
        
        self.show_subset = ctk.CTkButton(self.database, border_width=1, border_color="black", text_color="black", text = "Show Subset", command = self.show_data)
        self.show_subset.grid(row=3, column=0, padx=20, pady=(10, 15))
        
        self.deleting = ctk.CTkButton(self.database, border_width=1, border_color="black", text_color="black", text = "Delete Dataset", command = self.delete_data_for_folder)
        self.deleting.grid(row=4, column=0, padx=20, pady=(10, 15))
        
        self.textbox2 = ctk.CTkTextbox(self.database, width=250)
        self.textbox2.grid(row=5, column=0, padx=(5, 5), pady=(15, 15), sticky="nsew")
        self.textbox2.insert("0.0", "Database Status:" + "\n")
        
        #self.importing = ctk.CTkButton(self.database, border_width=2, border_color="black", text_color="black", text = "Upload to Boku-Database")
        #self.importing.grid(row=6, column=0, padx=20, pady=(15, 10))
        
        
        
        
        self.projectfolder = "default"
        
        self.sample_loci_matrix = {}
        self.loci_dict = {}
        self.matrix_allele_dict = {}
        
        self.unique_prefix_name = ""
        
        self.var_list_dir = []
        self.var_list_file = []
        self.var_list_calc = []
        self.var_list_special = []
        
        self.param_list_dirs, self.param_list_files, self.param_list_special, self.calc_param_list = self.parse_parameter_file()
        #Order for self.param_list: [rawdata_dir, primerfile, sample_sheet, hyploid_diploid_type, outputfolder, previous_allele_list, operating_system, unique_id]
        #Order for self.calc_param_list: [max_mismatches, min_length, lengthWindow, consThreashold, min_count]
        
        
        self.rawdata_dir = self.param_list_dirs[0]
        self.outputfolder = self.param_list_dirs[1]
        
        self.primerfile = self.param_list_files[0]
        self.sample_sheet = self.param_list_files[1]
        self.allele_list_path = self.param_list_files[2]
        self.reference_file = self.param_list_files[3]
        self.r_script_path = os.path.normpath(self.param_list_files[4])
        
        self.ploidy = self.param_list_special[0]
        self.os_param = self.param_list_special[1]
        if (self.os_param == ""):
            self.os_param = "windows"
        if (self.ploidy == ""):
            self.ploidy = "diploid"
        self.unique_prefix_name = self.param_list_special[2]
        
        self.max_mismatches = self.calc_param_list[0]
        self.minCount = self.calc_param_list[4]
        self.min_length = self.calc_param_list[1]   
        self.alpha = self.calc_param_list[3]
        self.lengthWindow = self.calc_param_list[2]
        
        #self.r_script_path = os.path.normpath("C:\Program Files\R\R-4.3.1\bin\Rscript.exe")
        
        self.is_working = False
        self.is_working2 = False
        
        self.checkbox_states = []
        
        self.checkbox_states_dict = {}
        
        self.checkbox_states_include_dict = {}
        self.checkbox_states_exclude_dict = {}
        
        self.checkbox_states_organism_genalex = {}
        self.checkbox_states_organism_genalex2 = {}
        
        self.checkbox_states_dict2 = {}
        self.checkbox_states_dict2["Organism"] = {}
        self.checkbox_states_dict2["Organism"]["Include"] = []
        self.checkbox_states_dict2["Organism"]["Exclude"] = []
        self.checkbox_states_dict2["Loci"] = {}
        self.checkbox_states_dict2["Loci"]["Include"] = []
        self.checkbox_states_dict2["Loci"]["Exclude"] = []
        self.checkbox_states_dict2["Sample"] = {}
        self.checkbox_states_dict2["Sample"]["Include"] = []
        self.checkbox_states_dict2["Sample"]["Exclude"] = []
        
        self.checkboxes_organism_include = []
        self.values_organism_include = []
        
        self.checkboxes_organism_exclude = []
        self.values_organism_exclude = []
        
        self.checkboxes_loci_include = []
        self.values_loci_include = []
        
        self.checkboxes_loci_exclude = []
        self.values_loci_exclude = []
        
        self.checkboxes_sample_include = []
        self.values_sample_include = []
        
        self.checkboxes_sample_exclude = []
        self.values_sample_exclude = []
        
        self.list_all = []
        
        
        #self.reference_filename = 
        
        self.workspace_location = os.getcwd()
        
        self.projectfolder = "default"
        self.projectfolder_path = os.path.normpath(self.workspace_location + "/" + self.projectfolder)
        
        self.reference_file_name = "reference.csv"
        self.reference_file_path = os.path.normpath(self.workspace_location + "/" + self.reference_file_name)
        
        self.local_database_name = "database_ssr_local.db"
        self.local_database_path = os.path.normpath(self.workspace_location + "/" + self.local_database_name)
        
        self.get_status()
        
        self.sample_loci_matrix = {}
        self.loci_dict = {}
        self.reference_dict = {}
        
        self.textbox_bool = False
        self.import_check = tk.IntVar()
        self.import_check.set(0)
        
        #self.genalex_parameter = ""
        self.genalex_parameter = tk.StringVar()
        
        
        
        self.param_type_mapping = {
            "rawdata_dir": "directory",
            "primerfile": "file",
            "sample_sheet": "file",
            "output_folder": "directory",
            "allelelist path": "file",
            "Reference file" : "file"
        }
        
        self.sample_params_import = ["", "", ""]
        self.sequence_params_import = ["", "", ""]
    
    def check_status(self):
        print("Check if workspace is ready:")
        bool_var = True
        
        if (not(os.path.exists(os.path.normpath(self.workspace_location)))):
            bool_var = False
        if (not(os.path.exists(os.path.normpath(self.workspace_location + '/' + self.rawdata_dir))) or self.rawdata_dir == ""):
            bool_var = False
        if (not(os.path.exists(os.path.normpath(self.workspace_location + '/' + self.primerfile))) or self.primerfile == ""):
            bool_var = False
        if (not(os.path.exists(os.path.normpath(self.workspace_location + '/' + self.sample_sheet))) or self.sample_sheet == ""):
            bool_var = False
        if (not(os.path.exists(os.path.normpath(self.reference_file))) and not(self.reference_file.split(".")[1] == "csv")):
            bool_var = False
        if (not(os.path.exists(os.path.normpath(self.r_script_path))) and not(self.r_script_path.split(".")[0] == "Rscript")):
            bool_var = False
        if (not(os.path.exists(os.path.normpath(self.workspace_location + '/bin')))):
            bool_var = False
        
        
        if (bool_var):
            new_window = tk.Toplevel(self)
            new_window.title("Workspace Status")
            
            # Long text (for demonstration purposes)
            long_text = f"Workspace is ready! \n\n" + "Rawdatafolder : " + self.rawdata_dir + " \n" + "Samplesheet : " + self.sample_sheet + "\n" + "PrimerFile : " + self.primerfile + "\n\n" + "\n" + "Reference file path : " + self.reference_file + "\n" + "R_script exe path : " + self.r_script_path + "\n\n"
            #long_text = f"You selected the file: {filepath}" + "\n\n" + "Lorem ipsum..." * 100
            
            # Add a Text widget with a vertical scrollbar
            text_widget = tk.Text(new_window, wrap=tk.WORD)
            scrollbar = tk.Scrollbar(new_window, command=text_widget.yview)
            text_widget.config(yscrollcommand=scrollbar.set)
            
            text_widget.insert(tk.END, long_text)
            text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=20, pady=20)
            scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        else:
            new_window = tk.Toplevel(self)
            new_window.title("Workspace Status")
            
            # Long text (for demonstration purposes)
            long_text = f"Workspace is NOT ready! Some Folders, Files or values are missing or are not correct. Check below! \n\n" + "Rawdatafolder : " + self.rawdata_dir + " \n" + "Samplesheet : " + self.sample_sheet + "\n" + "PrimerFile : " + self.primerfile + "\n\n" + "\n" + "Reference file path : " + self.reference_file + "\n" + "R_script exe path : " + self.r_script_path + "\n\n"
            #long_text = f"You selected the file: {filepath}" + "\n\n" + "Lorem ipsum..." * 100
            
            # Add a Text widget with a vertical scrollbar
            text_widget = tk.Text(new_window, wrap=tk.WORD)
            scrollbar = tk.Scrollbar(new_window, command=text_widget.yview)
            text_widget.config(yscrollcommand=scrollbar.set)
            
            text_widget.insert(tk.END, long_text)
            text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=20, pady=20)
            scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
            
    
    def set_sample_import_params(self):
        
        print("Setting Sample Import:")
        newWindow = tk.Toplevel(self)
        newWindow.title("Reference Import Parameters")
        newWindow.grid_rowconfigure(0, weight=1)
        newWindow.grid_rowconfigure(1, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        newWindow.grid_columnconfigure(1, weight=1)
        
        
        samplename = ctk.CTkFrame(newWindow, width=100, corner_radius=2)
        samplename.grid(row=0, column=0, rowspan=1, sticky="nsew")
        samplename.grid_columnconfigure(0, weight=1)
        #self.parameters.grid_columnconfigure(1, weight=1)
        samplename.grid_rowconfigure(0, weight=1)
        
        sequence = ctk.CTkFrame(newWindow, width=100, corner_radius=2)
        sequence.grid(row=1, column=0, rowspan=1, sticky="nsew")
        sequence.grid_columnconfigure(0, weight=1)
        #self.parameters.grid_columnconfigure(1, weight=1)
        sequence.grid_rowconfigure(0, weight=1)
        
        samplename_params = ctk.CTkLabel(samplename, text="SampleNameInfo:", font=ctk.CTkFont(size=20, weight="bold"))
        samplename_params.grid(row=0, column=0, sticky="nsew", padx=20, pady=(20, 20))
        
        param_list_sample = ["ColumnNr", "RowNrBegin", "RowNrEnd"]
        defaults_sample = [3, 2, 10]
        number_sample = 0
        sample_entries = []
        for i, param in enumerate(param_list_sample):
            ctk.CTkLabel(samplename, text=param + ":", width=20).grid(column=0, row=i+1, padx=5, pady=5)
            #text = tk.StringVar()
            param_entry = ctk.CTkEntry(samplename, corner_radius=5) #, textvariable=text
            if (self.sample_params_import[i] == ""):
                param_entry.insert(tk.END, str(defaults_sample[i]))
            else:
                param_entry.insert(tk.END, str(self.sample_params_import[i]))
            param_entry.grid(column=1, row=i+1, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            sample_entries.append(param_entry)
            number_sample += 1
        
        
        
        sequence_params = ctk.CTkLabel(sequence, text="SequenceInfo:", font=ctk.CTkFont(size=20, weight="bold"))
        sequence_params.grid(row=0, column=0, sticky="nsew", padx=20, pady=(20, 20))
        
        param_list_sequence = ["ColumnNr", "RowNrBegin", "RowNrEnd"]
        defaults_sequence = [11, 2, 10]
        number_sequence = 0
        sequence_entries = []
        for i, param in enumerate(param_list_sequence):
            ctk.CTkLabel(sequence, text=param + ":", width=20).grid(column=0, row=i+1, padx=5, pady=5)
            #text = tk.StringVar()
            param_entry = ctk.CTkEntry(sequence, corner_radius=5) #, textvariable=text
            if (self.sequence_params_import[i] == ""):
                param_entry.insert(tk.END, str(defaults_sequence[i]))
            else:
                param_entry.insert(tk.END, str(self.sequence_params_import[i]))
            param_entry.grid(column=1, row=i+1, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            sequence_entries.append(param_entry)
            number_sequence += 1
        
        write_samplesheet_button = tk.Button(newWindow, text="Write Samplesheet", command = lambda x = sample_entries, y = sequence_entries : self.writesamplesheet(x, y))
        write_samplesheet_button.grid(column=0, row=2, sticky="sw", padx=15, pady=4)
    
    def writesamplesheet(self, x, y):
        print("Writing Samplesheet:")
        
        sample_names = []
        sequences = []
        
        sample_columnr = int(x[0].get())
        sample_rownr_begin = int(x[1].get())
        sample_rownr_end = int(x[2].get())
        
        sequence_columnr = int(y[0].get())
        sequence_rownr_begin = int(y[1].get())
        sequence_rownr_end = int(y[2].get())
        
        rows_new = []
        
        with open(self.reference_file, "r") as reference:
            csv_reader = csv.reader(reference, delimiter=';')
            for i, line in enumerate(csv_reader, 1):
                if (i < sample_rownr_begin):
                    continue
                if (i > sample_rownr_end):
                    break
                sample_names.append(line[sample_columnr-1])
                sequences.append(line[sequence_columnr-1])
                rows_new.append([line[sequence_columnr-1], line[sample_columnr-1]])
        
        print(sample_names)
        print(sequences)
        
        with open("samplesheettext.csv", 'w', newline='') as new_samplesheet:
            writer = csv.writer(new_samplesheet, delimiter =';', quoting=csv.QUOTE_MINIMAL)
            for row in rows_new:
                row = [element.replace("\n", "") for element in row]
                print(row)
                writer.writerow(row)
        
        
    def set_sample_import_params2(self):
        print("Setting Sample Import:")
        newWindow = tk.Toplevel(self)
        newWindow.title("Reference Import Parameters")
        newWindow.grid_rowconfigure(0, weight=1)
        newWindow.grid_rowconfigure(1, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        newWindow.grid_columnconfigure(1, weight=1)
        
        entries = []
        
        param_list = ["ColumnNrSample", "IndexNrSample", "StartRow", "EndRow"]
        #default_values = ['Output3', 'hhprimers.txt', 'rawdata', 'samplsheet_hedgehog.csv', 'diploid', 'linux', '',  'Sebastian']
        number = 0
        for i, param in enumerate(param_list):
            ctk.CTkLabel(newWindow, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            #text = tk.StringVar()
            param_entry = ctk.CTkEntry(newWindow, corner_radius=5) #, textvariable=text
            if (len(self.var_list) > 0):
                param_entry.insert(tk.END, self.var_list[i])
            else:
                param_entry.insert(tk.END, self.param_list[i])
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            entries.append(param_entry)
            number += 1
        
        
    
    def get_analysis(self):
        print("Get the analysis file")
        
        #filepath = os.path.normpath("")
    
    def download_data(self, x):
       print("Downloading")
       filename = "Subset"
       
       os.chdir(self.workspace_location)
       
       with open(filename + ".fasta", 'w') as fasta:
            for row in x:
                fasta.write(">" + row[4] + row[5] + "\n")
                fasta.write(str(row[8]) + "\n")
                
       #if self.os_param == "windows":
           #os.startfile(filename + ".fasta")
       #elif self.os_param == "linux":
           #subprocess.run(["xdg-open", filename + ".fasta"], check=True)
    
    def get_reference_dict(self, sample_loci_matrix):
        
        dict_reference = {}
        os.chdir(self.workspace_location)
        
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

    
    def get_genalex_params(self):
        print("Getting Genalex Params:")
        
        newWindow = tk.Toplevel(self)
        newWindow.title("Genalex Parameter")
        newWindow.grid_rowconfigure(0, weight=1)
        newWindow.grid_rowconfigure(1, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        newWindow.grid_columnconfigure(1, weight=1)
        
        genalex_param = ctk.CTkFrame(newWindow, width=50, corner_radius=2)
        genalex_param.grid(row=0, column=0, rowspan=4, sticky="nsew")
        genalex_param.grid_columnconfigure(0, weight=1)
        #genalex_param.grid_rowconfigure(0, weight = 1)
        #genalex_param.grid_rowconfigure(1, weight = 1)
        #genalex_param.grid_rowconfigure(2, weight = 1)
        
        
        
        
        genalex_param_logo = ctk.CTkLabel(genalex_param, text="Genalex Parameter", font=ctk.CTkFont(size=15, weight="bold"))
        genalex_param_logo.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        rb1 = ttk.Radiobutton(genalex_param, text="Organism", variable=self.genalex_parameter, value="Organism")
        rb2 = ttk.Radiobutton(genalex_param, text="Country", variable=self.genalex_parameter, value="Country")
        rb3 = ttk.Radiobutton(genalex_param, text="Locality", variable=self.genalex_parameter, value="Locality")
        
        rb1.grid(row=1, column=0, rowspan=1, sticky="nsew")
        rb2.grid(row=2, column=0, rowspan=1, sticky="nsew")
        rb3.grid(row=3, column=0, rowspan=1, sticky="nsew")
        
        print(self.genalex_parameter.get())
    
    
    
    def import_genalex_latest(self):
        directory_path = filedialog.askdirectory(initialdir=self.workspace_location, title="Select folder")
        print(directory_path)
        directory_path = os.path.normpath(directory_path)
        print(directory_path)
        print()
        #file_path = os.path.join(directory_path, "")
        #print(file_path)
        #print(file_path.split('\\')[-1])
        if (self.os_param == "linux"):
            output_folder = directory_path.split('/')[-1]
            output_path = directory_path
        elif (self.os_param == "windows"):
            output_folder = directory_path.split('\\')[-1]
            output_path = directory_path
        
        loci_dict, sample_loci_matrix = self.get_allele_sample_data(output_path, output_folder)
        
        
        self.import_genalex_simple(loci_dict, sample_loci_matrix, output_folder)
        
    def get_status(self):
        db_instance = SQdatabase(self.local_database_path)
        current_size = db_instance.get_size_of_table()
        nr_unique_sequences = db_instance.get_count_unique_sequence()
        
        self.textbox2.delete(0.0, 'end')
        self.textbox2.insert('end-1c', "Database Status: \n")
        print()
        self.textbox2.insert('end-1c', "Nr. of Records: " + str(current_size) + "\n")
        print()
        self.textbox2.insert('end-1c', "Nr. of Unique Sequences: " + str(nr_unique_sequences))
        
    
        
    def delete_data_for_folder(self):
        
        self.is_working2 = False
        
        print("Choose the folder.")
        directory_path = filedialog.askdirectory(initialdir=self.workspace_location, title="Select folder")
        directory_path = os.path.normpath(directory_path)
        #print(directory_path)
        output = ""
        #print(directory_path.split('\\')[-1])
        
        if (self.os_param == "linux"):
            output = directory_path.split('/')[-1]
        elif (self.os_param == "windows"):
            output = directory_path.split('\\')[-1]
        else:
            print("This is no proper operating system. Either Linux or Windows")
        
        print("Deleting:")
        print(directory_path)
        print(output)
        db_instance = SQdatabase(self.local_database_path)
        db_instance.deleting_data_from_folder(output)
        db_instance.closing()
        
        self.get_status()
    
    def get_subset_data_genalex(self):
        query = "SELECT Project, Organism, Country, Locality, Sample, Loci, AlleleIdx, Length, AlleleSequence FROM ssr_table"

        conditions = []
        parameters = []
        
        organism_include = self.checkbox_states_dict2["Organism"]["Include"]
        organism_exclude = self.checkbox_states_dict2["Organism"]["Exclude"]
        
        loci_include = self.checkbox_states_dict2["Loci"]["Include"]
        loci_exclude = self.checkbox_states_dict2["Loci"]["Exclude"]
        
        sample_include = self.checkbox_states_dict2["Sample"]["Include"]
        sample_exclude = self.checkbox_states_dict2["Sample"]["Exclude"]
        
        # Organism conditions
        if organism_include:
            placeholders = ', '.join(['?'] * len(organism_include))
            conditions.append(f"Organism IN ({placeholders})")
            parameters.extend(organism_include)
        if organism_exclude:
            placeholders = ', '.join(['?'] * len(organism_exclude))
            conditions.append(f"Organism NOT IN ({placeholders})")
            parameters.extend(organism_exclude)
        
        # Loci conditions
        if loci_include:
            placeholders = ', '.join(['?'] * len(loci_include))
            conditions.append(f"Loci IN ({placeholders})")
            parameters.extend(loci_include)
        if loci_exclude:
            placeholders = ', '.join(['?'] * len(loci_exclude))
            conditions.append(f"Loci NOT IN ({placeholders})")
            parameters.extend(loci_exclude)
        
        # Sample conditions
        if sample_include:
            placeholders = ', '.join(['?'] * len(sample_include))
            conditions.append(f"Sample IN ({placeholders})")
            parameters.extend(sample_include)
        if sample_exclude:
            placeholders = ', '.join(['?'] * len(sample_exclude))
            conditions.append(f"Sample NOT IN ({placeholders})")
            parameters.extend(sample_exclude)
        
        #conditions.append("AlleleSequence != ?")
        #parameters.append("")
        
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
        
        #print()
        #print(query)
         
        db_instance = SQdatabase(self.local_database_path)
        rows = db_instance.get_specified_subset(query, parameters)
        
        db_instance.closing()
        
        return rows
        
    
    def show_data(self):
        print(self.checkbox_states_dict2)
        
        #query = "SELECT * FROM ssr_table"
        query = "SELECT Project, Organism, Country, Locality, Sample, Loci, AlleleIdx, Length, AlleleSequence FROM ssr_table"

        conditions = []
        parameters = []
        
        organism_include = self.checkbox_states_dict2["Organism"]["Include"]
        organism_exclude = self.checkbox_states_dict2["Organism"]["Exclude"]
        
        loci_include = self.checkbox_states_dict2["Loci"]["Include"]
        loci_exclude = self.checkbox_states_dict2["Loci"]["Exclude"]
        
        sample_include = self.checkbox_states_dict2["Sample"]["Include"]
        sample_exclude = self.checkbox_states_dict2["Sample"]["Exclude"]
        
        # Organism conditions
        if organism_include:
            placeholders = ', '.join(['?'] * len(organism_include))
            conditions.append(f"Organism IN ({placeholders})")
            parameters.extend(organism_include)
        if organism_exclude:
            placeholders = ', '.join(['?'] * len(organism_exclude))
            conditions.append(f"Organism NOT IN ({placeholders})")
            parameters.extend(organism_exclude)
        
        # Loci conditions
        if loci_include:
            placeholders = ', '.join(['?'] * len(loci_include))
            conditions.append(f"Loci IN ({placeholders})")
            parameters.extend(loci_include)
        if loci_exclude:
            placeholders = ', '.join(['?'] * len(loci_exclude))
            conditions.append(f"Loci NOT IN ({placeholders})")
            parameters.extend(loci_exclude)
        
        # Sample conditions
        if sample_include:
            placeholders = ', '.join(['?'] * len(sample_include))
            conditions.append(f"Sample IN ({placeholders})")
            parameters.extend(sample_include)
        if sample_exclude:
            placeholders = ', '.join(['?'] * len(sample_exclude))
            conditions.append(f"Sample NOT IN ({placeholders})")
            parameters.extend(sample_exclude)
        
        #conditions.append("AlleleSequence != ?")
        #parameters.append("")
        
        if conditions:
            query += " WHERE " + " AND ".join(conditions)
        
        print()
        #print(query)
        
        
        db_instance = SQdatabase(self.local_database_path)
        rows = db_instance.get_specified_subset(query, parameters)
        
        print()
        #print(rows)
        
        self.display_data(rows)
        
        db_instance.closing()
    
    def display_data(self, rows):
        print("Displaying")
        
        newWindow = tk.Toplevel(self)
    
        newWindow.grid_rowconfigure(0, weight=1)
        newWindow.grid_rowconfigure(1, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        newWindow.grid_columnconfigure(1, weight=1)
        
        columns = ('Project', 'Organism', 'Country', 'Locality', 'SampleName', 'LociName', 'AlleleIdx', 'AlleleLength', 'AlleleSequence')
        global tree_show
        tree_show = ttk.Treeview(newWindow, columns=columns, show='headings')
        tree_show.heading('Project', text='Project', anchor=tk.CENTER)
        tree_show.heading('Organism', text='Organism', anchor=tk.CENTER)
        tree_show.heading('Country', text='Country', anchor=tk.CENTER)
        tree_show.heading('Locality', text='Locality', anchor=tk.CENTER)
        tree_show.heading('SampleName', text='Sample', anchor=tk.CENTER)
        tree_show.heading('LociName', text='Loci', anchor=tk.CENTER)
        tree_show.heading('AlleleIdx', text='Index', anchor=tk.CENTER)
        tree_show.heading('AlleleLength', text='Length', anchor=tk.CENTER)
        tree_show.heading('AlleleSequence', text='Sequence', anchor=tk.CENTER)
        
        tree_show.grid_rowconfigure(0, weight=1)
        tree_show.grid_rowconfigure(1, weight=1)
        tree_show.grid_columnconfigure(0, weight=1)
        tree_show.grid_columnconfigure(1, weight=1)
        #tree.grid_columnconfigure(2, weight=3)
        
        scrollbar_y = ttk.Scrollbar(newWindow, orient=tk.VERTICAL)
        scrollbar_y.configure(command = tree_show.yview)
        scrollbar_y.grid(row=0, column=1, sticky=tk.NS)
        
        scrollbar_x = ttk.Scrollbar(newWindow, orient=tk.HORIZONTAL)
        scrollbar_x.configure(command = tree_show.xview)
        scrollbar_x.grid(row=1, column=0, sticky=tk.EW)
        
    
        for row in rows:
            #row_element = ""
            row_element = []
            for element in row:
                #row_element += str(element) + "  \t\t\t\t\t\t  "
                row_element.append(str(element)) # + "  \t\t\t\t\t\t  ")
            tree_show.insert('', tk.END, values=row_element)
            #tree.column("#0", stretch=False)
        
        
        tree_show.bind('<<TreeviewSelect>>', lambda event, x = tree_show : self.selecting(x))
        #tree_show.bind('<<TreeviewEnter>>', self.download)
        
        tree_show.grid(row = 0, column = 0, sticky=tk.NSEW, padx = 1, pady = 1)
        
        tree_show.configure(yscroll=scrollbar_y.set)
        tree_show.configure(xscroll=scrollbar_x.set)
        #scrollbar_y.grid(row=0, column=0, sticky=tk.NS)
        tree_show.configure(selectmode = "extended")
        print("This is it.")
        
        self.button_download = ttk.Button(newWindow, text = "Download data", command = lambda x = rows : self.download_data(x))
        self.button_download.grid(column = 0, row = 2, sticky = tk.W, padx=20, pady=20) #
        
        self.button_download_special = ttk.Button(newWindow, text = "Import for genalex", command = lambda x = rows : self.import_genalex(x))
        self.button_download_special.grid(column = 1, row = 2, sticky = tk.W, padx=20, pady=20) #
        
        
    def import_genalex_simple(self, loci_dict, sample_loci_matrix, outputfolder):
        print("Getting the genalex output for the outputfolder " + outputfolder + " \n")
        
        sample_number = len(sample_loci_matrix)
        loci_list_unique = sorted(list(set([locus for locus in loci_dict.keys()])))
        if (self.ploidy == "diploid"):
            loci_length = len(loci_list_unique)*2
            loci_list = [loci for loci in loci_list_unique for _ in range(2)]
        else:
            loci_length = len(loci_list_unique)
            loci_list = loci_list_unique
        
        print()
        #print(loci_list)
        print()
        #print(outputfolder)
        print()
        
        os.chdir(self.workspace_location)
        
        first_row = []
        second_row = ["", "", ""]
        
        first_row.append(int(len(loci_list_unique)))
        first_row.append(int(sample_number))
        
        
        reference_dict, list_info = self.get_reference_dict(sample_loci_matrix)
        
        print(reference_dict)
        print()
        
        #In this step get the parameter from the scrollframe
        #self.genalex_parameter = "Organism"
        
        param_length = 0
        param_name = ""
        
        full_list = []
        unique_list = []
        
        print(self.genalex_parameter.get())
        if (self.genalex_parameter.get() == "Organism"):
            full_list = list_info[1]
            param_name = "Organism"
        elif (self.genalex_parameter.get() == "Country"):
            full_list = list_info[2]
            param_name = "Country"
        elif (self.genalex_parameter.get() == "Locality"):
            full_list = list_info[3]
            param_name = "Locality"
        
        unique_list = list(set(full_list))
        param_length = len(unique_list)
        first_row.append(param_length)
        for element in unique_list:
            second_row.append(element)
            count_element = full_list.count(element)
            first_row.append(count_element)
            

        with open('genalex_' + outputfolder + '.csv', 'w', newline='') as file_csv:
            writer = csv.writer(file_csv, delimiter =',', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(first_row)
            writer.writerow(second_row)
            writer.writerow(["samples"] + [param_name] + loci_list)
            print(["samples"] + loci_list)
            #next_row = [0] * row_length
            
            for param in unique_list:
                for sample, loci in sample_loci_matrix.items(): # and (reference_dict[sample][param_name] == param)
                    if (sample in reference_dict and (reference_dict[sample][param_name] == param)): #not(reference_dict[sample] is None
                        #print("NOW:")
                        #print(reference_dict[sample]["Organism"])
                        next_row = [0] * loci_length
                        if (self.ploidy == "diploid"):
                            for i in range(0, len(loci_list), 2):
                                locus = loci_list[i]
                            #for i, locus in enumerate(loci_list_unique):
                                if (locus in sample_loci_matrix[sample] and len(list(sample_loci_matrix[sample][locus].keys())) > 1):
                                    print("Reached")
                                    next_row[i] = list(sample_loci_matrix[sample][locus].keys())[0]
                                    next_row[i+1] = list(sample_loci_matrix[sample][locus].keys())[1]
                                else:
                                    next_row[i] = list(sample_loci_matrix[sample][locus].keys())[0]
                                    next_row[i+1] = list(sample_loci_matrix[sample][locus].keys())[0]
                                print([sample] + next_row)
                        
                        elif (self.ploidy == "haploid"):
                            for i in range(0, len(loci_list), 1):
                                locus = loci_list[i]
                                if (locus in sample_loci_matrix[sample]):
                                    print("Reached")
                                    next_row[i] = list(sample_loci_matrix[sample][locus].keys())[0]
                        
                        
                        writer.writerow([sample] + [param] + next_row)
                            #else:
                                #writer.writerow([sample] + ["None"] + next_row)
            
            for sample, loci in sample_loci_matrix.items(): # and (reference_dict[sample][param_name] == param)
                if (not(sample in reference_dict)):
                    #print("NOW:")
                    #print(reference_dict[sample]["Organism"])
                    next_row = [0] * loci_length
                    if (self.ploidy == "diploid"):
                        for i in range(0, len(loci_list), 2):
                            locus = loci_list[i]
                        #for i, locus in enumerate(loci_list_unique):
                            if (locus in sample_loci_matrix[sample] and len(list(sample_loci_matrix[sample][locus].keys())) > 1):
                                print("Reached")
                                next_row[i] = list(sample_loci_matrix[sample][locus].keys())[0]
                                next_row[i+1] = list(sample_loci_matrix[sample][locus].keys())[1]
                            else:
                                next_row[i] = list(sample_loci_matrix[sample][locus].keys())[0]
                                next_row[i+1] = list(sample_loci_matrix[sample][locus].keys())[0]
                            print([sample] + next_row)
                    
                    elif (self.ploidy == "haploid"):
                        for i in range(0, len(loci_list), 1):
                            locus = loci_list[i]
                            if (locus in sample_loci_matrix[sample]):
                                print("Reached")
                                next_row[i] = list(sample_loci_matrix[sample][locus].keys())[0]
                    
                    
                    writer.writerow([sample] + ["None"] + next_row)

        if self.os_param == "windows":
            os.startfile('genalex_' + outputfolder + '.csv')
        elif self.os_param == "linux":
            subprocess.run(["xdg-open", 'genalex_' + outputfolder + '.csv'], check=True)
        
     
    
    def get_subset_dict(self, rows):
        print("Getting the subset as a dict")
        
        sub_dict = {}
        
        #First get all loci and samples from the rows and make them unique:
        loci_list = list(set([row[5] for row in rows]))
        sample_list = list(set([row[4] for row in rows]))
        
        #Building the dict: {sample : {locus : {AlleleNumber : Sequence}}}
        for sample in sample_list:
            sub_dict[sample] = {}
            for locus in loci_list:
                if locus not in sub_dict[sample].keys():
                    sub_dict[sample][locus] = {}
                for row in rows:
                    if (row[5] == locus and row[4] == sample):
                        sub_dict[sample][locus][row[6]] = row[8]
        
        return sub_dict, loci_list, sample_list
    
    
    
    
    def import_genalex(self, x):
        
        #rows in the order: Project, Organism, Country, Locality, Sample, Loci, AlleleIdx, Length, AlleleSequence
        rows = self.get_subset_data_genalex()
        
        sub_dict, loci_list_unique, sample_list = self.get_subset_dict(rows)
        
        if (self.ploidy == "diploid"):
            loci_list = [loci for loci in loci_list_unique for _ in range(2)]
        else:
            loci_list = loci_list_unique
        
        loci_length_unique = len(loci_list_unique)
        loci_length = len(loci_list)
        sample_length = len(sub_dict)
        
        os.chdir(self.workspace_location)
        
        first_row = []
        second_row = ["", "", ""]
        
        first_row.append(int(len(loci_list_unique)))
        first_row.append(int(sample_length))
        
        reference_dict, list_info = self.get_reference_dict(sub_dict)
        
        
        param_length = 0
        param_name = ""
        
        full_list = []
        unique_list = []
        
        if (self.genalex_parameter.get() == "Organism"):
            full_list = list_info[1]
            param_name = "Organism"
        elif (self.genalex_parameter.get() == "Country"):
            full_list = list_info[2]
            param_name = "Country"
        elif (self.genalex_parameter.get() == "Locality"):
            full_list = list_info[3]
            param_name = "Locality"
        
        unique_list = list(set(full_list))
        param_length = len(unique_list)
        first_row.append(param_length)
        for element in unique_list:
            second_row.append(element)
            count_element = full_list.count(element)
            first_row.append(count_element)
            
        
        
        with open('genalex_subset.csv', 'w', newline='') as file_csv:
            writer = csv.writer(file_csv, delimiter =',', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(first_row)
            writer.writerow(second_row)
            writer.writerow(["samples"] + [param_name] + loci_list)
            print(["samples"] + loci_list)
            #next_row = [0] * row_length
            
            for param in unique_list:
                for sample, loci in sub_dict.items(): # and (reference_dict[sample][param_name] == param)
                    if (sample in reference_dict and (reference_dict[sample][param_name] == param)): #not(reference_dict[sample] is None
                        #print("NOW:")
                        #print(reference_dict[sample]["Organism"])
                        next_row = [0] * loci_length
                        if (self.ploidy == "diploid"):
                            for i in range(0, len(loci_list), 2):
                                locus = loci_list[i]
                            #for i, locus in enumerate(loci_list_unique):
                                if (locus in sub_dict[sample] and len(list(sub_dict[sample][locus].keys())) > 1):
                                    next_row[i] = list(sub_dict[sample][locus].keys())[0]
                                    next_row[i+1] = list(sub_dict[sample][locus].keys())[1]
                                elif (locus in sub_dict[sample] and len(list(sub_dict[sample][locus].keys())) == 1):
                                    next_row[i] = list(sub_dict[sample][locus].keys())[0]
                                    next_row[i+1] = list(sub_dict[sample][locus].keys())[0]
                                print([sample] + next_row)
                        
                        elif (self.ploidy == "haploid"):
                            for i in range(0, len(loci_list), 1):
                                locus = loci_list[i]
                                if (locus in sub_dict[sample]):
                                    print("Reached")
                                    next_row[i] = list(sub_dict[sample][locus].keys())[0]
                        
                        
                        writer.writerow([sample] + [param] + next_row)
            
            for sample, loci in sub_dict.items(): # and (reference_dict[sample][param_name] == param)
                if (not(sample in reference_dict)):
                    #print("NOW:")
                    #print(reference_dict[sample]["Organism"])
                    next_row = [0] * loci_length
                    if (self.ploidy == "diploid"):
                        for i in range(0, len(loci_list), 2):
                            locus = loci_list[i]
                        #for i, locus in enumerate(loci_list_unique):
                            if (locus in sub_dict[sample] and len(list(sub_dict[sample][locus].keys())) > 1):
                                print("Reached")
                                next_row[i] = list(sub_dict[sample][locus].keys())[0]
                                next_row[i+1] = list(sub_dict[sample][locus].keys())[1]
                            elif (locus in sub_dict[sample] and len(list(sub_dict[sample][locus].keys())) == 1):
                                next_row[i] = list(sub_dict[sample][locus].keys())[0]
                                next_row[i+1] = list(sub_dict[sample][locus].keys())[0]
                                
                            print([sample] + next_row)
                    
                    elif (self.ploidy == "haploid"):
                        for i in range(0, len(loci_list), 1):
                            locus = loci_list[i]
                            if (locus in sub_dict[sample]):
                                print("Reached")
                                next_row[i] = list(sub_dict[sample][locus].keys())[0]
                    
                    
                    writer.writerow([sample] + ["None"] + next_row)


        if self.os_param == "windows":
            os.startfile('genalex_subset.csv')
        elif self.os_param == "linux":
            subprocess.run(["xdg-open", 'genalex_subset.csv'], check=True)
        
        
        
    
    
    def import_genalex2(self, x):
        rows = self.get_data_genalex()
        sub_dict = {}
    
        print()
        print()
        print("For Genalex:")
        #print(rows)
        print()
        list_loci = []
        #list_loci.append("samples")
        for row in rows:
            list_loci.append(row[5])
        
        list_loci_unique = set(list_loci) #each loci will be unique
        print(list_loci_unique)
        if (self.ploidy == "diploid"):
            loci_list = [loci for loci in list_loci_unique for _ in range(2)]
        else:
            loci_list = list_loci_unique
            
        
        print()
        print("All loci:")
        loci_length = len(loci_list)
            
        #loci_list = ["samples"] + loci_list
        print(loci_list)
        
        sample_list = list(set([row[4] for row in rows]))
        print()
        print("All Samples:")
        print(sample_list)
        
        for sample in sample_list:
            sub_dict[sample] = {}
            for locus in loci_list:
                if locus not in sub_dict[sample].keys():
                    sub_dict[sample][locus] = {}
                for row in rows:
                    if (row[5] == locus and row[4] == sample):
                        sub_dict[sample][locus][row[6]] = row[8]
                        
        
        #Project, Organism, Country, Locality, Sample, Loci, AlleleIdx, Length, AlleleSequence
        print()
        print("Full dict:")
        print(sub_dict)
        print()
        
        
        with open('genalex_subset.csv', 'w', newline='') as file_csv:
            writer = csv.writer(file_csv, delimiter =',', quoting=csv.QUOTE_MINIMAL)
            writer.writerow("")
            writer.writerow("")
            writer.writerow(["samples"] + loci_list)
            print(["samples"] + loci_list)
            #next_row = [0] * row_length
            for sample, loci in sub_dict.items():
                next_row = [0] * loci_length
                if (self.ploidy == "diploid"):
                    for i in range(0, len(loci_list), 2):
                        locus = loci_list[i]
                    #for i, locus in enumerate(loci_list_unique):
                        if (locus in sub_dict[sample] and len(list(sub_dict[sample][locus].keys())) > 1):
                            print("Reached")
                            next_row[i] = list(sub_dict[sample][locus].keys())[0]
                            next_row[i+1] = list(sub_dict[sample][locus].keys())[1]
                        elif (locus in sub_dict[sample] and len(list(sub_dict[sample][locus].keys())) == 1):
                            next_row[i] = list(sub_dict[sample][locus].keys())[0]
                            next_row[i+1] = list(sub_dict[sample][locus].keys())[0]
                        print([sample] + next_row)
                
                elif (self.ploidy == "haploid"):
                    for i in range(0, len(loci_list), 1):
                        locus = loci_list[i]
                        if (locus in sub_dict[sample] and len(list(sub_dict[sample][locus].keys())) == 1):
                            print("Reached")
                            next_row[i] = list(sub_dict[sample][locus].keys())[0]
                            
                writer.writerow([sample] + next_row)
        
    
        
    def selecting(self, tree_object):
        print("Test")
        #print(tree_object.selection())
        for selected_item in tree_object.selection():
            #print(selected_item) #gives me the indizes of the item
            item = tree_object.item(selected_item)
            record = item['values']
            print(record)
            all_info = ""
            for element in record:
                all_info += str(element) + "\n\n"
                
            tk.messagebox.showinfo(title='Information', message=all_info)
    
    
    def add_data_to_local_database(self):
        
        self.textbox2.delete(0.0, 'end')
        self.textbox2.insert('end-1c', "Adding data to local Database! \n")
        
        thread = threading.Thread(target=self.adding_dataset)
        thread.start()
        
        self.is_working2 = True
        # Start appending periods to the textbox
        self.append_period2()
        
        #self.get_status()
    
    
    def adding_dataset(self):
        print("Adding data")
        directory_path = filedialog.askdirectory(initialdir=self.workspace_location, title="Select folder")
        #print(directory_path)
        
        #self.is_working2 = True
        directory_path = os.path.normpath(directory_path)
        #print(directory_path)
        print()
        #file_path = os.path.join(directory_path, "")
        #print(file_path)
        #print(file_path.split('\\')[-1])
        if (self.os_param == "linux"):
            output_folder = directory_path.split('/')[-1]
            output_path = directory_path
        elif (self.os_param == "windows"):
            output_folder = directory_path.split('\\')[-1]
            output_path = directory_path
        print()
        #print(output_folder)
        
        
        
        loci_dict, sample_loci_matrix = self.get_allele_sample_data(output_path, output_folder)
        
        #print(loci_dict)
        print()
        #print(sample_loci_matrix)
        print()
        
        self.add_to_local_database(output_folder, sample_loci_matrix, loci_dict)
        
        self.is_working2 = False
        
        self.get_status()
    
        
    def add_to_local_database(self, output_folder, sample_loci_matrix, loci_dict):
        print("Truly adding the data")    
        
        dict_reference = {}
        os.chdir(self.workspace_location)
        
        dict_reference, full_list = self.get_reference_dict(sample_loci_matrix)
        
        
        print()
        print("The dict:")
        print(dict_reference)
        
        db_instance = SQdatabase(self.local_database_path)
        
        for sample_name in list(sample_loci_matrix.keys()):
            for loci_name in list(sample_loci_matrix[sample_name].keys()):
                for number in list(sample_loci_matrix[sample_name][loci_name].keys()): #"WF-Project SSRGBAS" "Plant-oak"
                    if (not(len(dict_reference) == 0) and not(dict_reference.get(sample_name) == None)):
                        #print(new_dict[sample_name].keys()[0])
                        db_instance.insert_record([output_folder, self.projectfolder, dict_reference[sample_name]["Project"], dict_reference[sample_name]["Organism"], dict_reference[sample_name]["Country"], dict_reference[sample_name]["Locality"], sample_name, loci_name, number, len(sample_loci_matrix[sample_name][loci_name][number]), sample_loci_matrix[sample_name][loci_name][number], "Default"])
                    else:
                        db_instance.insert_record([output_folder, self.projectfolder, "None", "None", "None", "None", sample_name, loci_name, number, len(sample_loci_matrix[sample_name][loci_name][number]), sample_loci_matrix[sample_name][loci_name][number], "Default"])
        
        sequence_rows = db_instance.get_notempty_sequence_rows()
        print()
        print(sequence_rows)
        print()
        number = db_instance.get_size_of_table()
        print(number)
        number2 = db_instance.get_notempty_sequence_rows_number()
        print(number2)
        
        db_instance.closing()
        
        print()
        print("Added to dabatase")
        print()

    
    def add_to_local_database2(self, output_folder, sample_loci_matrix, loci_dict):
        print("Truly adding the data")    
        
        dict_reference = {}
        os.chdir(self.workspace_location)
        
        if os.path.exists(self.reference_file_path):
            with open(self.reference_file_path, "r") as reference:
                csv_reader = csv.reader(reference, delimiter=';')
                for i, line in enumerate(csv_reader):
                    print(line)
                    #line = re.split(',|\t|;', line)
                    #sample_ref = line_splitted[2]
                    sample_ref = line[2]
                    # Create a CSV reader object
                    
                    for sample in sample_loci_matrix.keys():
                        if (sample_ref == sample):
                            dict_reference[sample_ref] = {}
                            dict_reference[sample_ref][line[0]] = {}
                            dict_reference[sample_ref][line[0]][line[1]] = {}
                            dict_reference[sample_ref][line[0]][line[1]][line[5]] = line[6]
                    #self.sample_list.append(re.split(',|\t|;', line)[0])
                #reader = csv.reader(reference, delimiter=";")
            #print(dict_reference)
            
        self.reference_dict = dict_reference
        
        print()
        print("The dict:")
        print(dict_reference)
        
        db_instance = SQdatabase(self.local_database_path)
        
        for sample_name in list(sample_loci_matrix.keys()):
            for loci_name in list(sample_loci_matrix[sample_name].keys()):
                for number in list(sample_loci_matrix[sample_name][loci_name].keys()): #"WF-Project SSRGBAS" "Plant-oak"
                    if (not(len(dict_reference) == 0) and not(dict_reference.get(sample_name) == None)):
                        #print(new_dict[sample_name].keys()[0])
                        db_instance.insert_record([output_folder, self.projectfolder, list(dict_reference[sample_name].keys())[0], list(dict_reference[sample_name][list(dict_reference[sample_name].keys())[0]].keys())[0], list(dict_reference[sample_name][list(dict_reference[sample_name].keys())[0]][list(dict_reference[sample_name][list(dict_reference[sample_name].keys())[0]].keys())[0]].keys())[0], "", sample_name, loci_name, number, len(sample_loci_matrix[sample_name][loci_name][number]), sample_loci_matrix[sample_name][loci_name][number], "Default"])
                    else:
                        db_instance.insert_record([output_folder, self.projectfolder, "None", "None", "None", "", sample_name, loci_name, number, len(sample_loci_matrix[sample_name][loci_name][number]), sample_loci_matrix[sample_name][loci_name][number], "Default"])
        
        sequence_rows = db_instance.get_notempty_sequence_rows()
        print()
        print(sequence_rows)
        print()
        number = db_instance.get_size_of_table()
        print(number)
        number2 = db_instance.get_notempty_sequence_rows_number()
        print(number2)
        
        db_instance.closing()
        
        print()
        print("Added to dabatase")
        print()
        
    def get_allele_sample_data(self, output_path, outputfolder):
        print("Get the data.")
        loci_dict = {}
        sample_loci_matrix = {}
        
        current_working_path = os.getcwd()
        for subdir, dirs, files in os.walk(output_path):
            #print(subdir)
            splitted_subdir = ""
            if (self.os_param == "linux"):
                splitted_subdir = subdir.split('/')[-1]
            elif (self.os_param == "windows"):
                splitted_subdir = subdir.split('\\')[-1]
            else:
                print("Neither Linux nor Windows?")
                raise Exception("AlleleFolder cannot be found!")
            print(splitted_subdir)
            if (splitted_subdir == "AlleleCall"):
                print("Now")
                print(subdir)
                allele_call_folder_path = os.path.normpath(subdir)
                print(allele_call_folder_path)
                os.chdir(allele_call_folder_path)
                for file in os.listdir(allele_call_folder_path):
                    if ("allelle_list" in file):
                        with open(file, "r") as allele_list:
                            print("Found")
                            for line in allele_list:
                                #print("Reached")
                                led = line.rstrip('\r\n')
                                #print(led)
                                if len(led.split('_')) == 2: # This identify the marker name. All markers names have an underscore, thus if you slit it it will return a list of two ellements
                                    locus = '_'.join(led.split('_')[:2])
                                    print(locus)
                                    loci_dict[locus] = {}
                                elif (led == '\\' or led == '' or led == "/"):
                                    print("End")
                                else:
                                    led = led.split('\t')
                                    allele_number = led[0].rstrip(':')
                                    sequence = led[1].replace('-', '')
                                    if (locus in loci_dict.keys()):
                                        loci_dict[locus][int(allele_number)] = sequence
                    elif ("matrix" in file):
                        with open(file, "r") as matrix:
                            i = 0
                            for line in matrix:
                                print(line)
                                led = line.rstrip('\r\n')
                                if (i == 0):
                                    loci_names = led.split()
                                elif (i > 0):
                                    splitted_line = line.split()
                                    sample_name = splitted_line[0]
                                    sample_loci_matrix[sample_name] = {}
                                    for i in range(1, len(loci_names), 1):
                                        sample_loci_matrix[sample_name][loci_names[i]] = {}
                                    for j in range(1, len(splitted_line), 1):
                                        if (not(str(splitted_line[j]) == "")): #without len()-1
                                            if (not(int(splitted_line[j]) == 0)):
                                                sample_loci_matrix[sample_name][loci_names[j]][int(splitted_line[j])] = loci_dict[loci_names[j]][int(splitted_line[j])]
                                            else:
                                                sample_loci_matrix[sample_name][loci_names[j]][int(splitted_line[j])] = ""
                                        
                                i += 1
                            
                            print()
                            print(loci_names)
                            print()
                            print("Finished")
                    
                                    
                break
            os.chdir(self.workspace_location)
            self.sample_loci_matrix = sample_loci_matrix
            self.loci_dict = loci_dict
        return loci_dict, sample_loci_matrix
                        
    #Import Genalex        
    def import_data(self, directory, output):
        current_working_path = os.getcwd()
        for subdir, dirs, files in os.walk(directory):
            print(subdir.split('\\')[-1]) #subdir contains all possible folders in my directory, including directory itself
            #We split up the folder paths with '\\' so we only get the last name of the path, namely the folder name itself
            #The split function returns to us a list of elements splitted through the '\\' character -> \
            #[-1] indicates the last element of the list, namely the folder name itself
            if (subdir.split('\\')[-1] == output): #finding our correct output directory
                print("")
                direct = subdir + "\AlleleCall" #entering into our AlleleCall Directory of the output directory
                #global loci_dict
                #loci_dict = {}
                for subdir, dirs, files in os.walk(direct):
                    #print(files[0])
                    os.chdir(direct) #changing path to the AlleleCall folder
                    with open(files[0], "r") as file1: #We have 2 txt. files in our AlleleCall folder
                        lines = (line.rstrip() for line in file1)
                        for line in lines:
                            #print(line)
                            if (not(line == "\\") and line.split() and len(line.split()) == 1):  #different condition: not(line.split()[0].isdigit()) and not(line == "\\") and line.split() and len(line.split()) == 1
                                loci = line.split()[0]
                                #print(loci)
                                self.loci_dict[loci] = {}
                            elif (line[0].isdigit() and not(line == "\\") and line.split()):
                                line_splitted = line.split(":", 1)
                                number = line_splitted[0]
                                sequence = line_splitted[1].join(line_splitted[1].split())
                                self.loci_dict[loci][int(number)] = sequence      
                        #for key in loci_dict.copy(): #eliminating loci_names with no alleles from our dictionary
                            #if (len(loci_dict[key]) == 0):
                                #loci_dict.pop(key)
                    
                    with open(files[1], "r") as file2: #opening our second txt. file (matrix)
                        #global sample_loci_matrix
                        #sample_loci_matrix = {}
                        k = 0
                        for line in file2:
                            if (k == 0):
                                loci_names = line.split() #getting our first row of loci names
                            elif (k > 0):
                                splitted_line = line.split()
                                sample_name = splitted_line[0]
                                self.sample_loci_matrix[sample_name] = {}
                                for i in range(1, len(loci_names)-1, 1):
                                    self.sample_loci_matrix[sample_name][loci_names[i]] = {}
                                for j in range(1, len(splitted_line)-1, 1):
                                    if (not(str(splitted_line[j]) == "")):
                                        if (not(int(splitted_line[j]) == 0)):
                                            self.sample_loci_matrix[sample_name][loci_names[j]][int(splitted_line[j])] = self.loci_dict[loci_names[j]][int(splitted_line[j])]
                                        else:
                                            self.sample_loci_matrix[sample_name][loci_names[j]][int(splitted_line[j])] = ""
                                #for key in sample_loci_matrix[sample_name].copy():
                                    #if (len(sample_loci_matrix[sample_name][key]) == 0):
                                        #sample_loci_matrix[sample_name].pop(key)
                            k += 1
                break  #to avoid iterating over the other files, subdirs and dirs, we break once we reach the the first Path    
        os.chdir(current_working_path)
        print(self.loci_dict)
        print()
        #print(self.sample_loci_matrix)
        return self.loci_dict, self.sample_loci_matrix, current_working_path
    
    
    
    def add_dataset(self):
        print("Adding data")
    
    def updating_checkboxes(self, checkboxes, values, search_var, scrollframe):
        for checkbox in checkboxes:
            # Store the state of the checkbox before removing it
            self.checkbox_states_dict[checkbox.cget("text")] = checkbox.get()
            checkbox.grid_forget()
        
        checkboxes.clear()
    
        query = search_var.get().lower()
        filtered_values = [value for value in values if query in value.lower()]
        
        for i, value in enumerate(filtered_values):
            # Check if there's a stored state for this checkbox
            state = self.checkbox_states_dict.get(value, 0)
                
            checkbox_var = tk.IntVar(value=state)
            checkbox = ctk.CTkCheckBox(scrollframe, text=value, variable=checkbox_var)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(5, 5), sticky="w")
            checkboxes.append(checkbox)
    
    
    
        
    
    def on_closing(self):
        #self.checkbox_states_dict2["Organism"]["Include"] = []
        list_organism_include = []
        for checkbox in self.checkboxes_organism_include:
            checkbox_text = checkbox.cget("text")
            print(checkbox_text)
            if checkbox.get() == 1:
                list_organism_include.append(checkbox_text)
                #self.checkbox_states_dict2["Organism"]["Include"].append(checkbox_text)
            self.checkbox_states_include_dict[checkbox_text] = checkbox.get()
        print("After closing:")
        print(list_organism_include)
        print()
        self.checkbox_states_dict2["Organism"]["Include"] = list(set(list_organism_include))
        self.checkboxes_organism_include = []
        
        
        list_organism_exclude = []
        for checkbox in self.checkboxes_organism_exclude:
            checkbox_text = checkbox.cget("text")
            print(checkbox_text)
            if checkbox.get() == 1:
                list_organism_exclude.append(checkbox_text)
                #self.checkbox_states_dict2["Organism"]["Include"].append(checkbox_text)
            self.checkbox_states_exclude_dict[checkbox_text] = checkbox.get()
        print("After closing:")
        print(list_organism_exclude)
        print()
        self.checkbox_states_dict2["Organism"]["Exclude"] = list(set(list_organism_exclude))
        self.checkboxes_organism_exclude = []
        
        
        list_loci_include = []
        for checkbox in self.checkboxes_loci_include:
            checkbox_text = checkbox.cget("text")
            print(checkbox_text)
            if checkbox.get() == 1:
                list_loci_include.append(checkbox_text)
                #self.checkbox_states_dict2["Organism"]["Include"].append(checkbox_text)
            self.checkbox_states_include_dict[checkbox_text] = checkbox.get()
        print("After closing:")
        print(list_loci_include)
        print()
        self.checkbox_states_dict2["Loci"]["Include"] = list(set(list_loci_include))
        self.checkboxes_loci_include = []
        
        
        list_loci_exclude = []
        for checkbox in self.checkboxes_loci_exclude:
            checkbox_text = checkbox.cget("text")
            print(checkbox_text)
            if checkbox.get() == 1:
                list_loci_exclude.append(checkbox_text)
                #self.checkbox_states_dict2["Organism"]["Include"].append(checkbox_text)
            self.checkbox_states_exclude_dict[checkbox_text] = checkbox.get()
        print("After closing:")
        print(list_loci_exclude)
        print()
        self.checkbox_states_dict2["Loci"]["Exclude"] = list(set(list_loci_exclude))
        self.checkboxes_loci_exclude = []
        
        
        
        list_sample_include = []
        for checkbox in self.checkboxes_sample_include:
            checkbox_text = checkbox.cget("text")
            print(checkbox_text)
            if checkbox.get() == 1:
                list_sample_include.append(checkbox_text)
                #self.checkbox_states_dict2["Organism"]["Include"].append(checkbox_text)
            self.checkbox_states_include_dict[checkbox_text] = checkbox.get()
        print("After closing:")
        print(list_sample_include)
        print()
        self.checkbox_states_dict2["Sample"]["Include"] = list(set(list_sample_include))
        self.checkboxes_sample_include = []
        
        
        list_sample_exclude = []
        for checkbox in self.checkboxes_sample_exclude:
            checkbox_text = checkbox.cget("text")
            print(checkbox_text)
            if checkbox.get() == 1:
                list_sample_exclude.append(checkbox_text)
                #self.checkbox_states_dict2["Organism"]["Include"].append(checkbox_text)
            self.checkbox_states_exclude_dict[checkbox_text] = checkbox.get()
        print("After closing:")
        print(list_sample_exclude)
        print()
        self.checkbox_states_dict2["Sample"]["Exclude"] = list(set(list_sample_exclude))
        self.checkboxes_sample_exclude = []
        
        
        
        print()
        print("Whole dict:")
        print(self.checkbox_states_dict2)
        print()
        newWindowX.destroy()
    
    
    def extract_subset(self):
        print("Extracting")
        global newWindowX
        newWindowX = tk.Toplevel(self)
        newWindowX.title("Extract subset:")
        newWindowX.grid_rowconfigure(0, weight=1)
        #newWindow.grid_rowconfigure(1, weight=1)
        newWindowX.grid_columnconfigure(0, weight=1)
        newWindowX.grid_columnconfigure(1, weight=0)
        newWindowX.grid_columnconfigure(2, weight=1)
        newWindowX.grid_columnconfigure(3, weight=0)
        newWindowX.grid_columnconfigure(4, weight=1)
        #newWindow.grid_columnconfigure(1, weight=1)
        
        db_instance = SQdatabase(self.local_database_path)
        
        
        organism = ctk.CTkFrame(newWindowX, width=50, corner_radius=2)
        organism.grid(row=0, column=0, rowspan=4, sticky="nsew")
        organism.grid_columnconfigure(0, weight=1)
        organism.grid_rowconfigure(0, weight = 0)
        organism.grid_rowconfigure(1, weight = 1)
        
        
        
        
        logo_organism_include = ctk.CTkLabel(organism, text="Include Organism", font=ctk.CTkFont(size=20, weight="bold"))
        logo_organism_include.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        organism_scrollframe_include = ctk.CTkScrollableFrame(organism, width=50, corner_radius=2)
        organism_scrollframe_include.grid(row=1, column=0, rowspan=1, sticky="nsew")
        organism_scrollframe_include.grid_columnconfigure(0, weight=1)
        
        
        search_var_organism_include = tk.StringVar()
        search_entry_organism_include = ctk.CTkEntry(organism_scrollframe_include, textvariable=search_var_organism_include)
        search_entry_organism_include.grid(row=0, column=0, padx=10, pady=(5, 5), sticky="w")
        
        unique_organisms = db_instance.get_unique_organism()
        
        checkboxes_organism_include = []
        #self.values_organism_include = ["organism1", "organism2", "organism3", "organism4"] #call from database
        self.values_organism_include = [organism[0] for organism in unique_organisms]
        for i, value in enumerate(self.values_organism_include):
            state = self.checkbox_states_include_dict.get(value, 0)
    
            checkbox_var = tk.IntVar(value=state)
            checkbox = ctk.CTkCheckBox(organism_scrollframe_include, text=value, variable=checkbox_var)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(5, 5), sticky="w")
            self.checkboxes_organism_include.append(checkbox)
            if (checkbox.get() == 1):
                checkboxes_organism_include.append(checkbox.cget("text"))
            
        self.checkbox_states_dict2["Organism"]["Include"] = list(set(checkboxes_organism_include))
        
        
            
        search_var_organism_include.trace_add("write", lambda *args: self.updating_checkboxes(self.checkboxes_organism_include, self.values_organism_include, search_var_organism_include, organism_scrollframe_include))  
            
        
        
        
        
        logo_organism_exclude = ctk.CTkLabel(organism, text="Exclude Organism", font=ctk.CTkFont(size=20, weight="bold"))
        logo_organism_exclude.grid(row=2, column=0, padx=20, pady=(20, 20))
        
        
        organism_scrollframe_exclude = ctk.CTkScrollableFrame(organism, width=50, corner_radius=2)
        organism_scrollframe_exclude.grid(row=3, column=0, rowspan=1, sticky="nsew")
        organism_scrollframe_exclude.grid_columnconfigure(0, weight=1)
        
        search_var_organism_exclude = tk.StringVar()
        search_entry_organism_exclude = ctk.CTkEntry(organism_scrollframe_exclude, textvariable=search_var_organism_exclude)
        search_entry_organism_exclude.grid(row=0, column=0, padx=10, pady=(5, 5), sticky="w")
        
        checkboxes_organism_exclude = []
        #self.values_organism_exclude = ["organism1", "organism2", "organism3", "organism4"]
        self.values_organism_exclude = [organism[0] for organism in unique_organisms]
        for i, value in enumerate(self.values_organism_exclude):
            state = self.checkbox_states_exclude_dict.get(value, 0)
            
            checkbox_var = tk.IntVar(value=state)
            checkbox = ctk.CTkCheckBox(organism_scrollframe_exclude, text=value, variable=checkbox_var)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(5, 5), sticky="w")
            self.checkboxes_organism_exclude.append(checkbox)
            if (checkbox.get() == 1):
                checkboxes_organism_exclude.append(checkbox.cget("text"))
        
        
        
        search_var_organism_exclude.trace_add("write", lambda *args: self.updating_checkboxes(self.checkboxes_organism_exclude, self.values_organism_exclude, search_var_organism_exclude, organism_scrollframe_exclude)) 
        
    
            
        
        
        
        separator1 = tk.Canvas(newWindowX, width=2, bg='black')
        separator1.grid(row=0, column=1, sticky='ns')
        separator1.create_line(1, 0, 1, separator1.winfo_height())
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        loci = ctk.CTkFrame(newWindowX, width=50, corner_radius=2)
        loci.grid(row=0, column=2, rowspan=4, sticky="nsew")
        loci.grid_columnconfigure(0, weight=1)
        loci.grid_rowconfigure(0, weight = 0)
        loci.grid_rowconfigure(1, weight = 1)
        
        logo_loci_include = ctk.CTkLabel(loci, text="Include Loci", font=ctk.CTkFont(size=20, weight="bold"))
        logo_loci_include.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        loci_scrollframe_include = ctk.CTkScrollableFrame(loci, width=50, corner_radius=2)
        loci_scrollframe_include.grid(row=1, column=0, rowspan=1, sticky="nsew")
        loci_scrollframe_include.grid_columnconfigure(0, weight=1)
        
        search_var_loci_include = tk.StringVar()
        search_entry_loci_include = ctk.CTkEntry(loci_scrollframe_include, textvariable=search_var_loci_include)
        search_entry_loci_include.grid(row=0, column=0, padx=10, pady=(5, 5), sticky="w")
        
        
        unique_loci = db_instance.get_unique_loci()
        
        checkboxes_loci_include = []
        #self.values_loci_include = ["locus1", "locus2", "locus3", "locus4", "locus5", "locus6", "locus7"]
        self.values_loci_include = [loci[0] for loci in unique_loci]
        for i, value in enumerate(self.values_loci_include):
            state = self.checkbox_states_include_dict.get(value, 0)
            
            checkbox_var = tk.IntVar(value=state)
            checkbox = ctk.CTkCheckBox(loci_scrollframe_include, text=value, variable=checkbox_var)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(5, 5), sticky="w")
            self.checkboxes_loci_include.append(checkbox)
            if (checkbox.get() == 1):
                checkboxes_loci_include.append(checkbox.cget("text"))
        
        
        
        search_var_loci_include.trace_add("write", lambda *args: self.updating_checkboxes(self.checkboxes_loci_include, self.values_loci_include, search_var_loci_include, loci_scrollframe_include)) 
        
        
        
        
    
            
            
        logo_loci_exclude = ctk.CTkLabel(loci, text="Exclude loci", font=ctk.CTkFont(size=20, weight="bold"))
        logo_loci_exclude.grid(row=2, column=0, padx=20, pady=(20, 20))
        
        
        loci_scrollframe_exclude = ctk.CTkScrollableFrame(loci, width=50, corner_radius=2)
        loci_scrollframe_exclude.grid(row=3, column=0, rowspan=1, sticky="nsew")
        loci_scrollframe_exclude.grid_columnconfigure(0, weight=1)
        
        
        search_var_loci_exclude = tk.StringVar()
        search_entry_loci_exclude = ctk.CTkEntry(loci_scrollframe_exclude, textvariable=search_var_loci_exclude)
        search_entry_loci_exclude.grid(row=0, column=0, padx=10, pady=(5, 5), sticky="w")
        
        
        checkboxes_loci_exclude = []
        #self.values_loci_exclude = ["locus1", "locus2", "locus3", "locus4", "locus5", "locus6", "locus7"]
        self.values_loci_exclude = [loci[0] for loci in unique_loci]
        for i, value in enumerate(self.values_loci_exclude):
            state = self.checkbox_states_exclude_dict.get(value, 0)
            
            checkbox_var = tk.IntVar(value=state)
            checkbox = ctk.CTkCheckBox(loci_scrollframe_exclude, text=value, variable=checkbox_var)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(5, 5), sticky="w")
            self.checkboxes_loci_exclude.append(checkbox)
            if (checkbox.get() == 1):
                checkboxes_loci_exclude.append(checkbox.cget("text"))
        
        
        
        search_var_loci_exclude.trace_add("write", lambda *args: self.updating_checkboxes(self.checkboxes_loci_exclude, self.values_loci_exclude, search_var_loci_exclude, loci_scrollframe_exclude)) 
        
        
        separator2 = tk.Canvas(newWindowX, width=2, bg='black')
        separator2.grid(row=0, column=3, sticky='ns')
        separator2.create_line(1, 0, 1, separator1.winfo_height())

        
        
        sample = ctk.CTkFrame(newWindowX, width=50, corner_radius=2)
        sample.grid(row=0, column=4, rowspan=4, sticky="nsew")
        sample.grid_columnconfigure(0, weight=1)
        sample.grid_rowconfigure(0, weight = 0)
        sample.grid_rowconfigure(1, weight = 1)
        
        logo_sample_include = ctk.CTkLabel(sample, text="Include sample", font=ctk.CTkFont(size=20, weight="bold"))
        logo_sample_include.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        sample_scrollframe_include = ctk.CTkScrollableFrame(sample, width=50, corner_radius=2)
        sample_scrollframe_include.grid(row=1, column=0, rowspan=1, sticky="nsew")
        sample_scrollframe_include.grid_columnconfigure(0, weight=1)
        
        search_var_sample_include = tk.StringVar()
        search_entry_sample_include = ctk.CTkEntry(sample_scrollframe_include, textvariable=search_var_sample_include)
        search_entry_sample_include.grid(row=0, column=0, padx=10, pady=(5, 5), sticky="w")
        
        
        checkboxes_sample_include = []
        unique_samples = db_instance.get_unique_samples()
        self.values_sample_include = [sample[0] for sample in unique_samples]
        #self.values_sample_include = ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7"]
        for i, value in enumerate(self.values_sample_include):
            state = self.checkbox_states_include_dict.get(value, 0)
            
            checkbox_var = tk.IntVar(value=state)
            checkbox = ctk.CTkCheckBox(sample_scrollframe_include, text=value, variable=checkbox_var)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(5, 5), sticky="w")
            self.checkboxes_sample_include.append(checkbox)
            if (checkbox.get() == 1):
                checkboxes_sample_include.append(checkbox.cget("text"))
        
        
        
        search_var_sample_include.trace_add("write", lambda *args: self.updating_checkboxes(self.checkboxes_sample_include, self.values_sample_include, search_var_sample_include, sample_scrollframe_include))
        
        
        
            
            
        logo_sample_exclude = ctk.CTkLabel(sample, text="Exclude sample", font=ctk.CTkFont(size=20, weight="bold"))
        logo_sample_exclude.grid(row=2, column=0, padx=20, pady=(20, 20))
        
        
        sample_scrollframe_exclude = ctk.CTkScrollableFrame(sample, width=50, corner_radius=2)
        sample_scrollframe_exclude.grid(row=3, column=0, rowspan=1, sticky="nsew")
        sample_scrollframe_exclude.grid_columnconfigure(0, weight=1)
        
        
        search_var_sample_exclude = tk.StringVar()
        search_entry_sample_exclude = ctk.CTkEntry(sample_scrollframe_exclude, textvariable=search_var_sample_exclude)
        search_entry_sample_exclude.grid(row=0, column=0, padx=10, pady=(5, 5), sticky="w")
        
        
        checkboxes_sample_exclude = []
        unique_samples = db_instance.get_unique_samples()
        #self.values_sample_exclude = ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6", "sample7"]
        self.values_sample_exclude = [sample[0] for sample in unique_samples]
        for i, value in enumerate(self.values_sample_exclude):
            state = self.checkbox_states_exclude_dict.get(value, 0)
            
            checkbox_var = tk.IntVar(value=state)
            checkbox = ctk.CTkCheckBox(sample_scrollframe_exclude, text=value, variable=checkbox_var)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(5, 5), sticky="w")
            self.checkboxes_sample_exclude.append(checkbox)
            if (checkbox.get() == 1):
                checkboxes_sample_exclude.append(checkbox.cget("text"))
        
        
        
        search_var_sample_exclude.trace_add("write", lambda *args: self.updating_checkboxes(self.checkboxes_sample_exclude, self.values_sample_exclude, search_var_sample_exclude, sample_scrollframe_exclude))
        
        
        
        print()
        print("The entire dictionary:")
        print()
        print(self.checkbox_states_dict2)
        print()
        
        
        newWindowX.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        
        
        db_instance.closing()
    
        
    def return_values(self, x):
        print()
        print(x)
        # Sample function to return some values
        
    
    def extract_subset2(self):
        print("Extracting")
        
        newWindow = tk.Toplevel(self)
        newWindow.title("Extract subset:")
        newWindow.grid_rowconfigure(0, weight=1)
        newWindow.grid_rowconfigure(1, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        newWindow.grid_columnconfigure(1, weight=1)
        
        organism = ctk.CTkFrame(newWindow, width=50, corner_radius=2)
        organism.grid(row=0, column=0, rowspan=4, sticky="nsew")
        organism.grid_columnconfigure(0, weight=1)
        
        include_organism = ctk.CTkButton(organism, border_width=2, border_color="black", text_color="black", text = "Include Organism", command = self.including_organism)
        include_organism.grid(row=0, column=0, padx=20, pady=(10, 20))
        
        exclude_organism = ctk.CTkButton(organism, border_width=2, border_color="black", text_color="black", text = "Exclude Organism", command = self.choose_project)
        exclude_organism.grid(row=1, column=0, padx=20, pady=(10, 20))
        
        
        
        separator1 = tk.Canvas(newWindow, width=2, bg='black')
        separator1.grid(row=0, column=1, sticky='ns')
        separator1.create_line(1, 0, 1, separator1.winfo_height())
        
        
        
        loci = ctk.CTkFrame(newWindow, width=50, corner_radius=2)
        loci.grid(row=0, column=2, rowspan=4, sticky="nsew")
        loci.grid_columnconfigure(0, weight=1)
        #self.parameters.grid_columnconfigure(1, weight=1)
        #self.parameters.grid_rowconfigure(3, weight=1)
        
        include_loci = ctk.CTkButton(loci, border_width=2, border_color="black", text_color="black", text = "Include Loci", command = self.choose_project)
        include_loci.grid(row=0, column=0, padx=20, pady=(10, 20))
        
        exclude_loci = ctk.CTkButton(loci, border_width=2, border_color="black", text_color="black", text = "Exclude Loci", command = self.choose_project)
        exclude_loci.grid(row=1, column=0, padx=20, pady=(10, 20))
        
        
        separator2 = tk.Canvas(newWindow, width=2, bg='black')
        separator2.grid(row=0, column=3, sticky='ns')
        separator2.create_line(1, 0, 1, separator2.winfo_height())
        
        
        sample = ctk.CTkFrame(newWindow, width=50, corner_radius=2)
        sample.grid(row=0, column=4, rowspan=4, sticky="nsew")
        sample.grid_columnconfigure(0, weight=1)
        #self.parameters.grid_columnconfigure(1, weight=1)
        #self.parameters.grid_rowconfigure(3, weight=1)
        
        include_sample = ctk.CTkButton(sample, border_width=2, border_color="black", text_color="black", text = "Include Sample", command = self.choose_project)
        include_sample.grid(row=0, column=0, padx=20, pady=(10, 20))
        
        exclude_sample = ctk.CTkButton(sample, border_width=2, border_color="black", text_color="black", text = "Exclude Sample", command = self.choose_project)
        exclude_sample.grid(row=1, column=0, padx=20, pady=(10, 20))
        
        
        
    
    def including_organism(self):
        print("Include Organism:")
        
        newWindow = tk.Toplevel(self)
        newWindow.title("Include Organism:")
        newWindow.grid_rowconfigure(0, weight=1)
        #newWindow.grid_rowconfigure(1, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        #newWindow.grid_columnconfigure(1, weight=1)
        
        organism_include = ctk.CTkScrollableFrame(newWindow, width=50, corner_radius=2)
        organism_include.grid(row=0, column=0, rowspan=3, sticky="nsew")
        organism_include.grid_columnconfigure(0, weight=1)
        
        
        organism_include.grid_rowconfigure(0, weight=1)
        organism_include.grid_rowconfigure(1, weight=1)
        organism_include.grid_columnconfigure(0, weight=1)
        organism_include.grid_columnconfigure(1, weight=1)
        
        checkboxes = []
        
        organism_values = ["organism1", "organism2", "organism3", "organism4"]
        for i, value in enumerate(organism_values):
            checkbox = ctk.CTkCheckBox(organism_include, text=value)
            checkbox.grid(row=i, column=0, padx=10, pady=(5, 5), sticky="w")
            checkboxes.append(checkbox)
        
        
    def change_appearance_mode_event(self, new_appearance_mode: str):
        ctk.set_appearance_mode(new_appearance_mode)
        
        
    def parse_parameter_file(self):
        print("Parsing:")
        
        calc_param_list = []
        param_list_dirs = []
        param_list_files = []
        param_list_special = []
        
        paramter_file_path = os.path.normpath(os.getcwd() + '/Parameter.txt')
        if (os.path.exists(paramter_file_path)):
            with open(paramter_file_path) as parameters:
                for line in parameters:
                    if (not line) or line.startswith("#"):
                        continue
                    if (line.startswith("DataDir")):
                        raw_data_dir = line.split("=")[1].rstrip("\r\n").strip()
                        if (not raw_data_dir):
                            raise Exception("No name for data_dir in parameter file!")
                        param_list_dirs.append(os.path.basename(os.path.normpath(raw_data_dir)))
                        continue
                    if (line.startswith("PrimerFile")):
                        primerfile = line.split("=")[1].rstrip("\r\n").strip()
                        if (not primerfile):
                            raise Exception("No name for primerfile in parameter file!")
                        param_list_files.append(os.path.basename(os.path.normpath(primerfile)))
                        continue
                    if (line.startswith("SampleFile")):
                        sample_sheet = line.split("=")[1].rstrip("\r\n").strip()
                        if (not sample_sheet):
                            raise Exception("No name for sample file in parameter file!")
                        param_list_files.append(os.path.basename(os.path.normpath(sample_sheet)))
                        continue
                    if (line.startswith("WorkingDir")):
                        output_folder = line.split("=")[1].rstrip("\r\n").strip()
                        if (not output_folder):
                            output_folder = "Output_default"
                        param_list_dirs.append(output_folder)
                        continue
                    if (line.startswith("Ploidy")):
                        ploidy = line.split("=")[1].rstrip("\r\n").strip()
                        if (not ploidy):
                            ploidy = "diploid"
                        param_list_special.append(ploidy)
                        continue
                    
                    if (line.startswith("numeberMM")):
                        max_mismatches = line.split("=")[1].rstrip("\r\n").strip()
                        if (not max_mismatches):
                            max_mismatches = 2
                        else:
                            max_mismatches = int(max_mismatches)
                        calc_param_list.append(max_mismatches)
                        continue
                    if (line.startswith("minCount")):
                        minCount = line.split("=")[1].rstrip("\r\n").strip()
                        if (not minCount):
                            minCount = 10
                        else:
                            minCount = int(minCount)
                        calc_param_list.append(minCount)
                        continue
                    if (line.startswith("minLength")):
                        min_length = line.split("=")[1].rstrip("\r\n").strip()
                        if (not min_length):
                            min_length = 290
                        else:
                            min_length = int(min_length)
                        calc_param_list.append(min_length)
                        continue
                    if (line.startswith("consThreashold")):
                        alpha = line.split("=")[1].rstrip("\r\n").strip()
                        if (not alpha):
                            alpha = 0.7
                        else:
                            alpha = float(alpha)
                        calc_param_list.append(alpha)
                        continue
                    if (line.startswith("LengthWindow")):
                        values = line.split("=")[1].rstrip("\r\n").strip().split(",")
                        if (not values):
                            lengthWindow = ['310', '600']
                        else:
                            lengthWindow = []
                            lengthWindow.append(values[0])
                            lengthWindow.append(values[1])
                        calc_param_list.append(lengthWindow)
                        continue
                    if (line.startswith("AlleleList")):
                        allele_list = line.split("=")[1].rstrip("\r\n").strip()
                        param_list_files.append(os.path.normpath(allele_list))
                        continue
                    
                    if (line.startswith("operating_system")):
                        par_os = line.split("=")[1].rstrip("\r\n").strip()
                        if (not par_os):
                            par_os = "linux"
                        param_list_special.append(par_os.lower())
                        continue
                    if (line.startswith("unique_identifier")):
                        unique_name = line.split("=")[1].rstrip("\r\n").strip()
                        if (not unique_name):
                            unique_name = "Default"
                        param_list_special.append(unique_name)
                        continue
                    if (line.startswith("ReferenceFile")):
                        Reference = line.split("=")[1].rstrip("\r\n").strip()
                        if (not Reference):
                            Reference = ""
                        param_list_files.append(os.path.normpath(Reference))
                        continue
                    if (line.startswith("RExecutable")):
                        Rexecutable = line.split("=")[1].rstrip("\r\n").strip()
                        if (not Rexecutable):
                            Rexecutable = "R"
                        param_list_files.append(os.path.normpath(Rexecutable))
                        continue
        else:
            print("There is no paramter file with the name Parameter.txt")
            param_list_dirs = ["", ""]
            param_list_files = ["", "", "", "", ""]
            calc_param_list = ["", "", "", "", ""]
            param_list_special = ["", "", ""]
        
        return param_list_dirs, param_list_files, param_list_special, calc_param_list
    
    
    def get_general_params(self):
        
        newWindow = tk.Toplevel(self)
        newWindow.transient(self)
        newWindow.title("General Parameters:")
        newWindow.grid_rowconfigure(0, weight=0)
        newWindow.grid_rowconfigure(1, weight=1)
        newWindow.grid_rowconfigure(2, weight=0)
        newWindow.grid_rowconfigure(3, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        newWindow.grid_columnconfigure(1, weight=1)
        newWindow.grid_columnconfigure(2, weight=1)
        
        logo_dir_params = ctk.CTkLabel(newWindow, text="Folders:", font=ctk.CTkFont(size=20, weight="bold"))
        logo_dir_params.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        dir_params = ctk.CTkFrame(newWindow, width=100, corner_radius=2)
        dir_params.grid(row=1, column=0, rowspan=1, sticky="nsew")
        dir_params.grid_columnconfigure(0, weight=1)
        dir_params.grid_columnconfigure(1, weight=1)
        dir_params.grid_columnconfigure(2, weight=1)
        dir_params.grid_rowconfigure(0, weight=1)
        dir_params.grid_rowconfigure(1, weight=1)
        
        
        logo_file_params = ctk.CTkLabel(newWindow, text="Files:", font=ctk.CTkFont(size=20, weight="bold"))
        logo_file_params.grid(row=2, column=0, padx=20, pady=(20, 20))
        
        file_params = ctk.CTkFrame(newWindow, width=100, corner_radius=2)
        file_params.grid(row=3, column=0, rowspan=1, sticky="nsew")
        file_params.grid_columnconfigure(0, weight=1)
        file_params.grid_columnconfigure(1, weight=1)
        file_params.grid_columnconfigure(2, weight=1)
        file_params.grid_rowconfigure(0, weight=1)
        file_params.grid_rowconfigure(1, weight=1)
        file_params.grid_rowconfigure(2, weight=1)
        file_params.grid_rowconfigure(3, weight=1)
        
        
        logo_calc_params = ctk.CTkLabel(newWindow, text="Calculation Params:", font=ctk.CTkFont(size=20, weight="bold"))
        logo_calc_params.grid(row=0, column=1, padx=20, pady=(20, 20))
        
        calc_params = ctk.CTkFrame(newWindow, width=100, corner_radius=2)
        calc_params.grid(row=1, column=1, rowspan=1, sticky="nsew")
        calc_params.grid_columnconfigure(0, weight=1)
        calc_params.grid_columnconfigure(1, weight=1)
        calc_params.grid_rowconfigure(0, weight=1)
        calc_params.grid_rowconfigure(1, weight=1)
        calc_params.grid_rowconfigure(2, weight=1)
        calc_params.grid_rowconfigure(3, weight=1)
        calc_params.grid_rowconfigure(4, weight=1)
        
        
        logo_special_params = ctk.CTkLabel(newWindow, text="Additional Params:", font=ctk.CTkFont(size=20, weight="bold"))
        logo_special_params.grid(row=2, column=1, padx=20, pady=(20, 20))
        
        special_params = ctk.CTkFrame(newWindow, width=100, corner_radius=2)
        special_params.grid(row=3, column=1, rowspan=1, sticky="nsew")
        special_params.grid_columnconfigure(0, weight=1)
        special_params.grid_columnconfigure(1, weight=1)
        special_params.grid_rowconfigure(0, weight=1)
        special_params.grid_rowconfigure(1, weight=1)
        special_params.grid_rowconfigure(2, weight=1)
        
        
        entries_dirs = []
        entries_files = []
        entries_calc = []
        entries_special = []
        entries_total = []
        
        param_list_dirs = ["RawdataDir", "OutputDir"]
        number1 = 0
        for i, param in enumerate(param_list_dirs):
            ctk.CTkLabel(dir_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(dir_params, corner_radius=5) #, textvariable=text
            if not((param == "OutputDir")):
                browse_button = tk.Button(dir_params, text="Browse", command=lambda entry = param_entry: self.browse_dir(entry))
                browse_button.grid(column=2, row=i, padx=5, pady=5)
            if (len(self.var_list_dir) > 0):
                param_entry.insert(tk.END, self.var_list_dir[i])
            else:
                param_entry.insert(tk.END, self.param_list_dirs[i])
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            entries_dirs.append(param_entry)
            number1 += 1
        entries_total.append(entries_dirs)
            
        
        
        param_list_files = ["Primers", "SampleSheet", "AlleleList"]
        number2 = 0
        for i, param in enumerate(param_list_files):
            ctk.CTkLabel(file_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(file_params, corner_radius=5) #, textvariable=text
            
            browse_button = tk.Button(file_params, text="Browse", command=lambda entry=param_entry: self.browse_file_folder(entry))
            browse_button.grid(column=2, row=i, padx=5, pady=5)
            if (len(self.var_list_file) > 0):
                param_entry.insert(tk.END, os.path.basename(os.path.normpath(self.var_list_file[i])))
            else:
                param_entry.insert(tk.END, os.path.basename(os.path.normpath(self.param_list_files[i])))
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            entries_files.append(param_entry)
            number2 += 1
            
        ctk.CTkLabel(file_params, text="Reference:", width=20).grid(column=0, row=number2, padx=5, pady=5)
        param_entry = ctk.CTkEntry(file_params, corner_radius=5) #, textvariable=text
        browse_button = tk.Button(file_params, text="Browse", command=lambda entry=param_entry: self.browse_file_folder(entry))
        browse_button.grid(column=2, row=number2, padx=5, pady=5)
        if (len(self.var_list_file) > 0):
            param_entry.insert(tk.END, os.path.normpath(self.var_list_file[number2]))
        else:
            param_entry.insert(tk.END, self.param_list_files[number2])
        param_entry.grid(column=1, row=number2, padx=5, pady=5)
        #param_entry.insert(0, default_values[i])
        entries_files.append(param_entry)
        number2 += 1
        
        ctk.CTkLabel(file_params, text="RExecutable:", width=20).grid(column=0, row=number2, padx=5, pady=5)
        param_entry = ctk.CTkEntry(file_params, corner_radius=5) #, textvariable=text
        browse_button = tk.Button(file_params, text="Browse", command=lambda entry=param_entry: self.browse_file_folder(entry))
        browse_button.grid(column=2, row=number2, padx=5, pady=5)
        if (len(self.var_list_file) > 0):
            param_entry.insert(tk.END, os.path.normpath(self.var_list_file[number2]))
        else:
            param_entry.insert(tk.END, os.path.normpath(self.param_list_files[number2]))
        param_entry.grid(column=1, row=number2, padx=5, pady=5)
        #param_entry.insert(0, default_values[i])
        entries_files.append(param_entry)
        number2 += 1
        
        entries_total.append(entries_files)
        
        
        
        param_list_calc = ["max_mismatches", "minLength", "lengthWindow", "consThreashold", "minCount"]
        placeholder = ["2", "290", "310 600", "0.7", "20"]
        number3 = 0
        for i, param in enumerate(param_list_calc):
            ctk.CTkLabel(calc_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(calc_params, corner_radius=5) #, textvariable=text
            
            if (len(self.var_list_calc) > 0):
                param_entry.insert(tk.END, self.var_list_calc[i])
            else:
                if (self.calc_param_list[i] == ""):
                    param_entry.insert(tk.END, placeholder[i])
                else:
                    param_entry.insert(tk.END, self.calc_param_list[i])
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            entries_calc.append(param_entry)
            number3 += 1
        entries_total.append(entries_calc)
        
        param_list_special = ["Ploidy", "OperatingSystem", "UniqueID"]
        placeholder_special = ["diploid", "windows", "default"]
        number4 = 0
        for i, param in enumerate(param_list_special):
            ctk.CTkLabel(special_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(special_params, corner_radius=5) #, textvariable=text
            
            if (len(self.var_list_special) > 0):
                param_entry.insert(tk.END, self.var_list_special[i])
            else:
                if (self.param_list_special[i] == ""):
                    param_entry.insert(tk.END, placeholder_special[i])
                else:
                    param_entry.insert(tk.END, self.param_list_special[i])
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            entries_special.append(param_entry)
            number4 += 1
        entries_total.append(entries_special)
        
        update_params_button = tk.Button(newWindow, text="Write ParameterFile", command = lambda x = entries_total : self.write_params(x))
        update_params_button.grid(column=1, row=4, sticky="se", padx=15, pady=4)
        
        #, command = lambda x = entries_total : self.set_general_params(x)
        write_params_button = tk.Button(newWindow, text="Update Params", command = lambda x = entries_total : self.set_params(x))
        write_params_button.grid(column=0, row=4, sticky="sw", padx=15, pady=4)
        
    
    #directory_path = filedialog.askdirectory(initialdir=self.workspace_location, title="Select folder")
    
    def browse_dir(self, entry):
        directory_path = filedialog.askdirectory(initialdir=self.workspace_location, title="Select folder")
        directory_path = os.path.normpath(directory_path)
        if (self.os_param.lower() == "linux"):
            output_folder = directory_path.split('/')[-1]
            output_path = directory_path
        elif (self.os_param.lower() == "windows"):
            output_folder = directory_path.split('\\')[-1]
            output_path = directory_path
        else:
            new_window = tk.Toplevel(self)
            new_window.title("Operating System!")
            long_text = f"FIRST CHOOSE WHAT OPERATING SYSTEM YOU ARE ON! \n\n"
            text_widget = tk.Text(new_window, wrap=tk.WORD)
            scrollbar = tk.Scrollbar(new_window, command=text_widget.yview)
            text_widget.config(yscrollcommand=scrollbar.set)
            
            text_widget.insert(tk.END, long_text)
            text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=20, pady=20)
            scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
            
        entry.delete(0, tk.END)  # Clear the entry
        entry.insert(0, output_folder)  # Insert the file/folder path into the entry
        
    
    def browse_file_folder(self, entry):
        # Check if you want to select a file or folder (for simplicity, I'm using file selection here)
        file_path = filedialog.askopenfilename(initialdir=self.workspace_location, title="Select folder")
        file_path = os.path.normpath(file_path)
        if (self.os_param.lower() == "linux"):
            output_folder = file_path.split('/')[-1]
            output_path = file_path
        elif (self.os_param.lower() == "windows"):
            output_folder = file_path.split('\\')[-1]
            output_path = file_path
        if file_path:
            entry.delete(0, tk.END)  # Clear the entry
            entry.insert(0, output_path)  # Insert the file/folder path into the entry
    
    def checkbox_toggle(self):
        if (self.import_check.get() == 1):
            self.textbox_bool = True
            print("Set to True")
        else:
            self.textbox_bool = False
            print("Set to false")
            
        
    def set_params(self, entries):
        print("Setting them.")
        
        list_dir_params = entries[0]
        list_file_params = entries[1]
        list_calc_params = entries[2]
        list_special_params = entries[3]
        
        var_params_dir = []
        for i, entry in enumerate(list_dir_params):
            var_params_dir.append(entry.get())
        
        var_params_file = []
        for i, entry in enumerate(list_file_params):
            var_params_file.append(entry.get())
        
        var_params_calc = []
        for i, entry in enumerate(list_calc_params):
            var_params_calc.append(entry.get())
        
        var_params_special = []
        for i, entry in enumerate(list_special_params):
            var_params_special.append(entry.get())
        
        self.rawdata_dir = var_params_dir[0]
        self.outputfolder = var_params_dir[1]
        
        self.primerfile = var_params_file[0]
        self.sample_sheet = var_params_file[1]
        self.allele_list_path = var_params_file[2]
        self.reference_file = var_params_file[3]
        self.r_script_path = var_params_file[4]
        
        self.max_mismatches = var_params_calc[0]
        self.min_length = var_params_calc[1] 
        self.lengthWindow = var_params_calc[2]
        self.alpha = var_params_calc[3]
        self.minCount = var_params_calc[4]
        
        self.ploidy = var_params_special[0]
        self.os_param = var_params_special[1]
        self.unique_prefix_name = var_params_special[2]
        
        self.var_list_dir = []
        self.var_list_file = []
        self.var_list_calc = []
        self.var_list_special = []
    
        self.var_list_dir.append(self.rawdata_dir)
        self.var_list_dir.append(self.outputfolder)
        
        self.var_list_file.append(self.primerfile)
        self.var_list_file.append(self.sample_sheet)
        self.var_list_file.append(self.allele_list_path)
        self.var_list_file.append(self.reference_file)
        self.var_list_file.append(self.r_script_path)
        
    
        self.var_list_calc.append(self.max_mismatches)
        self.var_list_calc.append(self.min_length)
        self.var_list_calc.append(self.lengthWindow)
        self.var_list_calc.append(self.alpha)
        self.var_list_calc.append(self.minCount)
        
        self.var_list_special.append(self.ploidy)
        self.var_list_special.append(self.os_param)
        self.var_list_special.append(self.unique_prefix_name)
        
        print()
        print(self.rawdata_dir)
        print(self.primerfile)
        print(self.sample_sheet)
        print(self.outputfolder)
        print(self.ploidy)
    
       
    def write_params(self, entries):
        
        self.set_params(entries)
        
        with open("Parameter.txt", 'w') as file:
            file.write("###Parameters for Amplicon_Length_Allele_Call \n")
            file.write("#Diretory containing raw fastq/fastq.gz files \n")
            file.write("DataDir = " + os.path.normpath(self.rawdata_dir) + "\n")
            file.write("#File containing primer infromation \n")
            file.write("PrimerFile = " + os.path.basename(os.path.normpath(self.primerfile)) + " \n")
            file.write('#File containing sample correspondance \n')
            file.write("SampleFile = " + os.path.basename(os.path.normpath(self.sample_sheet)) + "\n")
            file.write('#Maximum number of mismatches between primer and amplicon sequences (integer; default=2) \n')
            file.write("numeberMM = " + str(self.max_mismatches) + " \n")
            file.write('#Minimum amplicon length to be considered (integer; default=250) \n')
            file.write("minLength= " + str(self.min_length) + " \n")
            file.write('#Window interval to be shown in the amplicon length plots (format: minimum_lengths,maximum_length; default=300,600) \n')
            self.lengthWindowx = self.lengthWindow.split()
            file.write("LengthWindow= " + str(self.lengthWindowx[0]) + "," + str(self.lengthWindowx[1]) + "\n")
            file.write("\n")
            file.write("###General parameters \n")
            file.write("#Ploidy level of the data (haploid/diploid; default = diploid) \n")
            file.write("Ploidy= " + str(self.ploidy) + "\n")
            file.write('#Directory to save all output files (default = "OutDir_" + Current_Date_and_Time) \n')
            file.write("WorkingDir = \n")
            file.write("\n")
            file.write('#Matrix conating the selected length based alleles (default = WorkingDir + "/MarkerPlots/markermatrix.csv") \n')
            file.write("LengthAllelels= \n")
            file.write('#Directory containg demultiplexed (default = WorkingDir + "/SeparatOut") \n')
            file.write("SeparatOut = \n")
            file.write('#Consensus threashold (float; default = 0.7) \n')
            file.write("consThreashold = " + str(self.alpha) + " \n")
            file.write('#Minimum number of read count for an allele to be considered (integer; default = 20) \n')
            file.write("minCount = " + str(self.minCount) + " \n")
            file.write('#List of allleles priduced in previous run (if left empty it will call the alleles denovo) \n')
            file.write("AlleleList = " + os.path.normpath(self.allele_list_path) + " \n")
            file.write('#Prefix to save the output matrix and the allele listfiles (default =  "Final_" + Current_Date_and_Time) \n')
            file.write("out_pre = \n")
            file.write("\n")
            file.write("operating_system = " + str(self.os_param) + " \n")
            file.write("unique_identifier = " + str(self.unique_prefix_name) + " \n")
            file.write("ReferenceFile = " + str(self.reference_file) + " \n")
            file.write("RExecutable = " + os.path.normpath(self.r_script_path) + " \n")
            file.write("\n")
        
    def choose_project(self):
        print("Setting Project")
        
        newWindow = tk.Toplevel(self)
        newWindow.title("Calculation Parameters:")
        newWindow.grid_rowconfigure(0, weight=1)
        #newWindow.grid_rowconfigure(1, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        #newWindow.grid_columnconfigure(1, weight=1)
        
        self.logo_pipeline = ctk.CTkLabel(newWindow, text="Choose a ProjectName", font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_pipeline.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        project_entry = ctk.CTkEntry(newWindow)
        project_entry.grid(row = 1, column = 0, padx=5, pady=(20, 20))
        
        enter = ctk.CTkButton(newWindow, border_width=1, border_color="black", text_color="black", text = "Enter", command = lambda x = project_entry : self.set_project_name(x))
        enter.grid(row=2, column=0, padx=20, pady=(20, 10))
        
    def set_project_name(self, x):
        print("Project Name: " + x.get())
        #self.textbox.insert('end-1c', "Project Name: " + x.get() + "\n")
        
        project_path = os.path.normpath(os.getcwd() + "/" + x.get())
        if os.path.exists(project_path):
            newWindow = tk.Toplevel(self)
            textmess = ctk.CTkLabel(newWindow, text="This project already exists. Saving Outputs in existing Project", font=ctk.CTkFont(size=20, weight="bold"))
            textmess.grid(row=0, column=0, padx=20, pady=(20, 20))
            
            self.projectfolder = x.get()
            self.projectfolder_path = os.path.normpath(self.workspace_location + "/" + self.projectfolder)
            #self.textbox.insert('end-1c', "Project Name: " + x.get() + "\n")
            #self.primerfile = primerfile
            #print("Primerfile: " + self.primerfile + "\n")
        else:
            os.mkdir(os.path.normpath(os.getcwd() + "/" + x.get()))
            self.textbox.insert('end-1c', "Project Name: " + x.get() + "\n")
            self.projectfolder = x.get()
            self.projectfolder_path = os.path.normpath(self.workspace_location + "/" + self.projectfolder)
        
        
    def run_first_script(self):
        
        self.textbox.delete(0.0, 'end')
        self.textbox.insert('end-1c', "Start First Script! \n")
        
        thread = threading.Thread(target=self.run_amplicon_methods)
        thread.start()
    
        # Start appending periods to the textbox
        self.append_period()
        
        
        
    def run_amplicon_methods(self):
        
        self.textbox.delete(0.0, 'end')
        
        self.textbox.insert('end-1c', "Start First Script! \n")
        Amplicon_inst = amplicon_class(self.projectfolder, self.outputfolder, self.primerfile, self.rawdata_dir, self.sample_sheet, self.max_mismatches, self.min_length, self.minCount, self.alpha, self.lengthWindow, self.ploidy, self.os_param)
        project, output, samplesheet, primerfile, rawdata = Amplicon_inst.return_params()
        
        #global is_working
        self.is_working = True
        
        self.textbox.insert('end-1c', "\n")
        self.textbox.insert('end-1c', "ProjectName : " + project + " \n")
        self.textbox.insert('end-1c', "OutputFolder : " + output + " \n")
        self.textbox.insert('end-1c', "SampleSheet : " + samplesheet + " \n")
        self.textbox.insert('end-1c', "Primerfile : " + primerfile + " \n")
        self.textbox.insert('end-1c', "RawDataFolder : " + rawdata + " \n")
        self.textbox.insert('end-1c', "\n")
        
        self.append_period()
        
        Amplicon_inst.set_executables()
        #Amplicon_inst.print_executables()
        
        Amplicon_inst.set_samplelist()
        Amplicon_inst.checkInputDir()
        print()
        self.textbox.insert('end-1c', "Trimming started \n")
        errorcode = Amplicon_inst.runTrimomatic() #trimmed fastq files in /QC/
        print()
        message = "Success"
        if (errorcode != 0):
            message = "Failure"
        self.textbox.insert('end-1c', "Trimming sequences finished with " + message + " \n")
        
        print()
        self.textbox.insert('end-1c', " \n")
        self.textbox.insert('end-1c', "Merging started \n")
        errorcode = Amplicon_inst.runUsearchMergePairs() #merged fastq files in /MergedOut/
        print()
        message = "Success"
        if (errorcode != 0):
            message = "Failure"
        self.textbox.insert('end-1c', "Merging sequences finished with " + message + " \n")
        #self.textbox.insert('end-1c', "Merging sequences finished with \n")
        
        print()
        self.textbox.insert('end-1c', "Demultiplexing started \n")
        count = Amplicon_inst.runDemultiplex() #demultiplexed fasta files in /Separatout/ with name convention: newSampleName - Locus/PrimerName
        print()
        self.textbox.insert('end-1c', "Demultiplexing finished with " + str(count) + " demultiplexed files \n")
        
        Amplicon_inst.getLengthStatistics()
        Amplicon_inst.write_to_analysis_lengths()
        
        print()
        self.textbox.insert('end-1c', "Finding Length-based Genotypes started \n")
        error_code = Amplicon_inst.GenotypeLength(self.r_script_path)
        print()
        message = "Success"
        if (errorcode != 0):
            message = "Failure"
        self.textbox.insert('end-1c', "Getting Genotype Lengths finished with " + message + " \n")
        #lengths = Amplicon_inst.get_unique_lengths()
        #self.textbox.insert('end-1c', "Found  " + str(lengths) + " unique lengths. \n")#get_unique_lengths
        #self.textbox.insert('end-1c', "Getting Genotype Lengths finished with \n")
        
        
        
        self.is_working = False
        
    def append_period2(self):
        # Append a period to the Text widget
        self.textbox2.insert('end-1c', ".")
        
        # If the working function is still running, schedule this function to run again
        if self.is_working2:
            self.after(1000, self.append_period2)
        
    def append_period(self):
        # Append a period to the Text widget
        self.textbox.insert('end-1c', ".")
        
        # If the working function is still running, schedule this function to run again
        if self.is_working:
            self.after(1000, self.append_period)
        
        
    
    def run_second_script(self):
        
        self.textbox.delete(0.0, 'end')
        self.textbox.insert('end-1c', "Start Second Script! \n")
        
        thread = threading.Thread(target=self.run_sequence_methods)
        thread.start()
    
        # Start appending periods to the textbox
        self.append_period()
    
    
    
    def run_sequence_methods(self):
        
        self.textbox.delete(0.0, 'end')
        self.textbox.insert('end-1c', "Start Second Script! \n")
        
        
        Sequence_inst = class_quality(self.projectfolder, self.outputfolder, self.primerfile, self.sample_sheet, self.minCount, self.alpha, self.ploidy, self.unique_prefix_name, self.os_param, self.allele_list_path, self.textbox_bool)
        #Sequence_inst.print_input_paths()
        self.textbox.insert('end-1c', " \n")
        
        #global is_working
        self.is_working = True
        
        print()
        self.textbox.insert('end-1c', "Get all Sample-Locus combo for each Genotype Length... \n")
        print()
        #self.textbox.insert('end-1c', " \n")
        allele_length_dict = Sequence_inst.get_Marker_Sample_Alleles_length_diploid()
        print()
        for sample_locus, length in allele_length_dict.items():
            for l, count in length.items():
                self.textbox.insert('end-1c', sample_locus + ":  \n")
                self.textbox.insert('end-1c', "length : " + str(l) + ", count : "  + str(count) + " \n")
        
        self.textbox.insert('end-1c', " \n")
        self.textbox.insert('end-1c', "Making ConsensusSequences... \n")
        Sequence_inst.RunConsensusAll()
        print()
        Sequence_inst.joinSamplesSameMarker()
        print()
        self.textbox.insert('end-1c', " \n")
        self.textbox.insert('end-1c', "Correcting Sequences... \n")
        Sequence_inst.correctAllSeqs()
        print()
        
        self.textbox.insert('end-1c', " \n")
        self.textbox.insert('end-1c', "Calling Alleles... \n")
        Sequence_inst.AlleleCall()
        if (self.allele_list_path == "" or os.path.basename(os.path.normpath(self.allele_list_path)) == "." or not(os.path.exists(self.allele_list_path))):
            self.textbox.insert('end-1c', "No AlleleList, Denovo Call!\n")
            self.textbox.insert('end-1c', "\n")
        else:
            self.textbox.insert('end-1c', "ALleleList is there, ListBased Call!\n")
        
        if (self.textbox_bool):
            self.textbox.insert('end-1c', "Output written into Genalex Format \n")
        else:
            self.textbox.insert('end-1c', "Output NOT written into Genalex Format \n")
        
        self.is_working = False


if __name__ == "__main__":
    main()