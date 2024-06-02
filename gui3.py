# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 11:42:37 2023

@author: Sebastian
"""

import tkinter as tk
from tkinter import ttk
import tkinter.messagebox
import customtkinter as ctk
import os, re, csv, subprocess, platform, time
import multiprocessing
from multiprocessing import Process, set_start_method, freeze_support

from tkinter import filedialog 

from lengths import length_class
from quality import class_quality
from database_class import SQdatabase

import threading, codecs





ctk.set_appearance_mode("System")  # Modes: "System" (standard), "Dark", "Light"
ctk.set_default_color_theme("blue")

def main():
    
    # try:
    #     set_start_method('spawn')
    # except RuntimeError:
    #     pass
    # print("Testing")
    freeze_support()
    p = Process(target=run_other_process_task)
    p.start()
    
    root = root_window()
    #frame = ssr_window(root)
    root.mainloop()
    
    p.join()
    
    return 0

def run_other_process_task():
    # Function to be run by multiprocessing
    print("This runs in a separate process.")


class root_window(ctk.CTk):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.title("SSR GUI")
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
        
        self.workspace = os.getcwd()
        self.build_mainframe()
        
        #["PrimerFile", "SampleFile", "RExecutable", "AlleleList"]
        #["MaxMismatch", "MinCount", "MinLength", "ConsensusThreshold", "lengthWindow"]
        #["Ploidy", "Operatingsystem", "Uniqueidentifier", "Indexcomboposition"]
        
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
            "SampleDuplicates" : False,
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
            },
            "Reference" : ""
        }
        
        self.dict_check = {"OutputFolder" : False, "Bin" : False, "RawData" : False, "PrimerFile" : False, "SampleFile" : False, "Reference" : False, "RExecutable" : False}
        
        
        self.parse_params()
        
        if (self.params["Operatingsystem"] != "windows" or self.params["Operatingsystem"] != "linux"):
            if (os.name == "nt"):
                self.params["Operatingsystem"] = "windows"
            elif (os.name == "posix"):
                self.params["Operatingsystem"] = "linux"
            else:
                self.params["Operatingsystem"] = "NoLinuxNoWindows"
        
        
        self.sample_number = 0
        
        #Parsing relevant data:
        
        
        
        
        
        
        
        
        
        
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
        
        
        self.output_dir = self.param_list_dirs[0]
        self.bin_dir = self.param_list_dirs[1]
        self.rawdata_dir = self.param_list_dirs[2]
        
        self.primerfile = self.param_list_files[0]
        self.sample_sheet = self.param_list_files[1]
        self.Rexe = self.param_list_files[2]
        self.allelelist = self.param_list_files[3]
        
        self.ploidy = self.param_list_special[0]
        self.os_param = self.param_list_special[1]
        self.unique_prefix_name = self.param_list_special[2]
        self.indexcombo_position = self.param_list_special[3]
        if (self.os_param == "" or not(self.os_param == "windows") or not(self.os_param == "linux")):
            self.os_name = platform.system()
            if (self.os_name == 'Windows'):
                self.os_param = "windows"
            elif (self.os_name == 'Linux'):
                self.os_param = "linux"
        if (self.ploidy == "" or not(self.ploidy == "diploid") or not(self.os_param == "haploid")):
            self.ploidy = "diploid"
        
        self.param_list_special[0] = self.ploidy
        self.param_list_special[1] = self.os_param
        
        print(self.os_param)
        
        self.max_mismatches = self.calc_param_list[0]
        self.minCount = self.calc_param_list[1]
        self.min_length = self.calc_param_list[2]   
        self.alpha = self.calc_param_list[3]
        self.lengthWindow = self.calc_param_list[4]
        print(self.lengthWindow)
        
        #self.r_script_path = os.path.normpath("C:\Program Files\R\R-4.3.1\bin\Rscript.exe")
        
        self.is_working = False
        self.is_working2 = False
        self.is_working_advanced1 = False
        self.is_working_advanced2 = False
        
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
        
        self.checkbox_advanced = []
        
        #self.reference_file_name = "reference.csv"
        #self.reference_file_path = os.path.normpath(self.workspace_location + "/" + self.reference_file_name)
        
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
        
        self.sample_names = []
        
        self.sample_params_import = ["3", "2", "10"]
        self.sequence_params_import = ["11", "2", "10"]
        self.samplesheetname = ""
        
        self.reference_dict = {}
        self.reference_file = ""
        self.columnr_params, self.rownr_params = self.parse_reference_parameterfile()
        
        print("ColumnPositions:" + "\n")
        print(self.columnr_params)
        print("RowPositions:" + "\n")
        print(self.rownr_params)
        print("Referencefilepath:")
        print(self.reference_file)
        print()
        
        self.adapterfile = "TrueSeqAdaptersInUsage.fa"
        
        self.check_var_advanced_pipeline_settings = [tk.IntVar(value=0) for _ in range(5)]
        
        #self.checkbox_states_pipeline_advanced = {}
        self.advanced_pipeline_options = ["Trimmomatic", "Usearch", "Demultiplexing", "Markerstatistics", "Markerplots+Markermatrix",  "Extract AlleleLengths", "ConsensusSequence", "ConsensusTogether", "Corrected", "AlleleCall"]
        self.checkbox_states_pipeline_advanced = {key : 0 for key in self.advanced_pipeline_options}
        self.checkboxes_pipeline1 = []
        self.checkboxes_pipeline2 = []
        
        self.report_dict = {}
        self.report_dict["Trimmomatic"] = []
        self.report_dict["Usearch"] = []
        self.report_dict["Demultiplexing"] = []
        self.report_dict["MarkerStatistics"] = []
        self.report_dict["Markerplots+Matrix"] = []
        self.report_dict["AllelesOut"] = []
        self.report_dict["ConsensusOut"] = []
        self.report_dict["ConsensusTogether"] = []
        self.report_dict["Corrected"] = []
        self.report_dict["AlleleCall"] = []
        
        
        self.newWindow_advanced = None
        
        self.performance = True
        
        self.checkbox_performance = 1
        
    def build_mainframe(self):
        self.columnconfigure(0, weight=1) 
        self.columnconfigure(1, weight=0)
        self.columnconfigure(2, weight=1)
        self.columnconfigure(3, weight=0)
        self.columnconfigure(4, weight=1)
        self.columnconfigure(5, weight=0)
        self.columnconfigure(6, weight=0)
        self.rowconfigure(0, weight=1)
        
        
        self.parameters = ctk.CTkFrame(self, width=100, corner_radius=2)
        self.parameters.grid(row=0, column=0, rowspan=6, sticky="nsew")
        self.parameters.grid_columnconfigure(0, weight=1)
        self.parameters.grid_rowconfigure(5, weight=1)
        
        self.logo_parameters = ctk.CTkLabel(self.parameters, text="Preparation", font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_parameters.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        self.reference_import_param = ctk.CTkButton(self.parameters, border_width=1, border_color="black", text_color="black", text = "Reference Parameters", command = self.set_reference_import_params)
        self.reference_import_param.grid(row=1, column=0, padx=20, pady=(10, 10))
        
        self.general_param = ctk.CTkButton(self.parameters, border_width=1, border_color="black", text_color="black", text = "Pipeline Parameters", command = self.get_general_params)
        self.general_param.grid(row=2, column=0, padx=20, pady=(10, 10))
        
        self.genalex_param = ctk.CTkButton(self.parameters, border_width=1, border_color="black", text_color="black", text = "Genalex Parameter", command = self.get_genalex_params)
        self.genalex_param.grid(row=3, column=0, padx=20, pady=(10, 10))
        
        self.checkit = ctk.CTkButton(self.parameters, border_width=1, border_color="black", text_color="black", text = "Workspace Status", command = self.check_status)
        self.checkit.grid(row=4, column=0, padx=20, pady=(10, 10))
        
        self.doc_param = ctk.CTkButton(self.parameters, border_width=1, border_color="black", text_color="black", text = "Instructions", command = self.get_instructions)
        self.doc_param.grid(row=6, column=0, padx=20, pady=(10, 5))
        
        self.appearance = ctk.CTkLabel(self.parameters, text="Appearance:", font=ctk.CTkFont(size=10, weight="bold"))
        self.appearance.grid(row=7, column=0, padx=20, pady=(5, 5))
        
        self.appearance_mode_optionemenu = ctk.CTkOptionMenu(self.parameters, values=["Dark", "Light"], command=self.change_appearance_mode_event)
        self.appearance_mode_optionemenu.grid(row=7, column=0, padx=20, pady=(5, 10))
        
        
        
        self.separator1 = tk.Canvas(self, width=2, bg='black')
        self.separator1.grid(row=0, column=1, sticky='ns')
        self.separator1.create_line(1, 0, 1, self.separator1.winfo_height())
        
        
        
        self.pipeline = ctk.CTkFrame(self, width=100, corner_radius=2)
        self.pipeline.grid(row=0, column=2, rowspan=4, sticky="nsew")
        self.pipeline.grid_columnconfigure(0, weight=1)
        self.pipeline.grid_rowconfigure(4, weight=1)
        
        self.logo_pipeline = ctk.CTkLabel(self.pipeline, text="Pipeline", font=ctk.CTkFont(size=20, weight="bold"))
        self.logo_pipeline.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        self.advanced_settings = ctk.CTkButton(self.pipeline, border_width=1, border_color="black", text_color="black", text = "Advanced Settings", command = self.advanced_setting)
        self.advanced_settings.grid(row=1, column=0, padx=20, pady=(10, 40))
        
        self.first_script = ctk.CTkButton(self.pipeline, border_width=1, border_color="black", text_color="black", text = "Run Length Detection", command = self.run_first_script)
        self.first_script.grid(row=2, column=0, padx=20, pady=(10, 10))
        
        self.second_script = ctk.CTkButton(self.pipeline, border_width=1, border_color="black", text_color="black", text = "Run SNP Detection", command = self.run_second_script)
        self.second_script.grid(row=3, column=0, padx=20, pady=(10, 20))
        
        
        self.textbox = ctk.CTkTextbox(self.pipeline, width=250)
        self.textbox.grid(row=4, column=0, padx=(5, 5), pady=(5, 5), sticky="nsew")
        self.textbox.insert("0.0", "Pipeline:" + "\n\n" * 2)
        #self.textbox.insert("0.0", "Pipeline:" + "\n\n" * 2)
        
        
        self.report = ctk.CTkButton(self.pipeline, border_width=1, border_color="black", text_color="black", text = "Get Report", command = self.get_report)
        self.report.grid(row=5, column=0, padx=20, pady=(20, 5))
        
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
        
        

    def parse_params(self):
        params = os.path.normpath(self.workspace + '/parameters.txt')
        if os.path.exists(params):
            with open(params, 'r') as paramfile:
                for line in paramfile:
                    line = line.rstrip("\r\n")
                    if (not line) or line.startswith("#"):
                        continue
                    #print(line)
                    self.set_param(line)
        else:
            print("\33[31;47m" + "No parameter file found with the name 'parameters.txt'!")
        
    def set_param(self, line):
        elements = line.split("=")
        key = elements[0].rstrip("\r\n").strip()
        value = elements[1].rstrip("\r\n").strip(".").strip() #.strip(".")
        if (value and not(value == ".") and not(value == "")):
            self.params[key] = value
            
    
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
        
        
        # entries_dirs = []
        # entries_files = []
        # entries_calc = []
        # entries_special = []
        # entries_total = []
        entries = []
        
        param_list_dirs = ["OutputFolder", "Bin", "RawData"]
        number1 = 0
        for i, param in enumerate(param_list_dirs):
            ctk.CTkLabel(dir_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(dir_params, corner_radius=5) #, textvariable=text
            if not((param == "OutputFolder")):
                browse_button = tk.Button(dir_params, text="Browse", command=lambda entry = param_entry: self.browse_dir(entry))
                browse_button.grid(column=2, row=i, padx=5, pady=5)
            elif (param == "OutputFolder"): #browse_dir_outputfolder
                browse_button = tk.Button(dir_params, text="Browse", command=lambda entry = param_entry: self.browse_dir_outputfolder(entry))
                browse_button.grid(column=2, row=i, padx=5, pady=5)
            param_entry.insert(tk.END, os.path.normpath(self.params[param]))
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            entries.append((param, param_entry))
            #param_entry.insert(0, default_values[i])
            
        
        
        param_list_files = ["PrimerFile", "SampleFile", "RExecutable", "AlleleList"]
        number2 = 0
        for i, param in enumerate(param_list_files):
            ctk.CTkLabel(file_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(file_params, corner_radius=5) #, textvariable=text
            browse_button = tk.Button(file_params, text="Browse", command=lambda entry=param_entry: self.browse_file_folder(entry))
            browse_button.grid(column=2, row=i, padx=5, pady=5)
            param_entry.insert(tk.END, os.path.normpath(self.params[param]))
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            entries.append((param, param_entry))
            
        
        
        #entries_total.append(entries_files)
        
        
        
        param_list_calc = ["MaxMismatch", "MinCount", "MinLength", "ConsensusThreshold", "LengthWindow"]
        number3 = 0
        for i, param in enumerate(param_list_calc):
            ctk.CTkLabel(calc_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(calc_params, corner_radius=5) #, textvariable=text
            param_entry.insert(tk.END, self.params[param])
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            entries.append((param, param_entry))
        
        param_list_special = ["Ploidy", "Operatingsystem", "Uniqueidentifier", "Indexcomboposition"]
        number4 = 0
        for i, param in enumerate(param_list_special):
            ctk.CTkLabel(special_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(special_params, corner_radius=5) #, textvariable=text
            param_entry.insert(tk.END, self.params[param])
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            entries.append((param, param_entry))
        
        write_params_button = tk.Button(newWindow, text="Write ParameterFile", command = lambda x = entries : self.write_params(x))
        write_params_button.grid(column=1, row=4, sticky="se", padx=15, pady=4)
        
        #, command = lambda x = entries_total : self.set_general_params(x)
        update_params_button = tk.Button(newWindow, text="Update Params", command = lambda x = entries : self.set_params(x))
        update_params_button.grid(column=0, row=4, sticky="sw", padx=15, pady=4)
    
    def set_params(self, entries):
        self.params = {entry[0] : entry[1].get() for entry in entries}       
        self.params["Reference"] = ""
        
        
    def write_params(self, entries):
        
        self.set_params(entries)
        
        with open("parameters.txt", 'w') as file:
            file.write("###Parameters for the Pipeline \n")
            file.write("\n")
            file.write("#Name of the Folder where results of the pipeline are stored: \n")
            file.write("OutputFolder = " + os.path.normpath(self.params["OutputFolder"]) + "\n")
            file.write("#Path to the Bin Folder where Trimmo, Usearch and R executables are stored: \n")
            file.write("Bin = " + os.path.normpath(self.params["Bin"]) + "\n")
            file.write("#Path to the Folder containing raw fastq/fastq.gz files: \n")
            file.write("RawData = " + os.path.normpath(self.params["RawData"]) + "\n")
            file.write("\n")
            file.write("#Path to the File containing Primer information (.txt): \n")
            file.write("PrimerFile = " + os.path.normpath(self.params["PrimerFile"]) + " \n")
            file.write('#Path to the File containing sample and indexcombo information (.csv): \n')
            file.write("SampleFile = " + os.path.normpath(self.params["SampleFile"]) + "\n")
            file.write('#Path to the Rexecutable (Only needed in windows) (.exe): \n')
            file.write("RExecutable = " + os.path.normpath(self.params["RExecutable"]) + "\n")
            file.write('#Path to the Allelelist (if user wants to do listbased Call) (.txt): \n')
            file.write("AlleleList = " + os.path.normpath(self.params["AlleleList"]) + " \n")
            file.write("\n")
            file.write('#Maximum number of mismatches between primer and amplicon sequences (integer; default=2) \n')
            file.write("MaxMismatch = " + str(self.params["MaxMismatch"]) + " \n")
            file.write('#Minimum number of read count for an allele to be considered (integer; default = 20) \n')
            file.write("MinCount = " + str(self.params["MinCount"]) + " \n")
            file.write('#Minimum amplicon length to be considered (integer; default=250) \n')
            file.write("MinLength= " + str(self.params["MinLength"]) + " \n")
            file.write('#Consensus threshold (float; default = 0.7) \n')
            file.write("ConsensusThreshold = " + str(self.params["ConsensusThreshold"]) + " \n")
            #self.lengthWindow = self.lengthWindow.split()
            liste = self.params["LengthWindow"].split(",")
            file.write("LengthWindow= " + str(liste[0].rstrip("\r\n")) + "," + str(liste[1].rstrip("\r\n")) + "\n")
            file.write("\n")
            file.write("#Ploidy level of the data (haploid/diploid; default = diploid): \n")
            file.write("Ploidy = " + str(self.params["Ploidy"]) + "\n")
            file.write("#Operating System of the current user (default = windows): \n")
            file.write("Operatingsystem = " + str(self.params["Operatingsystem"]) + " \n")
            file.write("#Unique Identifier for each Sequence: \n")
            file.write("Uniqueidentifier = " + str(self.params["Uniqueidentifier"])+ " \n")
            file.write("#Position of the Indexcombination in the raw Fastqfilename (default = 1): \n")
            file.write("Indexcomboposition = " + str(self.params["Indexcomboposition"]))
            #file.write("\n")
            #file.write("Reference = " + str(os.path.normpath(self.params["Reference"])))
            
        
        
        if self.os_param == "windows":
            os.startfile("parameters.txt")
        elif self.os_param == "linux":
            subprocess.run(["xdg-open", "parameters.txt"], check=True)
    
    def get_instructions(self):
        print("Get Instructions:")
        
        #["ColumnNrSample", "ColumnNrIndexComb", "ColumnNrOrganism", "ColumnNrProject", "ColumnNrCountry", "ColumnNrLocality"]
        
        with open("Instructions.txt", 'w') as file:
            file.write("###Instructions for handling the SSR-GBAS GUI: \n\n\n\n")
            file.write("#The GUI is separated into 3 columns: We prepare our workspace through the left column. \n\n")
            file.write("The 'Preparation' column handles the import of input parameter which are necessary in order to start the pipeline.\n")
            file.write("In the 'Pipeline' column the user can start the SSR-GBAS pipeline or use the advanced section to start individual parts of the pipeline.\n")
            file.write("The 'Database' column handles the storage of output data into a local database for further analysis if the user chooses to.\n\n\n\n")
            
            file.write("Preparation: \n \n")
            file.write("By clicking the 'Pipeline Parameter' button the user can set up all input parameters which are are necessary for running the pipeline.\n")
            file.write("These parameters are: \n")
            file.write("-) A bin folder which contains the executables Trimmomatic (along with a folder for adapters), Usearch and the R script 'Rscript_Markerlength_STUTTER_Color_BaryzentricMinsize20_notVerbose.R'. The bin folder also needs to be in the same location as the GUI executable in order for the GUI to recognize it. \n")
            file.write("-) An output folder where results will be stored \n")
            file.write("-) A folder with fastq files in GZ format which contain the microsatellite sequences \n")
            file.write("-) A folder with fastq files in GZ format \n")
            file.write("-) A txt file with primers (in the form of 3 columns: primername, forwardsequence, reversesequence) with delimiter being Comma (',') \n")
            file.write("-) A samplesheet which as a csv (in the form of 2 columns: indexcombination, sampleID) with delimiter being either ',' or ';' or '\t' \n")
            file.write("-) An allelelist which was created from a previous pipeline run and will be used as the base for the current run. \n")
            file.write("-) Furthermore you need to download R if you dont have it yet and if you work on windows you need to set the path to the executable Rscript.exe \n\n")
            
            file.write("After everything is set up use the Workspace Status Button to check if all input paths are read correctly and if the GUI can properly access the fastqfiles! (The GUI looks for fastqfiles with names containing the inxecombinations of the samplesheet) \n \n")
            
            file.write("The 'Genalex Parameter' button must be clicked and set on a certain parameter (Organism, Country, Locality) if the user wishes to later on have his output in genalex format. \n")
            file.write("For that the user must click the 'Reference Parameters' button to to import a reference file (in csv format and same delimiter as samplesheet) and then choose which columns represent which parameter. Furthermore the user can also create a samplesheet from the reference file if he does not have one yet. \n\n\n\n")
            
            
            file.write("Pipeline: \n \n")
            file.write("Once the workspace is set up with all necessary input parameters the user can start the first part of the pipeline by clicking the 'Run Length Detection' button. \n")
            file.write("This begins processing the fastqfiles by first trimming the sequences of their adapter sequeces and low quality regions using Trimmomatic, merging them using Usearch, cutting off sequences where depending on primer forward and reverse mismatches and filtering sequences based on their lengths. Lastly the Rscript chooses the most likely lengths (depending on count frequence) as the genotype for that microsatellite. \n")
            file.write("The final output results are in the 'Markerplot' folder which contains all the different sequence lengths for each marker-sample combination (counting frequency vs length plots) and a matrix (columns representing marker, rows representing samples) and each cell representing the genotype for that marker-sample combination (csv format) . \n\n")
            file.write("After the process is done the user can select the 'SNP Detection' to start the second part of the pipeline which involves finding all SNPs for each marker and for all marker-sample combinations. \n")
            file.write("This involves the extraction of the genotypes of the first output matrix, building consensus sequences per length and choosing the most likely variants as the SNPs for that particular marker. \n")
            file.write("The final output is stored in the AlleleCall folder and contains 2 txt files, the first (with the name allele_list.txt) showing all SNP sequences per marker and the second (matrix.txt) being a matrix (similar as the first outpt matrix) but now the cells representing an index, this index corresponds to a sequence in the allele_list txt filer \n\n")
            file.write("The 'Advanced Settings' button allows to use individual parts of the pipeline instead of running the entire process again. This is in case the user stumbles across storage problems. \n\n")
            file.write("The 'Get Report' button allows of gaining more information of the previously run pipeline process and 'Import as Genalex' allows the user to select an AlleleCall folder to write the output matrix into genalex format given a certain choosen parameter in the 'Preparation' column.  \n\n\n\n")
            
            
            file.write("Database: \n \n")
            file.write("The user has the option to store his output data into a personal database. Thi uses the SQlite Module from python to store data into a .db file using SQL. (Though in the current version the storing of data is very inefficient and should therefore for now be used on a smaller dataset.  \n")
            file.write("Clicking the 'Add Data' button allows the user to choose the an AlleleCall folder where his output is stored. There the SNPs are stored into a local database. If a reference file has been set then inforation regarding organism, country and locality will be stored as well, otherwise those fields will be shown as 'None'. \n")
            file.write("The 'Extract Subset' button allows the user to filter his output sequences according to the parameters Organism, Sample and Loci. This allows for a more customized output (both as a fasta or genalex format) \n")
            
            file.write("The final output results are in the 'Markerplot' folder which contains all the different sequence lengths for each marker-sample combination (counting frequency vs length plots) and a matrix (columns representing marker, rows representing samples) and each cell representing the genotype for that marker-sample combination (csv format) . \n\n")
            file.write("After the process is done the user can select the 'SNP Detection' to start the second part of the pipeline which involves finding all SNPs for each marker and for all marker-sample combinations. \n")
            file.write("This involves the extraction of the genotypes of the first output matrix, building consensus sequences per length and choosing the most likely variants as the SNPs for that particular marker. \n")
            file.write("The final output is stored in the AlleleCall folder and contains 2 txt files, one showing each genetic SNP variant per marker and the second being a matrix (similar as the first outpt matrix) but now the cells representing an index, this index corresponds to a sequence in the allele_list txt filer \n\n\n\n")
            file.write("The 'Advanced Settings' button allows to use individual parts of the pipeline instead of running the entire process again. This is in case the user stumbles across storage problems. \n\n")
            file.write("The 'Get Report' button allows of gaining more information of the previously run pipeline process and 'Import as Genalex' allows the user to select an AlleleCall folder to write the output matrix into genalex format given a certain choosen parameter in the 'Preparation' column.  \n\n")
            
            file.write("For more questions or problems please write a mail to sonnenbe@groupwise.boku.ac.at \n")
            
            
            
        
        if (self.os_param == "windows"):
            os.startfile("Instructions.txt")
        elif (self.os_param == "linux"):
            subprocess.call(["xdg-open", "Instructions.txt"])
    
    def get_error_window(self, title, message):
        newWindow = tk.Toplevel(self)
        newWindow.title(title)
        newWindow.grid_rowconfigure(0, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        logo_samplesheet = ctk.CTkLabel(newWindow, text=message, font=ctk.CTkFont(size=20, weight="bold"))
        logo_samplesheet.grid(row=0, column=0, padx=20, pady=(20, 20))
        
    def get_sample_number(self):
        print("Getting sample Number")
        number = 0
        samples = []
        
        if (not(os.path.exists(os.path.normpath(self.sample_sheet))) or self.sample_sheet == "."):
            return number, []
        
        with open(self.sample_sheet) as file:  #codecs.open(self.sample_sheet, encoding='utf8') with open(self.sample_sheet) as file
            for line in file:
                number += 1
                samples.append((re.split(',|\t|;', line)[0]))
        
        #with open(self.sample_sheet, "r") as samplesheet:
            #csv_reader = csv.reader(samplesheet, delimiter=';, ,')
            #for line in csv_reader:
                #number += 1
                #samples.append(line[0])
        
        return number, list(set(samples))
    

    def get_samples(self):
        print("Getting sample Number")
        number = 0
        samples = []
        samples_duplicates = {}
        indexcombos = []
        indexcombos_duplicates = {}
        incomplete_lines = 0
        incomplete_samples = []
        
        
        
        if (self.dict_check["SampleFile"]):
            with open(self.params["SampleFile"]) as samplefile:
                for line in samplefile:
                    if (len(re.split(',|\t|;', line)) == 1):
                        print("Reached")
                        incomplete_lines += 1
                        incomplete_samples.append(re.split(',|\t|;', line)[0])
                        continue
                    line = line.rstrip('\r\n')
                    print(len(re.split(',|\t|;', line)))
                    if (len(re.split(',|\t|;', line)) > 1):
                        number += 1
                        sample = re.split(',|\t|;', line)[1]
                        indexcombo = re.split(',|\t|;', line)[0]
                        if (sample in samples):
                            #self.params["SampleDuplicates"]
                            #samples_duplicates.append(sample)
                            if (sample not in samples_duplicates):
                                samples_duplicates[sample] = 1
                            samples_duplicates[sample] += 1
                        samples.append(sample)
                        if (indexcombo in indexcombos):
                            #self.params["SampleDuplicates"]
                            #samples_duplicates.append(sample)
                            if (indexcombo not in indexcombos_duplicates):
                                indexcombos_duplicates[indexcombo] = 1
                            indexcombos_duplicates[indexcombo] += 1
                        indexcombos.append(indexcombo)
                    # elif (len(re.split(' ', line)) == 1):
                    #     incomplete_lines += 1
                    #     incomplete_samples.append(re.split(',|\t|;', line)[0])
                    else:
                        incomplete_lines += 1
                        
            
        return number, samples, samples_duplicates, indexcombos, indexcombos_duplicates, incomplete_lines, incomplete_samples
        
    def get_rawdata(self):
        number = 0
        list_files = []
        files_duplicates = {}
        filetype_check = True
        
        if (not(os.path.exists(os.path.normpath(self.params["RawData"]))) or os.path.basename(os.path.normpath(self.params["RawData"])) == "None" or os.path.basename(os.path.normpath(self.params["RawData"])) == "" or os.path.basename(os.path.normpath(self.params["RawData"])) == "."):
            return number, [], files_duplicates, filetype_check
        
        for root, dirs, files in os.walk(self.params["RawData"]):
            number += len(files)
            for file in files:
                indexcombo = file.split("_")[int(self.params["Indexcomboposition"])-1]
                if (file.split(".")[-2] != "fastq" and file.split(".")[-1] != "gz"):
                    #print(file.split(".")[-2])
                    filetype_check = False
                if (indexcombo not in files_duplicates):
                    files_duplicates[indexcombo] = 1
                else:
                    files_duplicates[indexcombo] += 1
                #print(files_duplicates[indexcombo])
                list_files.append(file)
                
        
        return number, list_files, files_duplicates, filetype_check
    
    def get_compatibility(self, indexcombos, fastqfiles):
        number = 0
        for indexcombo in indexcombos:
            for fastqfile in fastqfiles:
                if (indexcombo in fastqfile):
                    number += 1
        return number    

    def check_status(self):
        new_window = tk.Toplevel(self)
        new_window.title("Workspace Status")
        new_window.geometry(f"{1100}x{700}")
        long_text = f"Input files and folders: \n\n\n"
        start_check = True
        #dict_check = {"OutputFolder" : False, "Bin" : False, "RawData" : False, "PrimerFile" : False, "SampleFile" : False, "Reference" : False, "RExecutable" : False}
        for param, val in self.dict_check.items():
            if (os.path.exists(os.path.normpath(self.params[param])) and os.path.basename(os.path.normpath(self.params[param])) != "None" and os.path.basename(os.path.normpath(self.params[param])) != "" and os.path.basename(os.path.normpath(self.params[param])) != "."):
                self.dict_check[param] = True
                long_text += param + " : " + os.path.normpath(self.params[param]) + "\n\n"
            else:
                self.dict_check[param] = False
                if (param != "Reference"):
                    start_check = False
                long_text += param + " : not set or not valid! \n\n"
        
        long_text += "\n\n\n"
        nrsamples, sampleIDs, sampleIDs_duplicates, indexcombos, indexcombos_duplicates, incomplete_lines, incomplete_samples = self.get_samples()
        number_rawdata, list_rawfiles, rawfiles_duplicates, filetype_check = self.get_rawdata()
        if (start_check):
            long_text += os.path.basename(os.path.normpath(self.params["SampleFile"])) + " " + str(nrsamples) + " IndexcomboNames and " + str(len(indexcombos_duplicates)) + " duplicate IndexcomboNames. \n"
            long_text += os.path.basename(os.path.normpath(self.params["SampleFile"])) + " " + str(nrsamples) + " SampleIDs and " + str(len(sampleIDs_duplicates)) + " duplicate SampleIDs. \n\n\n"
            if (len(indexcombos_duplicates) > 0):
                start_check = False
                long_text += "Duplicate IndexcomboNames: \n\n"
                #long_text += "Pipeline cannot start! \n\n"
                for duplicate, count  in indexcombos_duplicates.items():
                    long_text += duplicate + " appears " + str(count) + " times. \n"
                long_text += "\n"
            if (len(sampleIDs_duplicates) > 0):
                start_check = False
                long_text += "Duplicate SampleIDs: \n\n"
                #long_text += "Pipeline cannot start! \n\n"
                for duplicate, count  in sampleIDs_duplicates.items():
                    long_text += duplicate + " appears " + str(count) + " times. \n"
                long_text += "\n"
            if (incomplete_lines > 0):
                start_check = False
                long_text += "Found " + str(incomplete_lines) + " incomplete lines. \n"
                for indexcombo in incomplete_samples:
                    long_text += "Line with indexcombo " + str(indexcombo) + " has no SampleID. \n"
            # if (len(indexcombos_duplicates) > 0 or len(sampleIDs_duplicates) > 0 or incomplete_lines > 0):
            #     long_text += "Pipeline cannot start! \n\n"
            # else:
            #     long_text += "\n\n"
    
            long_text += "\n\n"
            long_text += os.path.basename(os.path.normpath(self.params["RawData"])) + " " + str(number_rawdata) + " zipped files "
            if (filetype_check):
                long_text += "with filetype fastq.gz. \n\n"
            else:
                long_text += ". \n"
                long_text += "Not all are of type fastq!! Check! \n"
                start_check = False
            for filename, count in rawfiles_duplicates.items():
                if (rawfiles_duplicates[filename] > 2):
                    long_text += "File with indexcombo " + str(filename) + " appears more than 2 times! Check! \n"
                    start_check = False
            
            if (not(number_rawdata % 2 == 0)):
                long_text += "CAREFUL! Non even number of zipped fastqfiles which means not the same number of forward and reverse files! \n\n"
            
        if (start_check):
            nr_sample_fastq_compatible = self.get_compatibility(indexcombos, list_rawfiles)
            if (nr_sample_fastq_compatible > 0):
                long_text += str(nr_sample_fastq_compatible) + " zipped fastqfiles will be processed in the pipeline! \n"
                if (nr_sample_fastq_compatible*2 < number_rawdata):
                    long_text += "Only " + str(nr_sample_fastq_compatible) + " zipped fastqfiles will be processed in the pipeline! Eventually check your samplesheet! \n"
                else:
                    long_text += str(nr_sample_fastq_compatible) + " zipped fastqfiles will be processed in the pipeline! \n"
            else:
                long_text += "No zipped fastqfiles will be processed. There is no connection between indexcombos of samplesheet and in filesnames of rawdata folder! \n\n"
            
        if (start_check):
            long_text += "Pipeline is save to start! \n\n"
        else:
            long_text += "Pipeline should not be started! \n\n"
        

            
            
        # nrfastq, fastqfiles = self.get_rawdata_number()
        # nr_sample_fastq_compatible = self.get_compatibility(samples, fastqfiles)
        
        text_widget = tk.Text(new_window, wrap=tk.WORD, height=20, width=50, font=("Helvetica", 16))
        scrollbar = tk.Scrollbar(new_window, command=text_widget.yview)
        text_widget.config(yscrollcommand=scrollbar.set)
        
        text_widget.insert(tk.END, long_text)
        text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=20, pady=20)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        
                
            
    
    def check_status2(self):
        #print("All Input Files and Folders:")
        dict_check = {"OutputFolder" : True, "Bin" : True, "RawData" : True, "Primerfile" : True, "SampleSheet" : True, "Reference" : True, "RScript" : True, "bin" : True}
        
        if (not(os.path.exists(os.path.normpath(self.workspace_location)))):
            dict_check["workspace"] = False
        if (not(os.path.exists(os.path.normpath(self.rawdata_dir))) or self.rawdata_dir == "" or self.rawdata_dir == "."):
            dict_check["rawdata"] = False
        if (not(os.path.exists(os.path.normpath(self.primerfile))) or self.primerfile == "" or self.primerfile == "."):
            dict_check["Primerfile"] = False
        if (not(os.path.exists(os.path.normpath(self.sample_sheet))) or self.sample_sheet == "" or self.sample_sheet == "."):
            dict_check["SampleSheet"] = False
        if (not(os.path.exists(os.path.normpath(self.reference_file))) and not(self.reference_file.split(".")[1] == "csv")):
            dict_check["Reference"] = False
        if (not(os.path.exists(os.path.normpath(self.Rexe))) and not(self.Rexe.split(".")[0] == "Rscript")):
            dict_check["RScript"] = False
        if (not(os.path.exists(os.path.normpath(self.bin_dir)))):
            dict_check["bin"] = False
        
        nrsamples, samples = self.get_sample_number()
        nrfastq, fastqfiles = self.get_rawdata_number()
        nr_sample_fastq_compatible = self.get_compatibility(samples, fastqfiles)
        #nrsamples = 9
        #nrfastq = int([file for file in os.listdir(os.path.normpath(self.workspace_location + '/' + self.rawdata_dir)) if ()][0])
        
        print()
        #print(samples)
        print()
        #print(fastqfiles)
        
        new_window = tk.Toplevel(self)
        new_window.title("Workspace Status")
        long_text = f"All input files and folders: \n\n\n"
        if (dict_check["workspace"] == True):
            long_text += "Workspace Location : " + str(os.path.normpath(self.workspace_location)) + "\n\n"
        else:
            long_text += "Workspace Location : Not a valid path! \n\n"
        if (dict_check["rawdata"]):
            long_text += "Rawdata : " + str(os.path.basename(os.path.normpath(self.rawdata_dir))) + "\n\n"
        else:
            long_text += "Rawdata : Not valid or missing! \n\n"
        if (dict_check["Primerfile"]):
            long_text += "Primerfile : " + str(os.path.basename(os.path.normpath(self.primerfile))) + "\n\n"
        else:
            long_text += "Primerfile : Not valid or missing! \n\n"
        if (dict_check["SampleSheet"]):
            long_text += "SampleSheet : " + str(os.path.basename(os.path.normpath(self.sample_sheet))) + "\n\n"
        else:
            long_text += "SampleSheet : Not valid or missing! \n\n"
        if (dict_check["Reference"]):
            long_text += "ReferenceFilePath : " + str(os.path.basename(os.path.normpath(self.reference_file))) + "\n\n"
        else:
            long_text += "ReferenceFilePath : Not valid or missing! \n\n"
        if (dict_check["RScript"]):
            long_text += "RScriptExePath : " + str(os.path.normpath(self.Rexe)) + "\n\n"
        else:
            long_text += "RScriptExePath : Not valid or missing! \n\n"
        if (dict_check["bin"]):
            long_text += "bin : exists! \n\n\n\n"
        else:
            long_text += "bin : Not valid or missing! \n\n\n\n"
        
        
        long_text += "Workspace Status: \n\n"
        long_text += "Samplesheet File provides " + str(nrsamples) + " SampleIds. \n"
        long_text += "Rawdata Folder provides " + str(nrfastq) + " zipped FastQfiles. \n"
        if (not(nrfastq % 2 == 0)):
            long_text += "CAREFUL! Non even number of fastqfiles which means not the same number of forward and reverse files! \n\n"
        else:
            long_text += "\n"
        
        long_text += str(nr_sample_fastq_compatible) + " fastqfiles of the raw datadir will be processed in the pipeline! \n"
        if (nr_sample_fastq_compatible == 0):
            long_text += "This means there is no fit between names of the samplesheet and the names of the raw datadir. Wrong samplesheet or raw datadir! \n"
        else:
            long_text += "\n"
        
        text_widget = tk.Text(new_window, wrap=tk.WORD, height=20, width=50, font=("Helvetica", 16))
        scrollbar = tk.Scrollbar(new_window, command=text_widget.yview)
        text_widget.config(yscrollcommand=scrollbar.set)
        
        text_widget.insert(tk.END, long_text)
        text_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=20, pady=20)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    
    def get_compatibility2(self, samples, fastqfiles):
        number = 0
        for sample in samples:
            for fastqfile in fastqfiles:
                if (sample in fastqfile):
                    number += 1
        return number
    
    def get_rawdata_number(self):
        #self.rawdata_dir
        print("Rawdata file Number")
        number = 0
        list_files = []
        
        if (not(os.path.exists(os.path.normpath(self.rawdata_dir)))):
            return number, []
        
        for root, dirs, files in os.walk(self.rawdata_dir):
            number += len(files)
            for file in files:
                list_files.append(file)
                
            
        return number, list_files
    
    def get_QC_files(self):
        #self.rawdata_dir
        #print("Rawdata file Number")
        number = 0
        list_files = []
        
        if (not(os.path.exists(os.path.normpath(self.output_folder_path + "/QC")))):
            return number, []
        
        for root, dirs, files in os.walk(self.rawdata_dir):
            number += len(files)
            for file in files:
                list_files.append(file)
                
            
        return number, list_files
            
    
    def set_reference_import_params(self):
        newWindow = tk.Toplevel(self)
        newWindow.title("Reference Import Parameters")
        newWindow.grid_rowconfigure(0, weight=1)
        newWindow.grid_rowconfigure(1, weight=1)
        newWindow.grid_rowconfigure(2, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        newWindow.grid_columnconfigure(1, weight=1)
        
        print("DA:")
        print(self.reference_file)
        
        info_ref = ctk.CTkFrame(newWindow, width=100, corner_radius=2)
        info_ref.grid(row=0, column=0, rowspan=1, sticky="nsew")
        info_ref.grid_columnconfigure(0, weight=1)
        info_ref.grid_rowconfigure(0, weight=1)
        
        info_columns_logo = ctk.CTkLabel(info_ref, text="Set Path for ReferenceFile:", font=ctk.CTkFont(size=20, weight="bold"))
        info_columns_logo.grid(row=0, column=0, sticky="nsew", padx=20, pady=(20, 20))
        
        ctk.CTkLabel(info_ref, text="Reference:", width=20).grid(column=0, row=1, padx=5, pady=5)
        param_entry_ref = ctk.CTkEntry(info_ref, corner_radius=5) #, textvariable=text
        browse_button = tk.Button(info_ref, text="Browse", command=lambda entry=param_entry_ref: self.browse_file_folder(entry))
        browse_button.grid(column=2, row=1, padx=5, pady=5)
        param_entry_ref.insert(tk.END, self.reference_file)
        param_entry_ref.grid(column=1, row=1, padx=5, pady=5)
        #param_entry.insert(0, default_values[i])
        #self.param_list_files[3] = os.path.normpath(param_entry.get())
        #self.reference_file = os.path.normpath(param_entry.get())
        
        info_columns = ctk.CTkFrame(newWindow, width=100, corner_radius=2)
        info_columns.grid(row=1, column=0, rowspan=1, sticky="nsew")
        info_columns.grid_columnconfigure(0, weight=1)
        info_columns.grid_rowconfigure(0, weight=1)
        
        info_columns_logo = ctk.CTkLabel(info_columns, text="Set ColumnIndex for Params:", font=ctk.CTkFont(size=20, weight="bold"))
        info_columns_logo.grid(row=0, column=0, sticky="nsew", padx=20, pady=(20, 20))
        
        param_list_columns = ["ColumnNrSample", "ColumnNrIndexComb", "ColumnNrOrganism", "ColumnNrProject", "ColumnNrCountry", "ColumnNrLocality"]
        number_column = 0
        column_entries = []
        for i, param in enumerate(param_list_columns):
            ctk.CTkLabel(info_columns, text=param + ":", width=20).grid(column=0, row=i+1, padx=5, pady=5)
            #text = tk.StringVar()
            param_entry = ctk.CTkEntry(info_columns, corner_radius=5) #, textvariable=text
            param_entry.insert(tk.END, str(self.columnr_params[i]))
            param_entry.grid(column=1, row=i+1, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            column_entries.append(param_entry)
            number_column += 1
        
        
        info_rows = ctk.CTkFrame(newWindow, width=100, corner_radius=2)
        info_rows.grid(row=2, column=0, rowspan=1, sticky="nsew")
        info_rows.grid_columnconfigure(0, weight=1)
        info_rows.grid_rowconfigure(0, weight=1)
        
        info_rows_logo = ctk.CTkLabel(info_rows, text="Set RowIndex for Params:", font=ctk.CTkFont(size=20, weight="bold"))
        info_rows_logo.grid(row=0, column=0, sticky="nsew", padx=20, pady=(20, 20))
        
        param_list_rows = ["RowNrBegin", "RowNrEnd"]
        number_rows = 0
        row_entries = []
        for i, param in enumerate(param_list_rows):
            ctk.CTkLabel(info_rows, text=param + ":", width=20).grid(column=0, row=i+1, padx=5, pady=5)
            #text = tk.StringVar()
            param_entry = ctk.CTkEntry(info_rows, corner_radius=5) #, textvariable=text
            param_entry.insert(tk.END, str(self.rownr_params[i]))
            param_entry.grid(column=1, row=i+1, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            row_entries.append(param_entry)
            number_rows += 1
        
        
        
        write_samplesheet_button = ttk.Button(newWindow, text="Write Samplesheet", command = lambda x = column_entries, y = row_entries, ref = param_entry_ref : self.set_samplesheetname_name(x, y, ref))
        write_samplesheet_button.grid(column=0, row=4, sticky="se", padx=15, pady=4)
        
        write_referenceparameters_button = ttk.Button(newWindow, text="Update and Write Reference Params", command = lambda x = column_entries, y = row_entries, ref = param_entry_ref : self.write_reference_parameterfile(x, y, ref))
        write_referenceparameters_button.grid(column=0, row=4, sticky="sw", padx=15, pady=4)
    
    def check_if_validparams(self, x):
        try:
            return [int(value.get()) for value in x]
        except:
            return None
            
    def set_reference_params(self, x, y, ref):
        
        if (self.check_if_validparams(x) == None):
            newWindow = tk.Toplevel(self)
            newWindow.title("!!")
            newWindow.grid_rowconfigure(0, weight=1)
            newWindow.grid_columnconfigure(0, weight=1)
            logo_samplesheet = ctk.CTkLabel(newWindow, text="A Non Valid Value has been entered for one of the Columns! Cannot Proceed! \n Please make sure that ALL ColumnNumbers are above Zero!", font=ctk.CTkFont(size=20, weight="bold"))
            logo_samplesheet.grid(row=0, column=0, padx=20, pady=(20, 20))
            raise Exception("One of the ColumnNr is not valid!")
        elif (self.check_if_validparams(y) == None):
            newWindow = tk.Toplevel(self)
            newWindow.title("!!")
            newWindow.grid_rowconfigure(0, weight=1)
            newWindow.grid_columnconfigure(0, weight=1)
            logo_samplesheet = ctk.CTkLabel(newWindow, text="A Non Valid Value has been entered for one of the Columns! Cannot Proceed! \n Please make sure that ALL RowNumbers are above Zero!", font=ctk.CTkFont(size=20, weight="bold"))
            logo_samplesheet.grid(row=0, column=0, padx=20, pady=(20, 20))
            raise Exception("One of the RowNr is not valid!")
        
        if (ref.get() == "" or ref.get() == "." or not(os.path.exists(os.path.normpath(ref.get())))):
            self.get_error_window("ReferenceFilePath Error", "The inputted Reference File Path is not valid!")
        else:
            self.reference_file = os.path.normpath(ref.get())
        
        self.columnr_params = [int(param.get()) for param in x]
        self.rownr_params = [int(param.get()) for param in y]
        #self.rawfastq_position_after = int(z.get())
        
        
        if (0 in self.columnr_params or 0 in self.rownr_params):
            self.get_error_window("!!", "We want NON-Zero Numbers! \n Please make sure that Numbers are ABOVE Zero!")
            raise Exception("ONLY NON-ZERO NUMBERS!")
        
        if (any(num < 0 for num in self.columnr_params) or any(num < 0 for num in self.rownr_params)):
            self.get_error_window("!!", "We want POSITIVE INTEGER Numbers! \n Please make sure that Numbers are ABOVE Zero!")
            raise Exception("ONLY Non-Zero POSITIVE NUMBERS!")
            
            
    
    def write_reference_parameterfile(self, x, y, ref):
        
        print("Write reference:")
        print(self.reference_file)
        self.set_reference_params(x, y, ref)
        
        #["ColumnNrSample", "ColumnNrIndexComb", "ColumnNrOrganism", "ColumnNrProject", "ColumnNrCountry", "ColumnNrLocality"]
        
        with open("Referenceparameters.txt", 'w') as file:
            file.write("###Column and Row-Parameters for the Importing from the Reference File: \n\n")
            file.write("#ColumnNumber for SAMPLE-ID: \n")
            file.write("ColumnNrSample = " + str(self.columnr_params[0]) + "\n")
            file.write("#ColumnNumber for INDEX-COMBINATION: \n")
            file.write("ColumnNrIndexCombo = " + str(self.columnr_params[1]) + "\n")
            file.write("#ColumnNumber for ORGANISM: \n")
            file.write("ColumnNrOrganism = " + str(self.columnr_params[2]) + "\n")
            file.write("#ColumnNumber for PROJECT: \n")
            file.write("ColumnNrProject = " + str(self.columnr_params[3]) + "\n")
            file.write("#ColumnNumber for Country: \n")
            file.write("ColumnNrCountry = " + str(self.columnr_params[4]) + "\n")
            file.write("#ColumnNumber for Locality: \n")
            file.write("ColumnNrLocality = " + str(self.columnr_params[5]) + "\n")
            file.write("\n")
            file.write("#RowNumber Beginning Index: \n")
            file.write("RowNrBegin = " + str(self.rownr_params[0]) + "\n")
            file.write("#RowNumber End Index: \n")
            file.write("RowNrEnd = " + str(self.rownr_params[1]) + "\n")
            file.write("\n\n")
            file.write("#Absolute Path to the Reference File (needs to be .csv): \n")
            file.write("ReferencefilePath = " + os.path.normpath(self.reference_file) + "\n")
            file.write("\n")
            
        
        if (self.os_param == "windows"):
            os.startfile("Referenceparameters.txt")
        elif (self.os_param == "linux"):
            subprocess.call(["xdg-open", "Referenceparameters.txt"])
            
    
    def parse_reference_parameterfile(self):
        column_param_list = []
        raw_param_list = []
        
        paramter_file_path = os.path.normpath(os.getcwd() + '/Referenceparameters.txt')
        if (os.path.exists(paramter_file_path)):
            with open(paramter_file_path) as parameters:
                for line in parameters:
                    if (not line) or line.startswith("#"):
                        continue
                    if (line.startswith("ColumnNrSample")):
                        ColumnNrSample = line.split("=")[1].rstrip("\r\n").strip()
                        if (not ColumnNrSample):
                            column_param_list.append(0)
                        else:
                            column_param_list.append(int(ColumnNrSample))
                        continue
                    if (line.startswith("ColumnNrIndexCombo")):
                        ColumnNrIndexCombo = line.split("=")[1].rstrip("\r\n").strip()
                        if (not ColumnNrIndexCombo):
                            column_param_list.append(0)
                        else:
                            column_param_list.append(int(ColumnNrIndexCombo))
                        continue
                    if (line.startswith("ColumnNrOrganism")):
                        ColumnNrOrganism = line.split("=")[1].rstrip("\r\n").strip()
                        if (not ColumnNrOrganism):
                            column_param_list.append(0)
                        else:
                            column_param_list.append(int(ColumnNrOrganism))
                        continue
                    if (line.startswith("ColumnNrProject")):
                        ColumnNrProject = line.split("=")[1].rstrip("\r\n").strip()
                        if (not ColumnNrProject):
                            column_param_list.append(0)
                        else:
                            column_param_list.append(int(ColumnNrProject))
                        continue
                    if (line.startswith("ColumnNrCountry")):
                        ColumnNrCountry = line.split("=")[1].rstrip("\r\n").strip()
                        if (not ColumnNrCountry):
                            column_param_list.append(0)
                        else:
                            column_param_list.append(int(ColumnNrCountry))
                        continue
                    if (line.startswith("ColumnNrLocality")):
                        ColumnNrLocality = line.split("=")[1].rstrip("\r\n").strip()
                        if (not ColumnNrLocality):
                            column_param_list.append(0)
                        else:
                            column_param_list.append(int(ColumnNrLocality))
                        continue
                    
                    if (line.startswith("RowNrBegin")):
                        RowNrBegin = line.split("=")[1].rstrip("\r\n").strip()
                        if (not RowNrBegin):
                            raw_param_list.append(0)
                        else:
                            raw_param_list.append(int(RowNrBegin))
                        continue
                    if (line.startswith("RowNrEnd")):
                        RowNrEnd = line.split("=")[1].rstrip("\r\n").strip()
                        if (not RowNrEnd):
                            raw_param_list.append(0)
                        else:
                            raw_param_list.append(int(RowNrEnd))
                        continue
                    if (line.startswith("ReferencefilePath")):
                        ReferencefilePath = line.split("=")[1].rstrip("\r\n").strip()
                        if (not ReferencefilePath or ReferencefilePath == "." or not(os.path.exists(os.path.normpath(ReferencefilePath)))):
                            self.reference_file = ""
                        else:
                            self.reference_file = os.path.normpath(ReferencefilePath)
                        continue
                    print("Hier:")
                    print(self.reference_file)
                    print()
            return column_param_list, raw_param_list
        
        else:
            return [3, 11, 2, 1, 6, 7], [2, 1000]
        
    
    
    def set_samplesheetname_name(self, x, y, ref):
        self.set_reference_params(x, y, ref)
        
        newWindow = tk.Toplevel(self)
        newWindow.title("SamplesheetName")
        newWindow.grid_rowconfigure(0, weight=1)
        newWindow.grid_columnconfigure(0, weight=1)
        
        logo_samplesheet = ctk.CTkLabel(newWindow, text="Enter a Samplesheet Name", font=ctk.CTkFont(size=20, weight="bold"))
        logo_samplesheet.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        samplesheet_entry = ctk.CTkEntry(newWindow)
        samplesheet_entry.grid(row = 1, column = 0, padx=5, pady=(20, 20))
        
        enter = ctk.CTkButton(newWindow, border_width=1, border_color="black", text_color="black", text = "Enter", command = lambda: self.writesamplesheet(samplesheet_entry.get())) #command = lambda x = rows : self.download_data(x)
        enter.grid(row=2, column=0, padx=20, pady=(20, 10))
        
        samplesheet_entry.bind("<Return>", lambda event: self.writesamplesheet(samplesheet_entry.get()))
        samplesheet_entry.focus_set()

    
    def writesamplesheet(self, samplesheet):
        
        sample_names = []
        sequences = []
        samplesheet = samplesheet + ".csv"
        
        if (not(os.path.exists(self.reference_file))):
            newWindow = tk.Toplevel(self)
            newWindow.title("!!")
            newWindow.grid_rowconfigure(0, weight=1)
            newWindow.grid_columnconfigure(0, weight=1)
            logo_samplesheet = ctk.CTkLabel(newWindow, text="The path to the ReferenceFile must be set in General Parameters!", font=ctk.CTkFont(size=20, weight="bold"))
            logo_samplesheet.grid(row=0, column=0, padx=20, pady=(20, 20))
            raise Exception("No Reference FilePath set!")
            
        sample_columnr = self.columnr_params[0]
        sequence_columnr = self.columnr_params[1]
        rownr_begin = self.rownr_params[0]
        rownr_end = self.rownr_params[1]
        
        rows_new = []
        
        
        with open(self.reference_file, "r") as reference:
            for i, line in enumerate(reference, 1):
                if (i < rownr_begin):
                    continue
                if (re.split(',|\t|;', line)[sample_columnr-1] == ""):
                    break
                if (i > rownr_end):
                    break
                sample_names.append(re.split(',|\t|;', line)[sample_columnr-1])
                sequences.append(re.split(',|\t|;', line)[sequence_columnr-1])
                rows_new.append([re.split(',|\t|;', line)[sequence_columnr-1], re.split(',|\t|;', line)[sample_columnr-1]])
        
        
        
        #with open(self.reference_file, "r") as reference:
            #csv_reader = csv.reader(reference, delimiter=';')
            #for i, line in enumerate(csv_reader, 1):
                #if (i < rownr_begin):
                    #continue
                #if (line[sample_columnr-1] == ""):
                    #break
                #if (i > rownr_end):
                    #break
                #sample_names.append(line[sample_columnr-1])
                #sequences.append(line[sequence_columnr-1])
                #rows_new.append([line[sequence_columnr-1], line[sample_columnr-1]])
        
        with open(samplesheet, 'w', newline='') as new_samplesheet:
            writer = csv.writer(new_samplesheet, delimiter =';', quoting=csv.QUOTE_MINIMAL)
            for row in rows_new:
                row = [element.replace("\n", "") for element in row]
                print(row)
                writer.writerow(row)
        
        if self.os_param == "windows":
            os.startfile(samplesheet)
        elif self.os_param == "linux":
            subprocess.run(["xdg-open", samplesheet], check=True)
        
        self.sample_names = sample_names
        
    def closeref(self, ref_path):
        self.reference_file = os.path.normpath(ref_path)
        
    
    def get_analysis(self):
        print("Get the analysis file")
        
        analysis_path = os.path.normpath(self.output_dir + "/Analysis/analysis.txt")
            
        
        if (self.os_param == "windows"):
            os.startfile(analysis_path)
        elif (self.os_param == "linux"):
            subprocess.call(["xdg-open", analysis_path])
        
        #filepath = os.path.normpath("")
    
    def download_data(self, x):
       print("Downloading")
       filename = "Subset.fasta"
       
       os.chdir(self.workspace_location)
       
       with open(filename, 'w') as fasta:
            for row in x:
                fasta.write(">Labcode: " + row[4] + ", Marker: " + row[5] + " \n")
                if (row[8] == ""):
                    fasta.write("No sequences \n")
                else:
                    fasta.write(str(row[8]) + "\n")
        
       if (self.os_param == "windows"):
           os.startfile(filename)
       elif (self.os_param == "linux"):
           subprocess.call(["xdg-open", filename])
                
    
    def get_reference_dict(self, sample_loci_matrix):
        
        dict_reference = {}
        os.chdir(self.workspace_location)
        
        list_all = []
        list_project = []
        list_organism = []
        list_country = []
        list_locality = []
        
        sample_columnr = self.columnr_params[0]
        sequence_columnr = self.columnr_params[1]
        rownr_begin = self.rownr_params[0]
        rownr_end = self.rownr_params[1]
        
        rows_new = []
         
        if os.path.exists(self.reference_file):
            with open(self.reference_file, "r") as reference:
                for i, line in enumerate(reference, 1):
                    line = re.split(',|\t|;', line)
                    if (i < rownr_begin):
                        continue
                    if (line[sample_columnr-1] == ""):
                        break
                    if (i > rownr_end):
                        break
                    
                    sample_ref = line[int(self.columnr_params[0])-1]
                    for sample in sample_loci_matrix.keys():
                        if (sample_ref == sample):
                            print("Sample: " + sample)
                            dict_reference[sample_ref] = {}
                            dict_reference[sample_ref]["Project"] = line[int(self.columnr_params[3])-1] #Project
                            list_project.append(line[0])
                            dict_reference[sample_ref]["Organism"] = line[int(self.columnr_params[2])-1] #Organism
                            list_organism.append(line[1])
                            dict_reference[sample_ref]["Country"] = line[int(self.columnr_params[4])-1] #Country
                            list_country.append(line[5])
                            dict_reference[sample_ref]["Locality"] = line[int(self.columnr_params[5])-1] #Locality
                            list_locality.append(line[6])

            #with open(self.reference_file, "r") as reference:
                #csv_reader = csv.reader(reference, delimiter=';')
                #for i, line in enumerate(csv_reader):
                    #print(line)
                    #sample_ref = line[int(self.columnr_params[0])-1]
                    #for sample in sample_loci_matrix.keys():
                        #if (sample_ref == sample):
                            #print("Sample: " + sample)
                            #dict_reference[sample_ref] = {}
                            #dict_reference[sample_ref]["Project"] = line[int(self.columnr_params[3])-1] #Project
                            #list_project.append(line[0])
                            #dict_reference[sample_ref]["Organism"] = line[int(self.columnr_params[2])-1] #Organism
                            #list_organism.append(line[1])
                            #dict_reference[sample_ref]["Country"] = line[int(self.columnr_params[4])-1] #Country
                            #list_country.append(line[5])
                            #dict_reference[sample_ref]["Locality"] = line[int(self.columnr_params[5])-1] #Locality
                            #list_locality.append(line[6])
        
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
        
        print("")
        print("HIIER")
        print(sample_loci_matrix)
        print()
        
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
        
        columns = ('Project', 'Organism', 'Country', 'Locality', 'SampleName', 'LociName', 'AlleleIdx', 'AlleleLength', 'Sequence')
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
        tree_show.heading('Sequence', text='Sequence', anchor=tk.CENTER)
        
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
        
        print("Hier Reference dict:")
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
        
        count_element = ""
        unique_list = list(set(full_list))
        param_length = len(unique_list)
        first_row.append(param_length)
        for element in unique_list:
            second_row.append(element)
            count_element = full_list.count(element)
            first_row.append(count_element)
        
        print(full_list)
        print("Count ist")
        print(count_element)
        
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
        #print(dict_reference)
        
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
        #print(sequence_rows)
        print()
        number = db_instance.get_size_of_table()
        #print(number)
        number2 = db_instance.get_notempty_sequence_rows_number()
        #print(number2)
        
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
                                    #print(locus)
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
                                #print(line)
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
                            #print(loci_names)
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
        
        paramter_file_path = os.path.normpath(os.getcwd() + '/parameters.txt')
        if (os.path.exists(paramter_file_path)):
            with open(paramter_file_path) as parameters:
                for line in parameters:
                    if (not line) or line.startswith("#"):
                        continue
                    if (line.startswith("OutputFolder")):
                        output_folder = line.split("=")[1].rstrip("\r\n").strip()
                        if (not output_folder):
                            output_folder = "OutputDefault"
                        param_list_dirs.append(output_folder)
                        continue
                    if (line.startswith("Bin")):
                        bin_dir = line.split("=")[1].rstrip("\r\n").strip()
                        if (not bin_dir):
                            print("No name for data_dir in parameter file!")
                        param_list_dirs.append(os.path.normpath(bin_dir))
                        continue
                    if (line.startswith("RawData")):
                        raw_data_dir = line.split("=")[1].rstrip("\r\n").strip()
                        if (not raw_data_dir):
                            print("No name for data_dir in parameter file!")
                        param_list_dirs.append(os.path.normpath(raw_data_dir))
                        continue
                    
                    if (line.startswith("PrimerFile")):
                        primerfile = line.split("=")[1].rstrip("\r\n").strip()
                        if (not primerfile):
                            print("No name for primerfile in parameter file!")
                        param_list_files.append(os.path.normpath(primerfile))
                        continue
                    if (line.startswith("SampleFile")):
                        sample_sheet = line.split("=")[1].rstrip("\r\n").strip()
                        if (not sample_sheet):
                            print("No name for sample file in parameter file!")
                        param_list_files.append(os.path.normpath(sample_sheet))
                        continue
                    if (line.startswith("RExecutable")):
                        Rexecutable = line.split("=")[1].rstrip("\r\n").strip()
                        if (not Rexecutable):
                            Rexecutable = "R"
                        param_list_files.append(os.path.normpath(Rexecutable))
                        continue
                    if (line.startswith("AlleleList")):
                        allele_list = line.split("=")[1].rstrip("\r\n").strip()
                        param_list_files.append(os.path.normpath(allele_list))
                        continue
                
                    
                    
                    if (line.startswith("MaxMismatch")):
                        max_mismatches = line.split("=")[1].rstrip("\r\n").strip()
                        if (not max_mismatches):
                            max_mismatches = 2
                        else:
                            max_mismatches = int(max_mismatches)
                        calc_param_list.append(max_mismatches)
                        continue
                    if (line.startswith("MinCount")):
                        minCount = line.split("=")[1].rstrip("\r\n").strip()
                        if (not minCount):
                            minCount = 10
                        else:
                            minCount = int(minCount)
                        calc_param_list.append(minCount)
                        continue
                    if (line.startswith("MinLength")):
                        min_length = line.split("=")[1].rstrip("\r\n").strip()
                        if (not min_length):
                            min_length = 290
                        else:
                            min_length = int(min_length)
                        calc_param_list.append(min_length)
                        continue
                    if (line.startswith("ConsensusThreshold")):
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
                    
                    if (line.startswith("Ploidy")):
                        ploidy = line.split("=")[1].rstrip("\r\n").strip()
                        if (not ploidy):
                            ploidy = "diploid"
                        param_list_special.append(ploidy)
                        continue
                    if (line.startswith("Operatingsystem")):
                        par_os = line.split("=")[1].rstrip("\r\n").strip()
                        if (not par_os):
                            par_os = "windows"
                        param_list_special.append(par_os.lower())
                        continue
                    if (line.startswith("Uniqueidentifier")):
                        unique_name = line.split("=")[1].rstrip("\r\n").strip()
                        if (not unique_name):
                            unique_name = "Default"
                        param_list_special.append(unique_name)
                        continue
                    if (line.startswith("Indexcomboposition")):
                        index_combo = line.split("=")[1].rstrip("\r\n").strip()
                        if (not index_combo):
                            index_combo = 1
                        else:
                            index_combo = int(index_combo)
                        param_list_special.append(index_combo)
                        continue
                    
                    
                    
        else:
            print("There is no paramter file with the name parameters.txt")
            param_list_dirs = ["OutputDefault", "", ""]
            param_list_files = ["", "", "", "None"]
            param_list_special = ["diploid", "windows", "Default", "1"]
            calc_param_list = ["2", "20", "290", "0.7", "310 600"]
        
        return param_list_dirs, param_list_files, param_list_special, calc_param_list
    
    
    def get_general_params2(self):
        
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
        
        param_list_dirs = ["OutputFolder", "Bin", "RawData"]
        number1 = 0
        for i, param in enumerate(param_list_dirs):
            ctk.CTkLabel(dir_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(dir_params, corner_radius=5) #, textvariable=text
            if not((param == "OutputFolder")):
                browse_button = tk.Button(dir_params, text="Browse", command=lambda entry = param_entry: self.browse_dir(entry))
                browse_button.grid(column=2, row=i, padx=5, pady=5)
            elif (param == "OutputFolder"): #browse_dir_outputfolder
                browse_button = tk.Button(dir_params, text="Browse", command=lambda entry = param_entry: self.browse_dir_outputfolder(entry))
                browse_button.grid(column=2, row=i, padx=5, pady=5)
            if (len(self.var_list_dir) > 0):
                param_entry.insert(tk.END, os.path.normpath(self.var_list_dir[i]))
            else:
                param_entry.insert(tk.END, os.path.normpath(self.param_list_dirs[i]))
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            entries_dirs.append(param_entry)
            number1 += 1
        entries_total.append(entries_dirs)
            
        
        
        param_list_files = ["Primers", "SampleSheet", "Rexecutable", "AlleleList"]
        number2 = 0
        for i, param in enumerate(param_list_files):
            ctk.CTkLabel(file_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(file_params, corner_radius=5) #, textvariable=text
            browse_button = tk.Button(file_params, text="Browse", command=lambda entry=param_entry: self.browse_file_folder(entry))
            browse_button.grid(column=2, row=i, padx=5, pady=5)
            if (len(self.var_list_file) > 0):
                param_entry.insert(tk.END, os.path.normpath(self.var_list_file[i]))
            else:
                param_entry.insert(tk.END, os.path.normpath(self.param_list_files[i]))
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            entries_files.append(param_entry)
            number2 += 1
            
        
        
        entries_total.append(entries_files)
        
        
        
        param_list_calc = ["MaxMismatches", "MinCount", "MinLength", "ConsThreashold", "lengthWindow"]
        number3 = 0
        for i, param in enumerate(param_list_calc):
            ctk.CTkLabel(calc_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(calc_params, corner_radius=5) #, textvariable=text
            if (len(self.var_list_calc) > 0):
                param_entry.insert(tk.END, self.var_list_calc[i])
            else:
                param_entry.insert(tk.END, self.calc_param_list[i])
            param_entry.grid(column=1, row=i, padx=5, pady=5)
            #param_entry.insert(0, default_values[i])
            entries_calc.append(param_entry)
            number3 += 1
        entries_total.append(entries_calc)
        
        param_list_special = ["Ploidy", "OperatingSystem", "UniqueID", "IndexComboPosition"]
        number4 = 0
        for i, param in enumerate(param_list_special):
            ctk.CTkLabel(special_params, text=param + ":", width=20).grid(column=0, row=i, padx=5, pady=5)
            param_entry = ctk.CTkEntry(special_params, corner_radius=5) #, textvariable=text
            if (len(self.var_list_special) > 0):
                param_entry.insert(tk.END, self.var_list_special[i])
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
    
    def browse_dir_outputfolder(self, entry):
        directory_path = filedialog.askdirectory(initialdir=self.workspace_location, title="Select folder")
        directory_path = os.path.normpath(directory_path)
        output_path = ""
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
        
        #outputfolder_name = entry.get()
        
        entry.delete(0, tk.END)  # Clear the entry
        entry.insert(0, output_path)  # Insert the file/folder path into the entry
        #os.path.mkdir(output_path + "/" + outputfolder_name)
        
        #output_entry.bind("<Return>", lambda event: self.writesamplesheet(output_entry.get()))
        #output_entry.focus_set()
        
    
    def changeoutputfolder(self, outputfolder_path):
        self.output_dir = outputfolder_path
    
    def browse_dir(self, entry):
        directory_path = filedialog.askdirectory(initialdir=self.workspace_location, title="Select folder")
        directory_path = os.path.normpath(directory_path)
        output_path = ""
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
        entry.insert(0, output_path)  # Insert the file/folder path into the entry
        
    
    def browse_file_folder(self, entry):
        # Check if you want to select a file or folder (for simplicity, I'm using file selection here)
        file_path = filedialog.askopenfilename(initialdir=self.workspace_location, title="Select folder")
        file_path = os.path.normpath(file_path)
        output_path = ""
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
            
        
    def set_params2(self, entries):
        print("Setting them.")
        
        self.params = {entry[0] : entry[1].get() for entry in entries}
        
        # list_dir_params = entries[0]
        # list_file_params = entries[1]
        # list_calc_params = entries[2]
        # list_special_params = entries[3]
        
        # var_params_dir = []
        # for i, entry in enumerate(list_dir_params):
        #     var_params_dir.append(entry.get())
        
        # var_params_file = []
        # for i, entry in enumerate(list_file_params):
        #     var_params_file.append(entry.get())
        
        # var_params_calc = []
        # for i, entry in enumerate(list_calc_params):
        #     var_params_calc.append(entry.get())
        
        # var_params_special = []
        # for i, entry in enumerate(list_special_params):
        #     var_params_special.append(entry.get())
        
        # self.output_dir = var_params_dir[0]
        # self.bin_dir = var_params_dir[1]
        # self.rawdata_dir = var_params_dir[2]
        
        # self.primerfile = var_params_file[0]
        # self.sample_sheet = var_params_file[1]
        # self.Rexe = var_params_file[2]
        # self.allelelist = var_params_file[3]
        
        # self.max_mismatches = var_params_calc[0]
        # self.minCount = var_params_calc[1]
        # self.min_length = var_params_calc[2] 
        # self.alpha = var_params_calc[3]
        # self.lengthWindow = []
        # length_window_params = var_params_calc[4].rstrip().split(" ")
        # for i, param in enumerate(length_window_params, 1):
        #     if (i > 2):
        #         break
        #     else:
        #         self.lengthWindow.append(param)
        # #self.lengthWindow = var_params_calc[4]
        # print(self.lengthWindow)
        
        # self.ploidy = var_params_special[0]
        # self.os_param = var_params_special[1]
        # self.unique_prefix_name = var_params_special[2]
        # self.indexcombo_position = var_params_special[3]
        
        # self.var_list_dir = []
        # self.var_list_file = []
        # self.var_list_calc = []
        # self.var_list_special = []
        
        # self.var_list_dir.append(self.output_dir)
        # self.var_list_dir.append(self.bin_dir)
        # self.var_list_dir.append(self.rawdata_dir)
        
        # self.var_list_file.append(self.primerfile)
        # self.var_list_file.append(self.sample_sheet)
        # self.var_list_file.append(self.Rexe)
        # self.var_list_file.append(self.allelelist)
        
    
        # self.var_list_calc.append(self.max_mismatches)
        # self.var_list_calc.append(self.minCount)
        # self.var_list_calc.append(self.min_length)
        # self.var_list_calc.append(self.alpha)
        # self.var_list_calc.append(self.lengthWindow)
        
        # self.var_list_special.append(self.ploidy)
        # self.var_list_special.append(self.os_param)
        # self.var_list_special.append(self.unique_prefix_name)
        # self.var_list_special.append(self.indexcombo_position)
        
        #print()
        #print(self.rawdata_dir)
        #print(self.primerfile)
        #print(self.sample_sheet)
        #print(self.outputfolder)
        #print(self.ploidy)
    
       
    def write_params2(self, entries):
        
        self.set_params(entries)
        
        with open("parameters.txt", 'w') as file:
            file.write("###Parameters for the Pipeline \n")
            file.write("\n")
            file.write("#Name of the Folder where results of the pipeline are stored: \n")
            file.write("OutputFolder = " + os.path.normpath(self.params[""]) + "\n")
            file.write("#Path to the Bin Folder where Trimmo, Usearch and R executables are stored: \n")
            file.write("Bin = " + os.path.normpath(self.bin_dir) + "\n")
            file.write("#Path to the Folder containing raw fastq/fastq.gz files: \n")
            file.write("RawData = " + os.path.normpath(self.rawdata_dir) + "\n")
            file.write("\n")
            file.write("#Path to the File containing Primer information (.txt): \n")
            file.write("PrimerFile = " + os.path.normpath(self.primerfile) + " \n")
            file.write('#Path to the File containing sample and indexcombo information (.csv): \n')
            file.write("SampleFile = " + os.path.normpath(self.sample_sheet) + "\n")
            file.write('#Path to the Rexecutable (Only needed in windows) (.exe): \n')
            file.write("RExecutable = " + os.path.normpath(self.Rexe) + "\n")
            file.write('#Path to the Allelelist (if user wants to do listbased Call) (.txt): \n')
            file.write("AlleleList = " + os.path.normpath(self.allelelist) + " \n")
            file.write("\n")
            file.write('#Maximum number of mismatches between primer and amplicon sequences (integer; default=2) \n')
            file.write("MaxMismatch = " + str(self.max_mismatches) + " \n")
            file.write('#Minimum number of read count for an allele to be considered (integer; default = 20) \n')
            file.write("MinCount = " + str(self.minCount) + " \n")
            file.write('#Minimum amplicon length to be considered (integer; default=250) \n')
            file.write("MinLength= " + str(self.min_length) + " \n")
            file.write('#Consensus threshold (float; default = 0.7) \n')
            file.write("ConsensusThreshold = " + str(self.alpha) + " \n")
            #self.lengthWindow = self.lengthWindow.split()
            file.write("LengthWindow= " + str(self.lengthWindow[0]) + "," + str(self.lengthWindow[1]) + "\n")
            file.write("\n")
            file.write("#Ploidy level of the data (haploid/diploid; default = diploid): \n")
            file.write("Ploidy = " + str(self.ploidy) + "\n")
            file.write("#Operating System of the current user (default = windows): \n")
            file.write("Operatingsystem = " + str(self.os_param) + " \n")
            file.write("#Unique Identifier for each Sequence: \n")
            file.write("Uniqueidentifier = " + str(self.unique_prefix_name)+ " \n")
            file.write("#Position of the Indexcombination in the raw Fastqfilename (default = 1): \n")
            file.write("Indexcomboposition = " + str(self.indexcombo_position))
            file.write("\n")
            
        
        
        if self.os_param == "windows":
            os.startfile("parameters.txt")
        elif self.os_param == "linux":
            subprocess.run(["xdg-open", "parameters.txt"], check=True)
        
    def advanced_setting(self):
        print("Additional Settings:")
        
        
        
        if self.newWindow_advanced is not None and tk.Toplevel.winfo_exists(self.newWindow_advanced):
            self.newWindow_advanced.destroy()  # Destroy the existing window
        
        self.checkboxes_pipeline1 = []
        self.checkboxes_pipeline2 = []
        self.checkbox_advanced = []
        
        self.newWindow_advanced = tk.Toplevel(self)
        self.newWindow_advanced.title("Advanced Parameters:")
        self.newWindow_advanced.grid_rowconfigure(0, weight=1)
        self.newWindow_advanced.grid_rowconfigure(1, weight=1)
        self.newWindow_advanced.grid_columnconfigure(0, weight=1)
        self.newWindow_advanced.grid_rowconfigure(6, weight=1)
        self.newWindow_advanced.grid_rowconfigure(7, weight=1)
        
        #newWindow_advanced.grid_columnconfigure(1, weight=1)
        #newWindow_advanced.grid_columnconfigure(2, weight=1)
        
        
        
        
        advanced_input_settings = ctk.CTkLabel(self.newWindow_advanced, text="Advanced Input Settings:", font=ctk.CTkFont(size=20, weight="bold"))
        advanced_input_settings.grid(row=0, column=0, padx=20, pady=(20, 20))
        
        advanced_input_settings = ctk.CTkFrame(self.newWindow_advanced, width=100, corner_radius=2)
        advanced_input_settings.grid(row=1, column=0, rowspan=1, sticky="nsew")
        advanced_input_settings.grid_columnconfigure(0, weight=1)
        advanced_input_settings.grid_columnconfigure(1, weight=1)
        advanced_input_settings.grid_columnconfigure(2, weight=1)
        advanced_input_settings.grid_rowconfigure(0, weight=1)
        advanced_input_settings.grid_rowconfigure(1, weight=1)
        
        
        ctk.CTkLabel(advanced_input_settings, text="Adapters : ", width=20).grid(column=0, row=0, padx=5, pady=5)
        param_entry = ctk.CTkEntry(advanced_input_settings, corner_radius=5) #, textvariable=text
        browse_button = tk.Button(advanced_input_settings, text="Browse", command=lambda entry=param_entry: self.browse_file_folder(entry))
        browse_button.grid(column=2, row=0, padx=5, pady=5)
        if (self.adapterfile == ""):
            param_entry.insert(tk.END, os.path.normpath("TrueSeqAdaptersInUsage.fa"))
        else:
            param_entry.insert(tk.END, os.path.normpath(self.adapterfile))
        param_entry.grid(column=1, row=0, padx=5, pady=5)
        
        checkbox_performance_var = tk.IntVar(value = self.checkbox_performance)
        checkbox_performance = ctk.CTkCheckBox(advanced_input_settings, text="HighPerformance", variable=checkbox_performance_var)
        checkbox_performance.grid(row=1, column=0, padx=10, pady=(5, 5), sticky="w")
        self.checkbox_advanced.append(checkbox_performance)
        
        advanced_pipeline1_settings_logo = ctk.CTkLabel(self.newWindow_advanced, text="Advanced Pipeline Settings (1):", font=ctk.CTkFont(size=20, weight="bold"))
        advanced_pipeline1_settings_logo.grid(row=2, column=0, padx=20, pady=(20, 20))
        
        advanced_pipeline1_settings = ctk.CTkFrame(self.newWindow_advanced, width=100, corner_radius=2)
        advanced_pipeline1_settings.grid(row=3, column=0, rowspan=1, sticky="nsew")
        advanced_pipeline1_settings.grid_columnconfigure(0, weight=1)
        #advanced_pipeline1_settings.grid_columnconfigure(1, weight=1)
        #advanced_pipeline1_settings.grid_columnconfigure(2, weight=1)
        advanced_pipeline1_settings.grid_rowconfigure(0, weight=1)
        advanced_pipeline1_settings.grid_rowconfigure(1, weight=1)
        
        advanced_pipeline2_settings_logo = ctk.CTkLabel(self.newWindow_advanced, text="Advanced Pipeline Settings (2):", font=ctk.CTkFont(size=20, weight="bold"))
        advanced_pipeline2_settings_logo.grid(row=4, column=0, padx=20, pady=(20, 20))
        
        advanced_pipeline2_settings = ctk.CTkFrame(self.newWindow_advanced, width=100, corner_radius=2)
        advanced_pipeline2_settings.grid(row=5, column=0, rowspan=1, sticky="nsew")
        advanced_pipeline2_settings.grid_columnconfigure(0, weight=1)
        #advanced_pipeline2_settings.grid_columnconfigure(1, weight=1)
        #advanced_input_settings.grid_columnconfigure(2, weight=1)
        advanced_pipeline2_settings.grid_rowconfigure(0, weight=1)
        advanced_pipeline2_settings.grid_rowconfigure(1, weight=1)
        
        
        scrollframe_pipeline1 = ctk.CTkFrame(advanced_pipeline1_settings, width=50, corner_radius=2)
        scrollframe_pipeline1.grid(row=0, column=0, rowspan=1, sticky="nsew")
        scrollframe_pipeline1.grid_columnconfigure(0, weight=1)
        
        #checkboxes_organism_include = []
        #self.values_organism_include = ["organism1", "organism2", "organism3", "organism4"] #call from database
        options = ["Trimmomatic", "Usearch", "Demultiplexing", "Markerstatistics", "Markerplots+Markermatrix"]
        for i, option in enumerate(options):
            #state = self.checkbox_states_include_dict.get(value, 0)
    
            checkbox_var = tk.IntVar(value = self.checkbox_states_pipeline_advanced[option])
            checkbox = ctk.CTkCheckBox(scrollframe_pipeline1, text=option, variable=checkbox_var)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(5, 5), sticky="w")
            self.checkboxes_pipeline1.append(checkbox)
            #if (checkbox.get() == 1):
                #checkboxes_organism_include.append(checkbox.cget("text"))
        
        scrollframe_pipeline2 = ctk.CTkFrame(advanced_pipeline2_settings, width=50, corner_radius=2)
        scrollframe_pipeline2.grid(row=0, column=0, rowspan=1, sticky="nsew")
        scrollframe_pipeline2.grid_columnconfigure(0, weight=1)
        
        #checkboxes_organism_include = []
        #self.values_organism_include = ["organism1", "organism2", "organism3", "organism4"] #call from database
        options = ["Extract AlleleLengths", "ConsensusSequence", "ConsensusTogether", "Corrected", "AlleleCall"]
        for i, option in enumerate(options):
            #state = self.checkbox_states_include_dict.get(value, 0)
    
            checkbox_var = tk.IntVar(value = self.checkbox_states_pipeline_advanced[option])
            checkbox = ctk.CTkCheckBox(scrollframe_pipeline2, text=option, variable=checkbox_var)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(5, 5), sticky="w")
            self.checkboxes_pipeline2.append(checkbox)
        #param_entry.bind("<Return>", process_and_close)
        #param_entry.focus_set()
        
        start_pipeline1 = ctk.CTkButton(self.newWindow_advanced, border_width=1, border_color="black", text_color="black", text = "Run Pipeline1 Advanced", command = self.run_pipeline1_advanced)
        start_pipeline1.grid(row=6, column=0, sticky="sw", padx=(5, 25), pady=(10, 5))
        
        start_pipeline2 = ctk.CTkButton(self.newWindow_advanced, border_width=1, border_color="black", text_color="black", text = "Run Pipeline2 Advanced", command = self.run_pipeline2_advanced)
        start_pipeline2.grid(row=6, column=0, sticky="se", padx=(25, 5), pady=(10, 5))
        
        #write_samplesheet_button = ttk.Button(newWindow_advanced, text="Run Pipeline1 Advanced", command = self.run_pipeline1_advanced) #, command = lambda x = column_entries, y = row_entries, ref = param_entry_ref : self.set_samplesheetname_name(x, y, ref)
        #write_samplesheet_button.grid(column=0, row=6, sticky="se", padx=15, pady=4)
        
        #write_referenceparameters_button = ttk.Button(newWindow_advanced, text="Run Pipeline2 Advanced", command = self.run_pipeline2_advanced) #command = lambda x = column_entries, y = row_entries, ref = param_entry_ref : self.write_reference_parameterfile(x, y, ref)
        #write_referenceparameters_button.grid(column=0, row=6, sticky="sw", padx=15, pady=4)
        
        self.newWindow_advanced.protocol("WM_DELETE_WINDOW", self.on_closing_advanced)
    
    def run_pipeline1_advanced(self):
        for checkbox in self.checkboxes_pipeline1:
            checkbox_text = checkbox.cget("text")
            #print(checkbox_text)
            self.checkbox_states_pipeline_advanced[checkbox_text] = checkbox.get()
        
        for checkbox in self.checkboxes_pipeline2:
            checkbox_text = checkbox.cget("text")
            self.checkbox_states_pipeline_advanced[checkbox_text] = checkbox.get()
        
        
        self.report_dict["Trimmomatic"] = []
        self.report_dict["Usearch"] = []
        self.report_dict["Demultiplexing"] = []
        self.report_dict["MarkerStatistics"] = []
        self.report_dict["Markerplots+Matrix"] = []
        
        #self.report.append("Dictionary")
        
        
        new_window = tk.Toplevel(self)
        new_window.title("Pipeline Advanced (Length Detection)")
        new_window.grid_rowconfigure(0, weight=1)
        new_window.grid_columnconfigure(0, weight=1)
        #long_text = f"Running the first script in advanced mode: \n\n\n"
        
        text_font = ("Helvetica", 20)
        
        textbox = tk.Text(new_window, width=250, font=text_font)
        textbox.grid(row=5, column=0, padx=(5, 5), pady=(15, 15), sticky="nsew")
        textbox.grid_rowconfigure(0, weight=1)
        textbox.grid_columnconfigure(0, weight=1)
        #long_text += "Samplesheet File provides " + str(nrsamples) + " SampleIds. \n"
        #long_text += "Rawdata Folder provides " + str(nrfastq) + " zipped FastQfiles. \n"
        textbox.insert("0.0", "Results are saved in the Outputfolder : " + os.path.normpath(self.output_dir) + " \n\n\n")
        
        textbox.insert('end-1c', "High Performance Mode: " + str(self.performance) + " \n\n")
        textbox.insert('end-1c', "Running the Length Detection Process in advanced mode: \n\n\n\n\n")
        
        self.is_working_advanced1 = True
        thread = threading.Thread(target=self.run_pipeline1_advanced_thread, args=(textbox, ), daemon=True)
        thread.start()
        
        self.append_period_advanced(textbox)
    
    def append_period_advanced(self, textbox):
        # Append a period to the Text widget
        textbox.insert('end-1c', ".")
        
        # If the working function is still running, schedule this function to run again
        if self.is_working_advanced1:
            #self.after(1000, self.append_period_advanced(textbox))
            self.after(1000, lambda: self.append_period_advanced(textbox))
    
    def run_pipeline1_advanced_thread(self, textbox):
        
        #Amplicon_inst = amplicon_class(self.output_dir, os.path.normpath(self.primerfile), os.path.normpath(self.rawdata_dir), os.path.normpath(self.bin_dir), self.adapterfile, os.path.normpath(self.sample_sheet), self.max_mismatches, self.min_length, self.minCount, self.alpha, self.lengthWindow, self.ploidy, self.os_param, self.indexcombo_position)
        Amplicon_inst = length_class(self.performance)
        Amplicon_inst.parse_params()
        Amplicon_inst.set_outputs()
        Amplicon_inst.check_outputs()
        Amplicon_inst.set_executables()
        Amplicon_inst.check_executables()
        Amplicon_inst.set_samplelist()
        Amplicon_inst.set_primers()
        
        Amplicon_inst.checkInputDir()
        start = time.time()
        
        
        #self.is_working = True
        if (self.checkbox_states_pipeline_advanced["Trimmomatic"] == 1):
            nrsamples, samples = self.get_sample_number()
            nrfastq, fastqfiles = self.get_rawdata_number()
            nr_sample_fastq_compatible = self.get_compatibility(samples, fastqfiles)
            textbox.insert('end-1c', "Samplesheet File provides " + str(nrsamples) + " SampleIds. \n")
            textbox.insert('end-1c', "Rawdata Folder provides " + str(nrfastq) + " zipped FastQfiles. \n")
            textbox.insert('end-1c', "Trimmomatic started with " + str(nr_sample_fastq_compatible) + " fastqs from the raw datadir! \n")
            
            #Amplicon_inst.set_executables()
            #Amplicon_inst.set_samplelist()
            #Amplicon_inst.checkInputDir()
            errorcode, number, number_filled, number_empty, total_size = Amplicon_inst.runTrimomatic() #trimmed fastq files in /QC/
            print()
            message = "Success"
            if (errorcode != 0):
                message = "Failure"
            textbox.insert('end-1c', "Trimming sequences finished with " + message + " \n")
            textbox.insert('end-1c', "Results in the QC Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + " non-empty trimmed files!  \n\n\n\n")
            self.report_dict["Trimmomatic"].append(nr_sample_fastq_compatible)
            self.report_dict["Trimmomatic"].append(errorcode)
            self.report_dict["Trimmomatic"].append(number)
            self.report_dict["Trimmomatic"].append(number_filled)
            self.report_dict["Trimmomatic"].append(number_empty)
            self.report_dict["Trimmomatic"].append(total_size)
        else:
            textbox.insert('end-1c', "Trimmomatic skipped! \n\n\n\n")
        
        if (self.checkbox_states_pipeline_advanced["Usearch"] == 1):
            if (os.path.exists(os.path.normpath(self.output_dir + "/QC")) and os.listdir(self.output_dir + "/QC")):
                number, number_filled, number_empty, total_size = Amplicon_inst.check_results(os.path.normpath(self.output_dir + "/QC"))
                textbox.insert('end-1c', "Input for Merging is QC Folder which contains " + str(number_empty) + " empty and " +str(number_filled) + " non-empty trimmed files! \n\n")
                textbox.insert('end-1c', "Merging started with Usearch! \n")
                
                errorcode, number, number_filled, number_empty, total_size = Amplicon_inst.runUsearchMergePairs()
                message = "Success"
                if (errorcode != 0):
                    message = "Failure"
                textbox.insert('end-1c', "Merging finished with " + message + " \n")
                textbox.insert('end-1c', "Results in the MergedOut Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + " non-empty merged files!  \n\n\n\n")
                self.report_dict["Usearch"].append(errorcode)
                self.report_dict["Usearch"].append(number)
                self.report_dict["Usearch"].append(number_filled)
                self.report_dict["Usearch"].append(number_empty)
                self.report_dict["Usearch"].append(total_size)
            else:
                textbox.insert('end-1c', "But there is NO QC Folder or the QC Folder has no files so Merging will be skipped! \n\n\n\n\n")
        else:
            textbox.insert('end-1c', "Merging skipped! \n\n\n\n\n")
        
        if (self.checkbox_states_pipeline_advanced["Demultiplexing"] == 1):
            if (os.path.exists(os.path.normpath(self.output_dir + "/MergedOut")) and os.listdir(self.output_dir + "/MergedOut")):
                number, number_filled, number_empty, total_size = Amplicon_inst.check_results(os.path.normpath(self.output_dir + "/MergedOut"))
                textbox.insert('end-1c', "Input for Demultiplexing is MergedOut Folder which contains " + str(number_empty) + " empty and " +str(number_filled) + " non-empty trimmed files! \n\n")
                textbox.insert('end-1c', "Demultiplexing started! \n\n")
                number, number_filled, number_empty, total_size, primer_number, number_surviving, number_abandoned = Amplicon_inst.runDemultiplex()
                textbox.insert('end-1c', "Demultiplexing finished! \n")
                textbox.insert('end-1c', "Results in the SeparatOut Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + " non-empty demultiplexed files!  \n\n\n\n")
                #self.report_dict["Demultiplexing"].append(errorcode)
                self.report_dict["Demultiplexing"].append(number)
                self.report_dict["Demultiplexing"].append(number_filled)
                self.report_dict["Demultiplexing"].append(number_empty)
                self.report_dict["Demultiplexing"].append(primer_number)
                self.report_dict["Demultiplexing"].append(number_surviving)
                self.report_dict["Demultiplexing"].append(number_abandoned)
                self.report_dict["Demultiplexing"].append(total_size)
            else:
                textbox.insert('end-1c', "But there is no MergedOut folder or the MergedOut Folder has no files so Demultiplexing will be skipped ! \n\n\n\n")
        else:
            textbox.insert('end-1c', "Demultiplexing skipped! \n\n\n\n\n")
            
        
        if (self.checkbox_states_pipeline_advanced["Markerstatistics"] == 1):
            if (os.path.exists(os.path.normpath(self.output_dir + "/SeparatOut")) and os.listdir(self.output_dir + "/SeparatOut")):
                number, number_filled, number_empty, total_size = Amplicon_inst.check_results(os.path.normpath(self.output_dir + "/SeparatOut"))
                textbox.insert('end-1c', "Input for Markerstatistics is SeparatOut Folder which contains " + str(number_empty) + " empty and " +str(number_filled) + " non-empty demultiplexed files! \n\n")
                textbox.insert('end-1c', "Finding sequence lengths and counts per sample-locus combination! \n\n")
                lowest_length, filename_lowest, largest_length, filename_largest, number, number_filled, number_empty, total_size = Amplicon_inst.getLengthStatistics()
                textbox.insert('end-1c', "Markerstatistics finished! \n")
                textbox.insert('end-1c', "Results in the MarkerStatistics Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + " non-empty statistics files!  \n")
                textbox.insert('end-1c', "Smallest length is " + str(lowest_length) + " from the file " + os.path.basename(os.path.normpath(filename_lowest)) + " and largest length is " + str(largest_length) + " from the file " + os.path.basename(os.path.normpath(filename_largest)) + "!  \n\n\n\n")
                #textbox.insert('end-1c', "The length " + str(length_highest_countfrequency) + " has the highest count frequency of " + str(highest_countfrequency) + "% in the file " + os.path.basename(os.path.normpath(filename_highestfrequence)) + " from all files!  \n\n\n\n")
                self.report_dict["MarkerStatistics"].append(number)
                self.report_dict["MarkerStatistics"].append(number_filled)
                self.report_dict["MarkerStatistics"].append(number_empty)
                self.report_dict["MarkerStatistics"].append(lowest_length)
                self.report_dict["MarkerStatistics"].append(os.path.basename(os.path.normpath(filename_lowest)))
                self.report_dict["MarkerStatistics"].append(largest_length)
                self.report_dict["MarkerStatistics"].append(os.path.basename(os.path.normpath(filename_largest)))
                self.report_dict["MarkerStatistics"].append(total_size)
                #if (0 in self.lengthWindow):
                    #self.lengthWindow = []
                    #self.lengthWindow.append(lowest_length)
                    #self.lengthWindow.append(largest_length)
            
            else:
                textbox.insert('end-1c', "But there is no SeparatOut folder or the SeparatOut Folder has no files so MarkerStatistics will be skipped ! \n\n\n\n")
        else:
            textbox.insert('end-1c', "MarkerStatistics skipped! \n\n\n\n\n")
            
        if (self.checkbox_states_pipeline_advanced["Markerplots+Markermatrix"] == 1):
            if (os.path.exists(os.path.normpath(self.output_dir + '/MarkerStatistics')) and os.listdir(self.output_dir + '/MarkerStatistics')):
                number, number_filled, number_empty, total_size = Amplicon_inst.check_results(os.path.normpath(self.output_dir + '/MarkerStatistics'))
                textbox.insert('end-1c', "Input for Markerplots is MarkerStatistics Folder which contains " + str(number_empty) + " empty and " +str(number_filled) + " non-empty .statistics files! \n\n")
                textbox.insert('end-1c', "Get Markerplots (Counting Frequency vs sequence length) and Sample-Marker Length Matrix ! \n\n\n\n")
                errorcode = Amplicon_inst.GenotypeLength()
                message = "Success"
                if (errorcode != 0):
                    message = "Failure"
                textbox.insert('end-1c', "Making Markerplots and Matrix ended with " + message + "! \n\n\n\n\n")
                
                self.report_dict["Markerplots+Matrix"].append(message)
                
                markerplots_path = os.path.normpath(self.output_dir + "/Markerplots/markerplots.pdf")
                if (message == "Success"):
                    if (self.os_param == "windows"):
                        os.startfile(markerplots_path)
                    elif (self.os_param == "linux"):
                        subprocess.call(["xdg-open", markerplots_path])
                        
            else:
                textbox.insert('end-1c', "But there is no Markerstatistics folder or the Markerstatistics Folder has no files so MarkerStatistics will be skipped ! \n\n")
        else:
            textbox.insert('end-1c', "Markerplots skipped! \n\n\n\n")
        
        end = time.time()
        textbox.insert('end-1c', "Time spent: " + str(end - start) + " \n\n")
        #"Markerstatistics", "Markerplots+Markermatrix"
        #lowest_length, filename_lowest, largest_length, filename_largest, number, number_filled, number_empty = Amplicon_inst.getLengthStatistics()
        
        self.is_working_advanced1 = False
    
    def run_pipeline2_advanced(self):
        print("Hei")
        for checkbox in self.checkboxes_pipeline1:
            checkbox_text = checkbox.cget("text")
            #print(checkbox_text)
            self.checkbox_states_pipeline_advanced[checkbox_text] = checkbox.get()
        
        for checkbox in self.checkboxes_pipeline2:
            checkbox_text = checkbox.cget("text")
            self.checkbox_states_pipeline_advanced[checkbox_text] = checkbox.get()
        
        
        
        
        self.report_dict["AllelesOut"] = []
        self.report_dict["ConsensusOut"] = []
        self.report_dict["ConsensusTogether"] = []
        self.report_dict["Corrected"] = []
        self.report_dict["AlleleCall"] = []
        
        
        new_window = tk.Toplevel(self)
        new_window.title("Pipeline Advanced (SNP Detection)")
        new_window.grid_rowconfigure(0, weight=1)
        new_window.grid_columnconfigure(0, weight=1)
        #long_text = f"Running the first script in advanced mode: \n\n\n"
        
        text_font = ("Helvetica", 20)
        
        textbox = tk.Text(new_window, width=250, font=text_font)
        textbox.grid(row=5, column=0, padx=(5, 5), pady=(15, 15), sticky="nsew")
        textbox.grid_rowconfigure(0, weight=1)
        textbox.grid_columnconfigure(0, weight=1)
        #long_text += "Samplesheet File provides " + str(nrsamples) + " SampleIds. \n"
        #long_text += "Rawdata Folder provides " + str(nrfastq) + " zipped FastQfiles. \n"
        textbox.insert("0.0", "Results are saved in the Outputfolder : " + os.path.normpath(self.output_dir) + " \n\n\n")
        
        textbox.insert('end-1c', "Running the SNP Detection Process in advanced mode: \n\n\n\n\n")
        
        self.is_working_advanced1 = True
        thread = threading.Thread(target=self.run_pipeline2_advanced_thread, args=(textbox, ), daemon=True)
        thread.start()
        
        self.append_period_advanced(textbox)
    
    
    def run_pipeline2_advanced_thread(self, textbox):
        
        Sequence_inst = class_quality(self.output_dir, os.path.normpath(self.primerfile), os.path.normpath(self.sample_sheet), self.minCount, self.alpha, self.ploidy, self.unique_prefix_name, self.os_param, self.allelelist, self.textbox_bool)
        
        #self.report_dict["AllelesOut"] = []
        #self.report_dict["ConcensusOut"] = []
        #self.report_dict["ConcensusTogether"] = []
        #self.report_dict["Corrected"] = []
        #self.report_dict["AlleleCall"] = []
        
        samples = []
    
        if (self.checkbox_states_pipeline_advanced["Extract AlleleLengths"] == 1):
            if (os.path.exists(os.path.normpath(self.output_dir + "/MarkerPlots")) and os.listdir(self.output_dir + "/MarkerPlots") or (os.path.exists(os.path.normpath(self.output_dir + "/SeparatOut")) and os.listdir(self.output_dir + "/SeparatOut"))):
                number, number_filled, number_empty = Sequence_inst.check_results(os.path.normpath(self.output_dir + '/SeparatOut'))
                textbox.insert('end-1c', "Extracting AlleleLengths from the Markermatrix in Markerplots Folder! \n")
                textbox.insert('end-1c', "Getting corresponding sequences from the SeparatOutFolder! \n")
                textbox.insert('end-1c', "SeparatOutFolder contains " + str(number_empty) + " empty and " + str(number_filled) + " non-empty demultiplexed files!  \n")
                
                if (self.ploidy == "diploid"):
                    allele_length_dict, number2, number_filled2, number_empty2, samples_without_lengths, samples, sample_duplicates = Sequence_inst.extract_AlleleLengths_diploid() 
                elif (self.ploidy == "haploid"):
                    allele_length_dict, number2, number_filled2, number_empty2, samples_without_lengths, samples, sample_duplicates = Sequence_inst.extract_AlleleLengths_haploid() #trimmed fastq files in /QC/
                print()
                textbox.insert('end-1c', "Extracting Allele Lengths and corresponding sequences finished! \n")
                textbox.insert('end-1c', "Results in the AllelesOut Folder which contains " + str(number_empty2) + " empty and " + str(number_filled2) + "non-empty files!  \n")
                textbox.insert('end-1c', "Number of Samples with 0 AlleleLengths: " + str(len(samples_without_lengths)) + "  \n")
                textbox.insert('end-1c', "Number of Sample Duplicates in Markermatrix: " + str(len(sample_duplicates)) + "  \n\n\n\n\n")
                self.report_dict["AllelesOut"].append(number2)
                self.report_dict["AllelesOut"].append(number_filled2)
                self.report_dict["AllelesOut"].append(number_empty2)
                self.report_dict["AllelesOut"].append(len(samples_without_lengths))
            else:
                textbox.insert('end-1c', "InputFolders are Markerplots and SeparatOut but there is NO Markerplots Folder or no SeparatOut Folder or SeparatOutFolder has no files so this process will be skipped! \n\n\n\n\n")
        else:
            textbox.insert('end-1c', "Extracting AlleleLengths skipped! \n\n\n\n\n")
        
        
        if (self.checkbox_states_pipeline_advanced["ConsensusSequence"] == 1):
            if (os.path.exists(os.path.normpath(self.output_dir + "/AllelesOut")) and os.listdir(self.output_dir + "/AllelesOut")):
                number, number_filled, number_empty = Sequence_inst.check_results(os.path.normpath(self.output_dir + '/AllelesOut'))
                textbox.insert('end-1c', "Input for the ConsensusSequence process is the AllelesOut Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + " files for each sample-locus-Lenght combination!  \n")
                textbox.insert('end-1c', "Making ConsensusSequences form each Sample-Locus-Length combination from the AllelesOut Folder! \n")
                
                number, number_filled, number_empty = Sequence_inst.RunConsensusAll()
                
                
                textbox.insert('end-1c', "Making ConsensusSequences finished! \n")
                textbox.insert('end-1c', "Results in the ConsensusOut Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + " non-empty files!  \n\n\n\n\n")
                self.report_dict["ConsensusOut"].append(number)
                self.report_dict["ConsensusOut"].append(number_filled)
                self.report_dict["ConsensusOut"].append(number_empty)
            else:
                textbox.insert('end-1c', "InputFolder is AllelesOut but there is NO AllelesOut Folder or the AllelesOut Folder is empty and the ConsensusSequences Process will be skipped! \n\n\n\n\n")
        else:
            textbox.insert('end-1c', "Making ConsensusSequences skipped! \n\n\n\n\n")
        
        
        if (self.checkbox_states_pipeline_advanced["ConsensusTogether"] == 1):
            if (os.path.exists(os.path.normpath(self.output_dir + "/ConsensusOut")) and os.listdir(self.output_dir + "/ConsensusOut")):
                number, number_filled, number_empty = Sequence_inst.check_results(os.path.normpath(self.output_dir + '/ConsensusOut'))
                textbox.insert('end-1c', "Input for the ConsensusTogether process is the AllelesOut Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + " files for each sample-locus-Lenght combination!  \n")
                textbox.insert('end-1c', "Combining all ConsensusSequences from each Sample-Locus-Length combination per Locus! \n")
                
                loci_number, loci = Sequence_inst.extract_marker_diploid()
                number, number_filled, number_empty = Sequence_inst.joinSamplesSameMarker()
                
                
                textbox.insert('end-1c', "ConsensusTogether process finished! \n")
                textbox.insert('end-1c', "Results in the ConsensuTogether Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + "non-empty files!  \n\n\n\n\n")
                self.report_dict["ConsensusTogether"].append(number)
                self.report_dict["ConsensusTogether"].append(number_filled)
                self.report_dict["ConsensusTogether"].append(number_empty)
            else:
                textbox.insert('end-1c', "InputFolder is ConsensusOut but there is NO ConsensusOut Folder or the ConsensusOut Folder is empty and the ConsensusTogether Process will be skipped! \n\n\n\n\n")
        else:
            textbox.insert('end-1c', "ConsensusTogether process skipped! \n\n\n\n\n")
            
        
        if (self.checkbox_states_pipeline_advanced["Corrected"] == 1):
            if (os.path.exists(os.path.normpath(self.output_dir + '/ConsensusTogether')) and os.listdir(self.output_dir + '/ConsensusTogether')):
                number, number_filled, number_empty = Sequence_inst.check_results(os.path.normpath(self.output_dir + '/ConsensusTogether'))
                textbox.insert('end-1c', "Input for the Correction process is the CosensusTogether Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + " files for each locus!  \n")
                textbox.insert('end-1c', "Correcting all sequences for each file in the ConsensusTogether Folder! \n")
                
                number, number_filled, number_empty = Sequence_inst.correctAllSeqs()
                
                
                textbox.insert('end-1c', "Correction process finished! \n")
                textbox.insert('end-1c', "Results in the Corrected Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + "non-empty files!  \n\n\n\n\n")
                self.report_dict["Corrected"].append(number)
                self.report_dict["Corrected"].append(number_filled)
                self.report_dict["Corrected"].append(number_empty)
            else:
                textbox.insert('end-1c', "InputFolder is ConsensusTogether folder but there is NO ConsensusTogether Folder or the ConsensusTogether folder is empty and the Corection Process will be skipped! \n\n\n\n\n")
        else:
            textbox.insert('end-1c', "Correction process skipped! \n\n\n\n\n")
            
        
        if (self.checkbox_states_pipeline_advanced["AlleleCall"] == 1):
            if (os.path.exists(os.path.normpath(self.output_dir + '/Corrected')) and os.listdir(self.output_dir + '/Corrected')):
                number, number_filled, number_empty = Sequence_inst.check_results(os.path.normpath(self.output_dir + '/Corrected'))
                textbox.insert('end-1c', "Input for the AlleleCall process is the Corrected Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + " files for each locus!  \n")
                if (self.allelelist == "" or os.path.basename(os.path.normpath(self.allelelist)) == "." or not(os.path.exists(self.allelelist))):
                    textbox.insert('end-1c', "Denovo Call for finding SNPs! \n")
                else:
                    textbox.insert('end-1c', "Listbased Call for finding SNPs! \n")
                
                sample_number, samples = Sequence_inst.extract_samples()
                print("Finished getting samples")
                Sequence_inst.extract_loci()
                print("Finished getting loci")
                additional_loci, markers_more_than2 = Sequence_inst.AlleleCall()
                print("Finished calling Alleles")
                
                
                textbox.insert('end-1c', "AlleleCall process finished! \n")
                textbox.insert('end-1c', "We have " + str(len(additional_loci)) + " Additional loci from previous list: \n")
                textbox.insert('end-1c', str(additional_loci) + " \n")
                textbox.insert('end-1c', "We have " + str(len(markers_more_than2)) + " markers with 'MC': \n")
                textbox.insert('end-1c', str(markers_more_than2) + " \n")
                textbox.insert('end-1c', "Results in the AlleleCall Folder!  \n\n\n\n\n")
                #self.report_dict["AlleleCall"].append(number)
                #self.report_dict["AlleleCall"].append(number_filled)
                #self.report_dict["AlleleCall"].append(number_empty)
                allelematrix = os.path.normpath(self.output_dir + "/AlleleCall//matrix.txt")
                if (self.os_param == "windows"):
                    os.startfile(allelematrix)
                elif (self.os_param == "linux"):
                    subprocess.call(["xdg-open", allelematrix])
            else:
                textbox.insert('end-1c', "InputFolder is Corrected folder but there is NO Corrected Folder or the Corrected folder is empty and the AlleleCall process will be skipped! \n\n\n\n\n")
        else:
            textbox.insert('end-1c', "AlleleCall process skipped! \n\n\n\n\n")
        
        
            
            
        #"Markerstatistics", "Markerplots+Markermatrix"
        #lowest_length, filename_lowest, largest_length, filename_largest, number, number_filled, number_empty = Amplicon_inst.getLengthStatistics()
        
        self.is_working_advanced1 = False
    
    
    def get_report(self):
        #self.report_dict["AllelesOut"] = []
        #self.report_dict["ConsensusOut"] = []
        #self.report_dict["ConsensusTogether"] = []
        #self.report_dict["Corrected"] = []
        #self.report_dict["AlleleCall"] = []
        #self.advanced_pipeline_options = ["Trimmomatic", "Usearch", "Demultiplexing", "Markerstatistics", "Markerplots+Markermatrix",  "Extract AlleleLengths", "ConsensusSequence", "ConsensusTogether", "Corrected", "AlleleCall"]
        with open("Report.txt", 'w') as file:
            file.write("###Report of the pipeline process: \n\n")
            if (self.checkbox_states_pipeline_advanced["Trimmomatic"] == 1):
                file.write("#Trimmomatic took as input " + str(self.report_dict["Trimmomatic"][0]) + " fastq files and ended with \n")
                file.write(str(self.report_dict["Trimmomatic"][3]) + " non-empty and " + str(self.report_dict["Trimmomatic"][4]) + " empty trimmed files. \n\n")
                file.write("Total size of all result files: " + str(self.report_dict["Trimmomatic"][5]) + " Bytes. \n\n")
            if (self.checkbox_states_pipeline_advanced["Usearch"] == 1):
                #file.write("#Usearch took as input " + str(self.report_dict["Trimmomatic"][3]) + " trimmed files and ended with \n")
                file.write("Usearch finished with " + str(self.report_dict["Usearch"][2]) + " non-empty and " + str(self.report_dict["Usearch"][3]) + " empty merged files. \n\n")
                file.write("Total size of all result files: " + str(self.report_dict["Usearch"][4]) + " Bytes. \n\n")
            if (self.checkbox_states_pipeline_advanced["Demultiplexing"] == 1):
                #file.write("#Usearch took as input " + str(self.report_dict["Trimmomatic"][3]) + " trimmed files and ended with \n")
                file.write("Demultiplexing finished with " + str(self.report_dict["Demultiplexing"][1]) + " non-empty and " + str(self.report_dict["Demultiplexing"][2]) + " empty demultiplexed files. \n")
                file.write("Furthermore " + str(self.report_dict["Demultiplexing"][4]) + " sequences survived and " + str(self.report_dict["Demultiplexing"][5]) + " sequences were abandoned. \n")
                file.write("Total size of all result files: " + str(self.report_dict["Demultiplexing"][6]) + " Bytes. \n\n")
            if (self.checkbox_states_pipeline_advanced["Markerstatistics"] == 1):
                file.write("MarkerStatistics finished with " + str(self.report_dict["MarkerStatistics"][1]) + " non-empty and " + str(self.report_dict["MarkerStatistics"][2]) + " empty demultiplexed files. \n")
                file.write("Smallest length is " + str(self.report_dict["MarkerStatistics"][3]) + " from the file " + str(self.report_dict["MarkerStatistics"][4]) + " and largest length is " + str(self.report_dict["MarkerStatistics"][5]) + " from the file " + str(self.report_dict["MarkerStatistics"][6]) + "!  \n\n")
                file.write("Total size of all result files: " + str(self.report_dict["MarkerStatistics"][7]) + " Bytes. \n\n")
            if (self.checkbox_states_pipeline_advanced["Markerplots+Markermatrix"] == 1):
                file.write("Generating Markerplots and Markermatrix finished with " + str(self.report_dict["Markerplots+Matrix"][0]) + "!  \n\n")
            if (self.checkbox_states_pipeline_advanced["Extract AlleleLengths"] == 1):
                file.write("Extracting Allelelengths and corresponding sequences for each locus-sample combination finished with " + str(self.report_dict["AllelesOut"][1]) + " non-empty and " + str(self.report_dict["AllelesOut"][2]) + " empty files! \n")
                file.write(str(self.report_dict["AllelesOut"][3]) + " Samples " + " have no Lengths! \n\n")
            if (self.checkbox_states_pipeline_advanced["ConsensusSequence"] == 1):
                file.write("ConsensesSequence process finished with " + str(self.report_dict["ConsensusOut"][1]) +  " non-empty and " + str(self.report_dict["ConsensusOut"][2]) + " empty files! \n\n")
            if (self.checkbox_states_pipeline_advanced["ConsensusTogether"] == 1):
                file.write("ConsensesTogether process finished with " + str(self.report_dict["ConsensusTogether"][1]) +  " non-empty and " + str(self.report_dict["ConsensusTogether"][2]) + " empty files! \n\n")
            if (self.checkbox_states_pipeline_advanced["Corrected"] == 1):
                file.write("Correction process finished with " + str(self.report_dict["Corrected"][1]) +  " non-empty and " + str(self.report_dict["Corrected"][2]) + " empty files! \n\n")
                
                
                
        if (self.os_param == "windows"):
            os.startfile("Report.txt")
        elif (self.os_param == "linux"):
            subprocess.call(["xdg-open", "Report.txt"])
    
    def on_closing_advanced(self):
        print()
        for checkbox in self.checkboxes_pipeline1:
            checkbox_text = checkbox.cget("text")
            #print(checkbox_text)
            self.checkbox_states_pipeline_advanced[checkbox_text] = checkbox.get()
        
        for checkbox in self.checkboxes_pipeline2:
            checkbox_text = checkbox.cget("text")
            self.checkbox_states_pipeline_advanced[checkbox_text] = checkbox.get()
        
        for checkbox in self.checkbox_advanced:
            #checkbox_text = checkbox.cget("text")
            self.checkbox_performance = checkbox.get()
        
        print(self.checkbox_performance)
        if (self.checkbox_performance):
            self.performance = True
        else:
            self.performance = False
        
        
        self.newWindow_advanced.destroy()
        
        
    def set_adapters(self, x):
        self.adapterfile = os.path.basename(os.path.normpath(x))
        print(self.adapterfile)
        
    def run_first_script(self):
        
        self.report_dict["Trimmomatic"] = []
        self.report_dict["Usearch"] = []
        self.report_dict["Demultiplexing"] = []
        self.report_dict["MarkerStatistics"] = []
        self.report_dict["Markerplots+Matrix"] = []
        
        self.checkbox_states_pipeline_advanced["Trimmomatic"] = 1
        self.checkbox_states_pipeline_advanced["Usearch"] = 1
        self.checkbox_states_pipeline_advanced["Demultiplexing"] = 1
        self.checkbox_states_pipeline_advanced["Markerstatistics"] = 1
        self.checkbox_states_pipeline_advanced["Markerplots+Markermatrix"] = 1
        
        self.textbox.delete(0.0, 'end')
        self.textbox.insert('end-1c', "Start First Script! \n")
        
        thread = threading.Thread(target=self.run_amplicon_methods)
        thread.start()
    
        # Start appending periods to the textbox
        self.append_period()
        
        
        
    def run_amplicon_methods(self):
        
        self.textbox.delete(0.0, 'end')
        
        start = time.time()
        self.textbox.insert('end-1c', "Start First Script! \n")
        #Amplicon_inst = amplicon_class(self.output_dir, os.path.normpath(self.primerfile), os.path.normpath(self.rawdata_dir), os.path.normpath(self.bin_dir), self.adapterfile, os.path.normpath(self.sample_sheet), self.max_mismatches, self.min_length, self.minCount, self.alpha, self.lengthWindow, self.ploidy, self.os_param, self.indexcombo_position)
        Amplicon_inst = length_class(self.performance)
        Amplicon_inst.parse_params()
        Amplicon_inst.set_outputs()
        Amplicon_inst.check_outputs()
        Amplicon_inst.set_executables()
        Amplicon_inst.check_executables()
        Amplicon_inst.set_samplelist()
        Amplicon_inst.set_primers()
        
        output, samplesheet, primerfile, rawdata = Amplicon_inst.return_params()
        
        Amplicon_inst.checkInputDir()
        
        #global is_working
        self.is_working = True
        
        self.textbox.insert('end-1c', "\n")
        self.textbox.insert('end-1c', "OutputFolder : " + output + " \n")
        self.textbox.insert('end-1c', "SampleSheet : " + samplesheet + " \n")
        self.textbox.insert('end-1c', "Primerfile : " + primerfile + " \n")
        self.textbox.insert('end-1c', "RawDataFolder : " + rawdata + " \n")
        self.textbox.insert('end-1c', "\n")
        
        self.append_period()
        
        
        
        nrsamples, samples = self.get_sample_number()
        nrfastq, fastqfiles = self.get_rawdata_number()
        nr_sample_fastq_compatible = self.get_compatibility(samples, fastqfiles)
        print()
        self.textbox.insert('end-1c', "Trimming started \n")
        errorcode, number, number_filled, number_empty, total_size = Amplicon_inst.runTrimomatic() #trimmed fastq files in /QC/
        print()
        message = "Success"
        if (errorcode != 0):
            message = "Failure"
        self.textbox.insert('end-1c', "Trimming sequences finished with " + message + " \n")
        self.report_dict["Trimmomatic"].append(nr_sample_fastq_compatible)
        self.report_dict["Trimmomatic"].append(errorcode)
        self.report_dict["Trimmomatic"].append(number)
        self.report_dict["Trimmomatic"].append(number_filled)
        self.report_dict["Trimmomatic"].append(number_empty)
        self.report_dict["Trimmomatic"].append(total_size)
        
        print()
        self.textbox.insert('end-1c', " \n")
        self.textbox.insert('end-1c', "Merging started \n")
        errorcode, number, number_filled, number_empty, total_size = Amplicon_inst.runUsearchMergePairs() #merged fastq files in /MergedOut/
        print()
        message = "Success"
        if (errorcode != 0):
            message = "Failure"
        self.textbox.insert('end-1c', "Merging sequences finished with " + message + " \n")
        self.report_dict["Usearch"].append(errorcode)
        self.report_dict["Usearch"].append(number)
        self.report_dict["Usearch"].append(number_filled)
        self.report_dict["Usearch"].append(number_empty)
        self.report_dict["Usearch"].append(total_size)
        
        print()
        self.textbox.insert('end-1c', "Demultiplexing started \n")
        number, number_filled, number_empty, total_size, primer_number, number_surviving, number_abandoned = Amplicon_inst.runDemultiplex() #demultiplexed fasta files in /Separatout/ with name convention: newSampleName - Locus/PrimerName
        print()
        self.textbox.insert('end-1c', "Demultiplexing finished with " + str(number_empty) + " empty and " + str(number_filled) + " non-empty demultiplexed files!  \n")
        #errorcode, number, number_filled, number_empty, primer_number, number_surviving, number_abandoned = Amplicon_inst.runDemultiplex()
        #self.report_dict["Demultiplexing"].append(errorcode)
        self.report_dict["Demultiplexing"].append(number)
        self.report_dict["Demultiplexing"].append(number_filled)
        self.report_dict["Demultiplexing"].append(number_empty)
        self.report_dict["Demultiplexing"].append(primer_number)
        self.report_dict["Demultiplexing"].append(number_surviving)
        self.report_dict["Demultiplexing"].append(number_abandoned)
        self.report_dict["Demultiplexing"].append(total_size)
        
        
        lowest_length, filename_lowest, largest_length, filename_largest, number, number_filled, number_empty, total_size = Amplicon_inst.getLengthStatistics()
        self.textbox.insert('end-1c', "Markerstatistics finished! \n")
        self.textbox.insert('end-1c', "Results in the MarkerStatistics Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + " non-empty statistics files!  \n")
        self.textbox.insert('end-1c', "Smallest length is " + str(lowest_length) + " from the file " + os.path.basename(os.path.normpath(filename_lowest)) + " and largest length is " + str(largest_length) + " from the file " + os.path.basename(os.path.normpath(filename_largest)) + "!  \n\n\n\n")
        #textbox.insert('end-1c', "The length " + str(length_highest_countfrequency) + " has the highest count frequency of " + str(highest_countfrequency) + "% in the file " + os.path.basename(os.path.normpath(filename_highestfrequence)) + " from all files!  \n\n\n\n")
        self.report_dict["MarkerStatistics"].append(number)
        self.report_dict["MarkerStatistics"].append(number_filled)
        self.report_dict["MarkerStatistics"].append(number_empty)
        self.report_dict["MarkerStatistics"].append(lowest_length)
        self.report_dict["MarkerStatistics"].append(filename_lowest)
        self.report_dict["MarkerStatistics"].append(largest_length)
        self.report_dict["MarkerStatistics"].append(filename_largest)
        self.report_dict["MarkerStatistics"].append(total_size)
        
        #if (0 in self.lengthWindow):
            #self.lengthWindow = []
            #self.lengthWindow.append(lowest_length)
            #self.lengthWindow.append(largest_length)
        
        
        #Amplicon_inst.write_to_analysis_lengths()
        
        
        
        print()
        self.textbox.insert('end-1c', "Finding Length-based Genotypes started \n")
        error_code = Amplicon_inst.GenotypeLength()
        print()
        message = "Success"
        if (errorcode != 0):
            message = "Failure"
        self.textbox.insert('end-1c', "Getting Genotype Lengths finished with " + message + " \n\n")
        self.report_dict["Markerplots+Matrix"].append(message)
        #lengths = Amplicon_inst.get_unique_lengths()
        #self.textbox.insert('end-1c', "Found  " + str(lengths) + " unique lengths. \n")#get_unique_lengths
        #self.textbox.insert('end-1c', "Getting Genotype Lengths finished with \n")
        
        markerplots_path = os.path.normpath(self.output_dir + "/Markerplots/markerplots.pdf")
        if (message == "Success"):
            if (self.os_param == "windows"):
                os.startfile(markerplots_path)
            elif (self.os_param == "linux"):
                subprocess.call(["xdg-open", markerplots_path])
        
        self.is_working = False
        
        end = time.time()
        self.textbox.insert('end-1c', "Time: " + str(end - start) + " \n")
        
        self.status_review_script1()
        
    def status_review_script1(self):
        print()
        
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
        
        self.report_dict["AllelesOut"] = []
        self.report_dict["ConsensusOut"] = []
        self.report_dict["ConsensusTogether"] = []
        self.report_dict["Corrected"] = []
        self.report_dict["AlleleCall"] = []
        
        self.checkbox_states_pipeline_advanced["Extract AlleleLengths"] = 1
        self.checkbox_states_pipeline_advanced["ConsensusSequence"] = 1
        self.checkbox_states_pipeline_advanced["ConsensusTogether"] = 1
        self.checkbox_states_pipeline_advanced["Corrected"] = 1
        self.checkbox_states_pipeline_advanced["AlleleCall"] = 1
        
        self.textbox.delete(0.0, 'end')
        self.textbox.insert('end-1c', "Start Second Script! \n")
        
        thread = threading.Thread(target=self.run_sequence_methods)
        thread.start()
    
        # Start appending periods to the textbox
        self.append_period()
    
    
    
    def run_sequence_methods(self):
        
        self.textbox.delete(0.0, 'end')
        self.textbox.insert('end-1c', "Start Second Script! \n")
        
        
        Sequence_inst = class_quality(self.output_dir, os.path.normpath(self.primerfile), os.path.normpath(self.sample_sheet), self.minCount, self.alpha, self.ploidy, self.unique_prefix_name, self.os_param, self.allelelist, self.textbox_bool)
        #Sequence_inst.print_input_paths()
        self.textbox.insert('end-1c', " \n")
        
        #global is_working
        self.is_working = True
        
        allele_length_dict = {}
        print()
        self.textbox.insert('end-1c', "Get all Sample-Locus combo for each Genotype Length... \n")
        print()
        #self.textbox.insert('end-1c', " \n")
        if (self.ploidy == "diploid"):
            allele_length_dict, number2, number_filled2, number_empty2, samples_without_lengths, samples, duplicates = Sequence_inst.extract_AlleleLengths_diploid() 
        elif (self.ploidy == "haploid"):
            allele_length_dict, number2, number_filled2, number_empty2, samples_without_lengths, samples, duplicates = Sequence_inst.extract_AlleleLengths_haploid() #trimmed fastq files in /QC/
        print()
        self.textbox.insert('end-1c', "Extracting Allele Lengths and corresponding sequences finished! \n")
        self.textbox.insert('end-1c', "Results in the AllelesOut Folder which contains " + str(number_empty2) + " empty and " + str(number_filled2) + "non-empty files!  \n")
        self.textbox.insert('end-1c', "Number of Samples with 0 AlleleLengths: " + str(len(samples_without_lengths)) + "  \n\n\n\n\n")
        self.report_dict["AllelesOut"].append(number2)
        self.report_dict["AllelesOut"].append(number_filled2)
        self.report_dict["AllelesOut"].append(number_empty2)
        self.report_dict["AllelesOut"].append(samples_without_lengths)
        print()
        for sample_locus, length in allele_length_dict.items():
            for l, count in length.items():
                self.textbox.insert('end-1c', sample_locus + ":  \n")
                self.textbox.insert('end-1c', "length : " + str(l) + ", count : "  + str(count) + " \n")
        
        self.textbox.insert('end-1c', " \n")
        self.textbox.insert('end-1c', "Making ConsensusSequences... \n")
        number, number_filled, number_empty = Sequence_inst.RunConsensusAll()
        self.textbox.insert('end-1c', "Making ConsensusSequences finished! \n")
        self.textbox.insert('end-1c', "Results in the ConsensusOut Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + "non-empty files!  \n\n\n\n\n")
        self.report_dict["ConsensusOut"].append(number)
        self.report_dict["ConsensusOut"].append(number_filled)
        self.report_dict["ConsensusOut"].append(number_empty)
        print()
        
        number, number_filled, number_empty = Sequence_inst.joinSamplesSameMarker()
        self.textbox.insert('end-1c', "ConsensusTogether process finished! \n")
        self.textbox.insert('end-1c', "Results in the ConsensuTogether Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + "non-empty files!  \n\n\n\n\n")
        self.report_dict["ConsensusTogether"].append(number)
        self.report_dict["ConsensusTogether"].append(number_filled)
        self.report_dict["ConsensusTogether"].append(number_empty)
        print()
        self.textbox.insert('end-1c', " \n")
        self.textbox.insert('end-1c', "Correcting Sequences... \n")
        #Sequence_inst.correctAllSeqs()
        number, number_filled, number_empty = Sequence_inst.correctAllSeqs()
        self.textbox.insert('end-1c', "Correction process finished! \n")
        self.textbox.insert('end-1c', "Results in the Corrected Folder which contains " + str(number_empty) + " empty and " + str(number_filled) + "non-empty files!  \n\n\n\n\n")
        self.report_dict["Corrected"].append(number)
        self.report_dict["Corrected"].append(number_filled)
        self.report_dict["Corrected"].append(number_empty)
        print()
        
        self.textbox.insert('end-1c', " \n")
        self.textbox.insert('end-1c', "Calling Alleles... \n")
        additional_loci, markers_more_than2 = Sequence_inst.AlleleCall()
        self.textbox.insert('end-1c', "AlleleCall process finished! \n")
        self.textbox.insert('end-1c', "Results in the AlleleCall Folder!  \n\n\n\n\n")
        self.report_dict["AlleleCall"].append(number)
        self.report_dict["AlleleCall"].append(number_filled)
        self.report_dict["AlleleCall"].append(number_empty)
        if (self.allelelist == "" or os.path.basename(os.path.normpath(self.allelelist)) == "." or not(os.path.exists(self.allelelist))):
            self.textbox.insert('end-1c', "No AlleleList, Denovo Call!\n")
            self.textbox.insert('end-1c', "\n")
        else:
            self.textbox.insert('end-1c', "ALleleList is there, ListBased Call!\n")
        
        allelematrix = os.path.normpath(self.output_dir + "/AlleleCall//matrix.txt")
        if (self.os_param == "windows"):
            os.startfile(allelematrix)
        elif (self.os_param == "linux"):
            subprocess.call(["xdg-open", allelematrix])
        
        #if (self.textbox_bool):
            #self.textbox.insert('end-1c', "Output written into Genalex Format \n")
        #else:
            #self.textbox.insert('end-1c', "Output NOT written into Genalex Format \n")
        
        self.is_working = False


if __name__ == "__main__":
    main()