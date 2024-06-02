# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 16:53:48 2022

@author: Sebastian
"""
import sqlite3 as sq
import os

def main():
    
    print("For Testing.")
    
    workspace = os.getcwd()
    database_location = os.path.normpath(workspace + "/database.db")
    
    db_instance = SQdatabase(database_location) #local
    
    #Testing:
    row = ("Folder", "Project6", "Organism0", "Country0", "Locality0", "Sample0", "loci0", "AlleleIdx0", "Length0", "ABABABABA0", "Default0")
    db_instance.insert_record(row)
    
    multiple_rows = []
    #Generate more rows:
    for i in range(0, 5, 1):
        row = ("Folder" + str(i), "Project" + str(i), "Organism" + str(i), "Country" + str(i), "Locality" + str(i), "Sample" + str(i), "loci" + str(i), "AlleleIdx" + str(i), "Length" + str(i), "ABABABABA" + str(i), "Default" + str(i))
        multiple_rows.append(row)
    
    #insert all rows:
    for row in multiple_rows:
        db_instance.insert_record(row)
    
    #Show all records:
    print()
    print("Show all records:")
    rows = db_instance.get_all_records()
    for row in rows:
        print(row)
        
    print()
    count = db_instance.get_count_sequence("ABABABABA0")
    print(count)
    
    db_instance.closing()
    
    
    #root = root_window() 
    
    #root.mainloop()
    
    return 0

class SQdatabase: #With SQdatabase class we can create objects/Instanzes of SQdatabase class (connecting to the same database)
    def __init__(self, my_db):
        self.connection = sq.connect(my_db)
        self.cursor = self.connection.cursor() #NOT NULL PRIMARY KEY
        self.cursor.execute("""CREATE TABLE IF NOT EXISTS ssr_table(
                    OutputFolderName VARCHAR NOT NULL,
                    ProjectFolder VARCHAR NOT NULL,
                    Project VARCHAR NOT NULL,
                    Organism VARCHAR NOT NULL,
                    Country VARCHAR NOT NULL,
                    Locality VARCHAR NOT NULL,
                    Sample VARCHAR NOT NULL,
                    Loci VARCHAR NOT NULL,
                    AlleleIdx VARCHAR NOT NULL,
                    Length INTEGER NOT NULL,
                    AlleleSequence VARCHAR[MAX] NOT NULL,
                    UniqueID VARCHAR NOT NULL,
                    CONSTRAINT PK_ssr_table PRIMARY KEY(ProjectFolder, Project, Organism, Country, Locality, Sample, Loci, AlleleIdx, Length, AlleleSequence, UniqueID)
                    );""")
        self.connection.commit()
        #PRIMARY KEY(SampleName, LociName, AlleleIdx, AlleleName)
        
    def closing(self):
        print("Closing the database.")
        self.connection.close()
        
    def insert_record(self, row):
        query = """INSERT OR IGNORE INTO ssr_table
                   (OutputFolderName, ProjectFolder, Project, Organism, Country, Locality, Sample, Loci, AlleleIdx, Length, AlleleSequence, UniqueID) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
        """
        self.cursor.execute(query, row)
        #self.cursor.execute("""INSERT INTO my_table (FolderName, Project, Organism, Country, Locality, Sample, Loci, AlleleIdx, Length, AlleleSequence, UniqueID) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) ;""", row)
        self.connection.commit()
    
    def get_all_records(self):
        query = """SELECT * 
                   FROM ssr_table
        """
        self.cursor.execute(query)
        rows = self.cursor.fetchall()
        return rows
        
    def get_count_sequence(self, sequence):
        query = """SELECT COUNT(*)
                   FROM ssr_table
                   WHERE AlleleSequence = ?"
        """
        self.cursor.execute("SELECT COUNT(*) FROM ssr_table WHERE AlleleSequence = ?", [sequence])
        number = self.cursor.fetchone()
        return number
    
    def get_count_unique_sequence(self):
        query = """SELECT COUNT(DISTINCT AlleleSequence)
                   FROM ssr_table
                   WHERE AlleleSequence != ''
        """
        self.cursor.execute(query)
        number = self.cursor.fetchone()
        return number[0]
    
    
        
        
    def get_notempty_sequence_rows(self):
        query = """SELECT * 
                   FROM ssr_table
                   WHERE AlleleSequence != ''
        """
        self.cursor.execute(query)
        rows = self.cursor.fetchall()
        return rows
        
    def get_size_of_table(self):
        query = """SELECT COUNT(*) 
                   FROM ssr_table
        """
        self.cursor.execute(query)  # Executing the query
        number = self.cursor.fetchone()  # Fetching the result
        return number[0]  # Returning the count
    
    def get_notempty_sequence_rows_number(self):
        query = """SELECT COUNT(*) 
                   FROM ssr_table
                   WHERE AlleleSequence != ''
        """
        self.cursor.execute(query)  # Executing the query
        number = self.cursor.fetchone()  # Fetching the result
        return number[0]  # Returning the count
        
        
    def get_specified_subset(self, query, parameters):
        self.cursor.execute(query, parameters)
        rows = self.cursor.fetchall()
        return rows
        
    
    def get_unique_organism(self):
        query = """SELECT DISTINCT Organism 
                   FROM ssr_table
        """
        self.cursor.execute(query)
        rows = self.cursor.fetchall()
        return rows
    
    def get_unique_loci(self):
        query = """SELECT DISTINCT Loci 
                   FROM ssr_table
        """
        self.cursor.execute(query)
        rows = self.cursor.fetchall()
        return rows

    def get_unique_samples(self):
        query = """SELECT DISTINCT Sample 
                   FROM ssr_table
        """
        self.cursor.execute(query)
        rows = self.cursor.fetchall()
        return rows
    
    def get_unique_sequences(self):
        query = """SELECT DISTINCT AlleleSequence
                   FROM ssr_table
        """
        self.cursor.execute(query)
        rows = self.cursor.fetchall()
        return rows
    
    
    
    def deleting_data_from_folder(self, foldername):
        self.cursor.execute("DELETE from ssr_table WHERE OutputFolderName = ?", [foldername])
        self.connection.commit()
    
    
    
    
    
    
    
    
    
    
    
    
        
    def allelename_unique(self, allelename):
        self.cursor.execute("""SELECT AlleleSequence = ?, COUNT(*) from my_table GROUP BY AlleleSequence""", [allelename])
        self.connection.commit()
        
    def deleting(self):
        self.cursor.execute("""DELETE from my_table""")
        self.connection.commit()
        
    def delete_data_from_folder(self, foldername):
        self.cursor.execute("DELETE from my_table WHERE FolderName = ?", [foldername])
        self.connection.commit()
        
    def loci_number(self):
        self.cursor.execute("SELECT COUNT(DISTINCT LociName) from my_table GROUP BY LociName")
        self.connection.commit()
        
    def insert_row(self, row):
        self.cursor.execute("INSERT INTO my_table (FolderName, Project, Organism, Country, Locality, SampleName, LociName, AlleleIdx, Length, AlleleSequence) values(?, ?, ?, ?, ?, ?, ?, ?, ?, ?) ;", row)
        self.connection.commit()
    
    def get_size(self):
        return self.cursor.execute("SELECT COUNT(*) FROM my_table")
    
    def insert_sample(self, sample):
        self.cursor.execute("INSERT INTO my_table (SampleName) values(?);", [sample])
        self.connection.commit()
        
    def insert_loci(self, loci):
        self.cursor.execute("INSERT INTO my_table (LociName) values(?);", [loci])  
        self.connection.commit()
        
    def insert_alleleIdx(self, alleleIdx):
        self.cursor.execute("INSERT INTO my_table (AlleleIdx) values(?);", [alleleIdx])  
        self.connection.commit()
        
    def insert_alleleName(self, allelename):
        self.cursor.execute("INSERT INTO my_table (AlleleSequence) values(?);", [allelename]) 
        self.connection.commit()
    
    
    def insert_row_if_sample_does_not_exist(self, row, sample):
        self.cursor.execute("INSERT INTO my_table (SampleName, LociName, AlleleIdx, AlleleSequence) values(?, ?, ?, ?) IF NOT EXISTS SampleName = ? ;", (row, sample))
        self.connection.commit()
        
    def insert_many_rows(self, multiple_rows):
        self.cursor.executemany('''INSERT INTO my_table (SampleName, LociName, AlleleIdx, AlleleSequence) values(?, ?, ?, ?);''', multiple_rows)
        self.connection.commit()
        
    def update_sample_at_loci(self, sample,  loci):
        self.cursor.execute("UPDATE my_table SET Samplename = ? WHERE LociName = ?", (sample, loci))
        self.connection.commit()
        
    def update_loci_at_sample(self, loci,  sample):
        self.cursor.execute("UPDATE my_table SET LociName = ? WHERE SampleName = ?", (loci, sample))
        self.connection.commit()
        
    def update_alleleIdx_at_sample_loci_pair(self, alleleIdx, sample,  loci):
        self.cursor.execute("UPDATE my_table SET alleleIdx = ? WHERE SampleName = ? AND LociName = ?", (alleleIdx, sample, loci))
        self.connection.commit()
        
    def delete_row_at_sample_loci_pair(self, sample,  loci):
        self.cursor.execute("DELETE FROM my_table WHERE SampleName = ? AND LociName = ?", (sample, loci))
        self.connection.commit()
        
            
    def get_notempty_records(self):
        self.cursor.execute("SELECT Project, Organism, Country, Locality, SampleName, LociName, AlleleIdx, Length, AlleleSequence FROM my_table WHERE SampleName IS NOT NULL and LociName IS NOT NULL and AlleleIdx IS NOT NULL and Length IS NOT NULL and AlleleSequence IS NOT NULL and AlleleSequence != '' ")
        rows = self.cursor.fetchall()
        return rows
    
    def get_records_without_indices(self):
        self.cursor.execute("SELECT Project, Organism, Country, Locality, SampleName, LociName, AlleleSequence FROM my_table WHERE SampleName IS NOT NULL and LociName IS NOT NULL and AlleleIdx IS NOT NULL and AlleleName IS NOT NULL")
        rows = self.cursor.fetchall()
        return rows
    
    def get_records_for_loci(self, loci):
        self.cursor.execute("SELECT DISTINCT Project, Organism, Country, Locality, SampleName, AlleleIdx, Length, AlleleName FROM my_table WHERE AlleleName != '' and LociName = ? GROUP BY AlleleIdx", [loci])
        rows = self.cursor.fetchall()
        return rows
    
    def get_records_for_sample(self, sample):
        self.cursor.execute("SELECT Project, Organism, Country, Locality, LociName, AlleleIdx, Length, AlleleName FROM my_table WHERE AlleleName != '' and SampleName = ?", [sample])
        rows = self.cursor.fetchall()
        return rows
    
    def get_record_for_sample_loci(self, sample, loci):
        self.cursor.execute("SELECT Project, Organism, Country, Locality, SampleName, LociName, AlleleIdx, Length, AlleleName FROM my_table WHERE SampleName = ? AND LociName = ?", (sample, loci))
        rows = self.cursor.fetchall()
        return rows
    
    def get_alleleIdx_for1loci(self, loci):
        self.cursor.execute("SELECT AlleleIdx FROM my_table WHERE LociName = ?", [loci])
        rows = self.cursor.fetchall()
        return rows
    

    
    def get_alleleSequence_for1loci(self, loci):
        self.cursor.execute("SELECT DISTINCT AlleleName FROM my_table WHERE LociName = ?", [loci])
        rows = self.cursor.fetchall()
        return rows
    
    def get_distinct_alleleSequence_for1loci(self, loci):
        self.cursor.execute("SELECT AlleleName FROM my_table WHERE LociName = ?", [loci])
        rows = self.cursor.fetchall()
        return rows
    
    def get_alleleSequence_for1loci2(self, loci):
        self.cursor.execute("SELECT * FROM my_table WHERE LociName = ?", [loci])
        rows = self.cursor.fetchall()
        return rows
        
    def show_sample_names(self):
        #self.cursor.execute("SELECT DISTINCT SampleName FROM my_table")
        #rows = self.cursor.fetchall()
        #for row in rows:
            #print(row)
        for element in self.cursor.execute("SELECT DISTINCT SampleName FROM my_table"):
            print(element[0])
            
            
if __name__ == "__main__":
    main()
            
    