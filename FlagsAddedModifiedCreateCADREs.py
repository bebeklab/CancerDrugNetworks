# CADRE
# This tool is realeased under GNU Lesser General Public License 
# by Gurkan Bebek, Maya Monroe, Douglas Brubaker

# the following will take input files from the user 
# it will turn them into the txt files with lists of dysregulated pathways for each drug
# it names all output files based on the input files 
# to run type name of script then the arguments and their associated files
# if handing more than one file to an argument separate the file names by a comma (no spaces)
# if using multiple files in an argument, the order of the files must correspond between all arguments with multiple files
# the following uses tTest as the metric instead of signal2noise, none of the other versions except in trial2 and trial3 are like this!!!

import argparse
import os
import sys
import pandas as pd
from pandas import ExcelWriter
import numpy as np
from compiler.ast import flatten

#creates the arguments as flags
parser = argparse.ArgumentParser()
parser.add_argument("-rna", help='The RNA file name', action='append') #append lets the user designate more than one file, it adds the files designated together as a list
parser.add_argument("-cnv", help='The CNV file name', action='append')
parser.add_argument("-ic50", help='The IC-50 file name', action='append')
parser.add_argument("-jar", help='The full path to gsea jar file')
parser.add_argument("-kegg", help='The full path to the downloaded KEGG file')
parser.add_argument("-chip", help='The full path to the downloaded Gene chip file')
parser.add_argument("-output", help='The full path to the folder to save the gsea results in', default=os.getcwd()) #this sets the folder to save to's default as the current working directory if none is classified 
parser.add_argument("-pathwaycutoff", help='Apply pathway threshold of 3, true or false', default='TRUE', choices=["TRUE", "FALSE"]) #allows the user to decide if they want to apply the pathway cutoff or not
parser.add_argument("-featurecutoff", help='Apply the feature threshold, true or false', default='TRUE', choices=["TRUE", "FALSE"]) #allows the user to decide if they want to apply the feature cutoff or not

args = parser.parse_args() #this saves the argument names in a dictionary args, with the user designated files under the keys


gsea_jar_filepath = args.jar
kegg_download_filepath = args.kegg
gene_chip_download_filepath = args.chip
folder_to_save_to_filepath = args.output
pathwaycutoff = args.pathwaycutoff
featurecutoff = args.featurecutoff
#this sets the user input file names stored in arg as variables for each argument that only accepts one file as an input 

list_significant_pathways = {} #this creates a dictionary to store the pathways for each drug from both databases in, we can't use a dataframe because the columns are different lengths
list_significant_genes_rna = {} #this creates a dictionary to store the genes for each drug from both databases in for the rna data
list_significant_genes_cnv = {} #this creates a dictionary to store the genes for each drug from both databases in for the cnv data

list_rna_basenames = []
list_cnv_basenames = []
list_ic50_basenames = []
#this creates empty lists to store the basenames of the files for the different databases 

num_databases = len(str(args.rna).strip('[]').strip("''").split(",")) #this lets us keep track of the number of databases being combined 

iterate = 0
for rna_filename in str(args.rna).strip('[]').strip("''").split(","): #this runs through all the rna files listed (need to split at the comma and strip the brackets and quotations, the input must be a string in order to do this)
	rna_filename = os.path.basename(rna_filename)#this records the entire file name for the rna data file
	rna_basename = rna_filename.split(".")[0] #this strips of the .txt part
	list_rna_basenames.append(rna_basename) #this stores the basename in the list of rna basenames for all databases
	cnv_filename = os.path.basename(str(args.cnv).strip('[]').strip("''").split(",")[iterate]) #this records the entire file name for the cnv data, the count enables it to iterate through the listed cnv files as we iterate through the rna file names
	cnv_basename = cnv_filename.split(".")[0] #this strips the filename of the .txt part
	list_cnv_basenames.append(cnv_basename) #this stores the basename in the list of cnv basenames for all databases
	ic50_filename = os.path.basename(str(args.ic50).strip('[]').strip("''").split(",")[iterate]) #this records the entire file name for the ic50 data, the count enables it to iterate through the listed ic50 files as we iterate through the rna file names
	ic50_basename = ic50_filename.split(".")[0] #this strips the filename of the .txt part
	list_ic50_basenames.append(ic50_basename) #this stores the basename in the list of ic50 basenames for all databases
	iterate = iterate + 1 #this adds one to the count every time we go through another rna filename, which lets us iterate through the correct cnv and ic50 files


	#now that we have the names of the user input files saved as variables, we can run them through the rest of programs



	# The following pulls out cell lines that aren't present in both RNA data and in cnv data 

	with open("%s" % (rna_filename),"U") as f:
		RNA_data = []
		for line in f:
			spl = line.split('\t')
			RNA_data.append(spl)
	f.close()
	#this converts the RNA data into a python array
	#note that cell lines are listed in the first row

	RNA_cell_lines = RNA_data[0] # this stores all RNA cell line names in a list 


	with open("%s" % (cnv_filename),"U") as f:
		cnv_data = []
		for line in f:
			spl = line.split()
			cnv_data.append(spl)
	f.close()
	#this converts the IC-50 data into a python array
	#note that cell lines are listed in the first row

	cnv_cell_lines = cnv_data[0]
	# this stores the cnv cell line names in a list

	cnv_cell_lines_to_save = [] #this creates an empty list where we can store the indexes of cnv cell lines that are present in the rna data
	for index, name in enumerate(cnv_cell_lines): #goes through each cell line name
		name_string = ''.join(name) #converts the cell line name into a string
		if name_string in RNA_cell_lines: #checks to see if the cnv cell line name is in the list of rna cell lines
			cnv_cell_lines_to_save.append(index) #if it is then the index of the cell line is saved in a list
	# This stores the indexes of the cnv cell lines that are present in the RNA data

	cnv_data_filtered_with_rna = [] #this creates an array where we can store all the cnv data for cell lines that are present in the rna data 
	count = 0
	i = len(cnv_data)
	while (i > count): #this iterates through the rows of the cnv data
		new_row = []
		old_row = cnv_data[count]
		for index in cnv_cell_lines_to_save: #this iterates through the cnv cell line columns of each row that are present in the rna data  
			new_row.append(str(old_row[index])) #this adds the columns present in rna data onto the new row of cnv data
		cnv_data_filtered_with_rna.append(new_row) #this adds each newly filtered row of cnv data onto the array of newly filtered cnv data 
		count = count + 1 
	
	RNA_cell_lines_to_save = [] #this creates an empty list where we can store the indexes of rna cell lines that are present in the cnv data
	for index, name in enumerate(RNA_cell_lines): #goes through each cell line name
		name_string = ''.join(name) #converts the cell line name into a string
		if name_string in cnv_cell_lines: #checks to see if the rna cell line name is in the list of cnv cell lines
			RNA_cell_lines_to_save.append(index) #if it is then the index of the cell line is saved in a list
	#this stores the indexes of the RNA cell lines that are present in the IC-50 data

	RNA_data_filtered_with_cnv = [] #this creates an array where we can store all the rna data for cell lines that are present in the cnv data 
	count = 0
	i = len(RNA_data)
	while (i > count): #this iterates through the rows of rna data
		new_row = []
		old_row = RNA_data[count]
		for index in RNA_cell_lines_to_save: #this iterates through the rna cell line columns of each row that are present in the cnv data  
			new_row.append(str(old_row[index])) #this adds the columns present in cnv data onto the new row of rna data
		RNA_data_filtered_with_cnv.append(new_row) #this adds each newly filtered row of rna data onto the array of newly filtered rna data 
		count = count + 1
	


	# The following pulls out cell lines that aren't present in both RNA data/CNV data and in IC-50 data 

	RNA_data = RNA_data_filtered_with_cnv #this renames the RNA filtered array created earlier

	cnv_data = cnv_data_filtered_with_rna #this renames the cnv filtered array created earlier

	RNA_cell_lines = RNA_data[0] # this stores all RNA cell line names in a list 


	with open("%s" % (ic50_filename),"U") as f:
		IC50_data = []
		for line in f:
			spl = line.split()
			IC50_data.append(spl)
	f.close()
	#this converts the IC-50 data into a python array
	#note that cell lines are listed in the first column

	IC50_cell_lines = [] #creates a list to save the cell line names to
	count = 0 
	i = len(IC50_data)
	while (i > count): #this iterates through the rows of the IC50 data
		row = IC50_data[count]
		name = row[0] #this pulls out the first column of each row (which contains the cell line names)
		IC50_cell_lines.append(name) #saves the cell line names in a list
		count = count + 1
	# this stores the IC-50 cell line names in a list

	IC50_cell_lines_to_save = [] #this creates an empty list where we can store the indexes of ic50 cell lines that are present in the rna and cnv data
	for index, name in enumerate(IC50_cell_lines): #goes through each cell line name
		name_string = ''.join(name) #converts the name into a string 
		if name_string in RNA_cell_lines: #checks to see if the ic50 cell line name is in the list of rna cell lines
			IC50_cell_lines_to_save.append(index) #saves the index of the cell line if it is
	# This stores the indexes of the IC-50 cell lines that are present in the RNA data (and thus the cnv data by default)

	IC50_data_filtered = [] #this creates an array to store the filtered ic50 data in
	for index in IC50_cell_lines_to_save: #goes through the indexes of columns that are present in both ic50 data and the rna/cnv data 
		row = IC50_data[index] #this records the data of the rows 
		row_string = '\t'.join(str(v) for v in row) #turns the columns of each row into a string and then adds a tab in between columns 
		IC50_data_filtered.append(row_string) #this adds the rows with tabs to the array
	#this creates an array that contains the IC50 values of only the cell lines present in RNA data 

	RNA_cnv_cell_lines_to_save = []
	for index, name in enumerate(RNA_cell_lines):
		name_string = ''.join(name)
		if name_string in IC50_cell_lines:
			RNA_cnv_cell_lines_to_save.append(index)
	#this stores the indexes of the RNA/cnv cell lines that are present in the IC-50 data

	RNA_data_filtered = []
	count = 0
	i = len(RNA_data)
	while (i > count):
		new_row = []
		old_row = RNA_data[count]
		for index in RNA_cnv_cell_lines_to_save:
			new_row.append(old_row[index])
		new_row_string = '\t'.join(str(v) for v in new_row)
		RNA_data_filtered.append(new_row_string)
		count = count + 1
	#this stores the rna data of the cell lines for which we have ic50 data

	cnv_data_filtered = []
	count = 0
	i = len(cnv_data)
	while (i > count):
		new_row = []
		old_row = cnv_data[count]
		for index in RNA_cnv_cell_lines_to_save:
			new_row.append(old_row[index])
		new_row_string = '\t'.join(str(v) for v in new_row)
		cnv_data_filtered.append(new_row_string)
		count = count + 1
	#this stores the cnv data of the cell lines for which we have ic50 data
	
	with open("%s_filtered.txt" % (ic50_basename),"w") as f:
		f.write('\n'.join(IC50_data_filtered) + '\n')
		f.close( )
	
	with open("%s_filtered.txt" % (rna_basename),"w") as f:
		f.write('\n'.join(RNA_data_filtered) + '\n')
		f.close( )

	with open("%s_filtered.txt" % (cnv_basename),"w") as f:
		f.write('\n'.join(cnv_data_filtered) + '\n')
		f.close( )
	#this saves the filtered cnv rna and ic50 files to text documents




	#This script converts the IC-50 CGP data into a .cls file for each drug for the RNA and cnv data
	#it is an elongated version of the RemoveCellLinesWithNoDataForEachCGPDrugAndClassify.py
	#it removes cell lines for which we don't have IC-50 data
	#it classifies each cell line as sensitive or resistant 
	#it makes the corresponding RNA and cnv expression data files for each drug 

	IC50_data = pd.read_csv("%s_filtered.txt" % (ic50_basename), delimiter='\t')
	#This stores the IC50_data in pandas dataframe

	cell_lines = IC50_data.Cell_Line
	#This stores the Cell lines in their own array

	file = ExcelWriter('%s_FILTERED_AND_CLASSIFIED.xlsx' %(ic50_basename)) #this creates the excel file we want to write to

	def definesensitivity(x):
		if x < 1:
			return 'S'
		else:
			return 'R'
	# this writes the function we'll use later to classify each cell line as sensitive or resistant

	drugs_with_enough_variation = pd.DataFrame() #this creates a dataframe that we will fill with keys that represent the drugs with enough variation in their drug response

	for drug_name in IC50_data.keys():
		if drug_name != 'Cell_Line': # this keeps us from making a Cell line sheet
			drug_data = pd.DataFrame()
			drug_data['Cell_Line'] = IC50_data.loc[:,('Cell_Line')] #this adds on the list of cell lines to the drug specific dataframe
			drug_data[drug_name] = IC50_data.loc[:,(drug_name)] #this adds on the data for each drug to the dataframe
			drug_data_filtered = drug_data.dropna() #this drops all cell lines for which we don't have data 
		
			drug_data_classified = pd.DataFrame() #this creates a dataframe that we can store the sensitivity classification in 
			drug_data_classified['Sensitivity'] = drug_data_filtered.loc[:,(drug_name)].map(definesensitivity) #this calls our sensitivity function and applies it to each drug
			#this applies the sensitivity function to the IC-50 data of each drug 
		
			drug_data_tofile = pd.DataFrame() #this creates a dataframe that will print correctly
			drug_data_tofile['Cell_Line'] = drug_data_filtered.loc[:,('Cell_Line')]
			drug_data_tofile[drug_name] = drug_data_filtered.loc[:,(drug_name)]
			drug_data_tofile['Sensitivity'] = drug_data_classified.loc[:,('Sensitivity')]
			#this copies over all the information to the file we will print
		
			drug_data_tofile.to_excel(file,'%s' % (drug_name)) #this creates a sheet in the file for each drug
		
			sensitive_lines = drug_data_tofile[drug_data_tofile['Sensitivity'] == 'S'] #pulls out the sensitive lines from the dataframe
			num_sensitive = len(sensitive_lines) #identifies the number of sensitive lines
			resistant_lines = drug_data_tofile[drug_data_tofile['Sensitivity'] == 'R'] #pulls out the resistant lines
			num_resistant = len(resistant_lines) #identifies the number of sensitive lines
		
			if num_resistant >= 3 and num_sensitive >= 3: #this applies the variation cutoff required by gsea
				cls_file_data = np.array(drug_data_tofile.Sensitivity) #creates a numpy array from the sensitivity info in the drug_data_tofile
				cls_file_data_condensed = [' '.join(cls_file_data)] #converts a list of characters into a string
				cls_file = [] #creates an array for the cls file data
				num_cell_lines = len(cls_file_data) #we need to include the number of cell lines in the cls file at the top of the cls file
				cls_file_first_line = [num_cell_lines, 2, 1] #writes the first line of the cls file as a list
				cls_file_first_line_condensed = (str(cls_file_first_line).replace(',','')).strip('[]') #this removes the commas and unnecessary brackets converts the integers to strings so that they can be joined
				cls_file.append(cls_file_first_line_condensed) #this appends the first line to the cls file array
				cls_file.append(['# Sensitive Resistant']) #this appends the second line to the cls file array
				cls_file.append(cls_file_data_condensed) #this appends the third line (sensitivity classification) to the cls file array
				cls_file_flatten = list(flatten(cls_file)) #condenses the lists into a single list to get rid of the quotation marks and brackets
				np.savetxt('%s_%s.cls' %(ic50_basename, drug_name), cls_file_flatten, fmt='%s') #saves the cls file array as a cls file, indicates that the data is to be saved as strings
				drugs_with_enough_variation[drug_name] = "" #creates a list of the drugs that pass the variation requirements
				#the above creates .cls files for only the drugs that have enough variation
			
				with open("%s_filtered.txt" %(rna_basename),"U") as f:
					RNA_data = []
					for line in f:
						spl = line.split('\t')
						RNA_data.append(spl)
				f.close()
				#this calls the RNA data
			
				RNA_cnv_cell_lines_to_save = drug_data_tofile.index.tolist() #this converts the list of dataframe indexes from the dataframe of ic50 data without the cell lines for which we have no data into a list
			
				RNA_data_drug = []
				count = 0 
				i = len(RNA_data)
				while (i > count):
					new_row = []
					old_row = RNA_data[count]
					new_row.append(old_row[0]) #this adds on the gene names, which are located in the first column but that aren't a cell line contained in the IC-50 data
					for index in RNA_cnv_cell_lines_to_save:
						new_row.append(old_row[(index + 1)]) #because the indexes of the IC-50 dataframe start with the data, not the keys where cell_line and such are contained, we must add one to the index gotten from the data frame in order to correct the shift that occurs, since the indexes of the rows of RNA data include the column that contains gene names
					new_row_string = '\t'.join(str(v) for v in new_row)
					RNA_data_drug.append(new_row_string)
					count = count + 1
				#this removes the cell lines that don't have IC-50 data for a particular drug
			
				with open("%s_%s.txt" %(rna_basename, drug_name),"w") as f:
					f.write('\n'.join(RNA_data_drug) + '\n')
					f.close( )
				#this creates an RNA file for each drug, but only for the drugs that have enough variation
			
			
				with open("%s_filtered.txt" %(cnv_basename),"U") as f:
					cnv_data = []
					for line in f:
						spl = line.split('\t')
						cnv_data.append(spl)
				f.close()
				#this calls the RNA data
			
				cnv_data_drug = []
				count = 0 
				i = len(cnv_data)
				while (i > count):
					new_row = []
					old_row = cnv_data[count]
					new_row.append(old_row[0]) #this adds on the gene names, which are located in the first column but that aren't a cell line contained in the IC-50 data
					for index in RNA_cnv_cell_lines_to_save:
						new_row.append(old_row[(index + 1)]) #because the indexes of the IC-50 dataframe start with the data, not the keys where cell_line and such are contained, we must add one to the index gotten from the data frame in order to correct the shift that occurs, since the indexes of the rows of RNA data include the column that contains gene names
					new_row_string = '\t'.join(str(v) for v in new_row)
					cnv_data_drug.append(new_row_string)
					count = count + 1
				#this removes the cell lines that don't have IC-50 data for a particular drug 
			
				with open("%s_%s.txt" %(cnv_basename, drug_name),"w") as f:
					f.write('\n'.join(cnv_data_drug) + '\n')
					f.close( )
				#this creates a cnv file for each drug, but only for the drugs that have enough variation
	
	drugs_with_enough_variation.to_excel('%s_Drugs_Run_Through_Gsea.xlsx' % (ic50_basename), index=False) #this stores the drugs with enough variation to be run by GSEA in a file for user reference




	#this script runs our RNA and IC-50 files for each drug through GSEA using java and saves the output files from gsea

	drugs_to_run = pd.read_excel('%s_Drugs_Run_Through_Gsea.xlsx' %(ic50_basename)) #creates the list of drugs to run
	gsea_folder_names_rna = [] #creates the array to store the folder names in

	for drug_name in drugs_to_run.keys(): #this iterates through the list of drugs 
		RNA_expression_file = os.path.abspath('%s_%s.txt' %(rna_basename, drug_name))
		cls_file = os.path.abspath('%s_%s.cls' %(ic50_basename, drug_name))
		#this gets the absolute paths for the RNA and IC50 files, which the javacode requires to work and which will vary from person to person depending on what folder they've saved in etc..
		output_file_name = 'RNA_%s' %(drug_name)

		java_script = "java -cp %s -Xmx1024m xtools.gsea.Gsea -res %s -cls %s#Sensitive_versus_Resistant -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 2000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label %s -metric tTest -sort real -order descending -chip %s -include_only_symbols false -make_sets false -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 1500 -set_min 5 -zip_report false -out %s -gui false" %(gsea_jar_filepath, RNA_expression_file, cls_file, kegg_download_filepath, output_file_name, gene_chip_download_filepath, folder_to_save_to_filepath)
		#this changes the java script for each drug so that it has the appropriate file paths and output file for each drug
	
		java_command = os.system(java_script) #this stores info about the execution of the java script in the command terminal
	
		if java_command != 0: #so long as the command is executed properly, the command terminal returns 0
			print "Error when calling GSEA for %s files" %(drug_name)
		# if the command terminal doesn't call GSEA properly an error is displayed telling you so
	
		all_subdirs = [d for d in os.listdir('.') if os.path.isdir(d)]
		latest_subdir = max(all_subdirs, key=os.path.getmtime)
		#this gives us the most recent folder created in our current working directory
		gsea_folder_names_rna.append(latest_subdir)
		#this saves the names of these folders into a python array

	with open("gsea_folder_names_%s.txt" %(rna_basename),"w") as f:
		f.write('\n'.join(gsea_folder_names_rna) + '\n')
		f.close( )
	#this saves the folders as a text document for user reference 




	#this script runs our CNV and IC-50 files for each drug through GSEA using java and saves the output files from gsea

	drugs_to_run = pd.read_excel('%s_Drugs_Run_Through_Gsea.xlsx' %(ic50_basename)) #creates the list of drugs to run
	gsea_folder_names_cnv = [] #creates the array to store the folder names in

	for drug_name in drugs_to_run.keys(): #this iterates through the list of drugs 
		cnv_expression_file = os.path.abspath('%s_%s.txt' %(cnv_basename, drug_name))
		cls_file = os.path.abspath('%s_%s.cls' %(ic50_basename, drug_name))
		#this gets the absolute paths for the CNV and IC50 files, which the javacode requires to work and which will vary from person to person depending on what folder they've saved in etc..
		output_file_name = 'CNV_%s' %(drug_name)

		java_script = "java -cp %s -Xmx1024m xtools.gsea.Gsea -res %s -cls %s#Sensitive_versus_Resistant -gmx %s -collapse false -mode Max_probe -norm meandiv -nperm 2000 -permute gene_set -rnd_type no_balance -scoring_scheme weighted -rpt_label %s -metric tTest -sort real -order descending -chip %s -include_only_symbols false -make_sets false -median false -num 100 -plot_top_x 20 -rnd_seed timestamp -save_rnd_lists false -set_max 1500 -set_min 5 -zip_report false -out %s -gui false" %(gsea_jar_filepath, cnv_expression_file, cls_file, kegg_download_filepath, output_file_name, gene_chip_download_filepath, folder_to_save_to_filepath)
		#this changes the java script for each drug so that it has the appropriate file paths and output file for each drug
	
		java_command = os.system(java_script) #this stores info about the execution of the java script in the command terminal
	
		if java_command != 0: #so long as the command is executed properly, the command terminal returns 0
			print "Error when calling GSEA for %s files" %(drug_name)
		# if the command terminal doesn't call GSEA properly an error is displayed telling you so
	
		all_subdirs = [d for d in os.listdir('.') if os.path.isdir(d)]
		latest_subdir = max(all_subdirs, key=os.path.getmtime)
		#this gives us the most recent folder created in our current working directory
		gsea_folder_names_cnv.append(latest_subdir)
		#this saves the names of these folders into a python array

	with open("gsea_folder_names_%s.txt" %(cnv_basename),"w") as f:
		f.write('\n'.join(gsea_folder_names_cnv) + '\n')
		f.close( )
	#this saves the folder as a text document for user reference 



	#takes the excel files produced by GSEA for RNA expression data
	#uses both sensitive and resistant results
	#removes pathways where p-value>0.05 or FDR q-value > 0.05
	#classifies pathways as up regulated or down regulated (based on sensitive or resistant results)
	#takes the enrichment score files
	#calculates z-scores
	#generates the list of up and down differentially expressed genes

	drugs_run = pd.read_excel('%s_Drugs_Run_Through_Gsea.xlsx' % (ic50_basename)) #creates the list of drugs that were run

	for folder_name in gsea_folder_names_rna: #this iterates through the list of folders 
	
		for drug_name in drugs_run.keys(): #this iterates through the list of drug names
			if drug_name in folder_name: #this lets us identify the drug name of the folder and use the variable drug_name
	
				trial_num = folder_name.split('.')[2] #this retrieves the gsea trial number 

				#first we deal with the results for sensitive gene sets
				if os.path.isfile("%s/gsea_report_for_S_%s.xls" %(folder_name, trial_num)): #sometimes gsea returns reports for S, sometimes it returns reports for sensitive, this checks to see which it returned for this trial 
					sensitive_data = pd.read_csv("%s/gsea_report_for_S_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the gsea report from the drugs results folder (titled drug_name) created for each drug in the folder where the user saved them using a relative file path
				else:
					sensitive_data = pd.read_csv("%s/gsea_report_for_Sensitive_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the gsea report from the drugs results folder (titled drug_name) created for each drug in the folder where the user saved them using a relative file path
				significant_sensitive_data = sensitive_data[(sensitive_data['NOM p-val'] < 0.05) & (sensitive_data['FDR q-val'] < 0.05)] #this removes the rows that don't meet the significance requirements
				significant_sensitive_pathways = significant_sensitive_data['NAME'].tolist() #this converts the column of the pandas dataframe containing the names of the pathways into a list 
				#we come up with a list of down regulated pathways
	
				significant_classified_sensitive_pathways = []	#this creates a list to store the pathways classified as down regulated in 
				i = len(significant_sensitive_pathways)
				count = 0
				while i > count: #this iterates through the list of unclassified pathways
					pathway = significant_sensitive_pathways[count]
					classified_pathway = str(pathway).replace('KEGG','DN') #this replaces KEGG with DN
					significant_classified_sensitive_pathways.append(classified_pathway) #this appends the classified pathway to the downregulated list 
					count = count + 1
				#this replaces KEGG with DN for every pathway
	
	
	
				#then we deal with the results for resistant gene sets
				if os.path.isfile("%s/gsea_report_for_R_%s.xls" %(folder_name, trial_num)): #sometimes gsea returns reports for S, sometimes it returns reports for sensitive, this checks to see which it returned for this trial 
					resistant_data = pd.read_csv("%s/gsea_report_for_R_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the gsea report from the drugs results folder (titled drug_name) created for each drug in the folder where the user saved them using a relative file path
				else:
					resistant_data = pd.read_csv("%s/gsea_report_for_Resistant_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the gsea report from the drugs results folder (titled drug_name) created for each drug in the folder where the user saved them using a relative file path
				significant_resistant_data = resistant_data[(resistant_data['NOM p-val'] < 0.05) & (resistant_data['FDR q-val'] < 0.05)] #this should remove the rows that don't meet the significance requirements
				significant_resistant_pathways = significant_resistant_data['NAME'].tolist() #this converts the column of the pandas dataframe containing the names of the pathways into a list 
				#we come up with a list of up regulated pathways
	
				significant_classified_resistant_pathways = [] #this creates a list to store the pathways classified as upregulated in 
				i = len(significant_resistant_pathways)
				count = 0
				while i > count: #this iterates through the list of unclassified pathways
					pathway = significant_resistant_pathways[count]
					classified_pathway = str(pathway).replace('KEGG','UP') #this replaces KEGG with DN
					significant_classified_resistant_pathways.append(classified_pathway) #this appends the classified pathway to the upregulated list
					count = count + 1
				#this replaces KEGG with UP for every pathway
	
	
	
	
				#then we add the list of downregulated pathways to the list of upregulated pathways, creating one list of dysregulated pathways
				significant_pathways = significant_classified_sensitive_pathways + significant_classified_resistant_pathways
	
				list_significant_pathways[drug_name] = significant_pathways
				#this creates a list of pathways for each drug in the dataframe


				#now we generate the list of significantly enriched genes
				if os.path.isfile("%s/ranked_gene_list_S_versus_R_%s.xls" %(folder_name, trial_num)): #this checks to see if the list is S vs R
					ranked_gene_list = pd.read_csv("%s/ranked_gene_list_S_versus_R_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the ranked gene list from the gsea report folder
					average = ranked_gene_list['SCORE'].mean(axis=0) #this computes the average of the enrichment scores
					stdev = ranked_gene_list['SCORE'].std(axis=0) #this computes the stdev of the enrichment scores
					subtracted_scores = ranked_gene_list['SCORE'].sub(average) #this subtracts the average from the scores
					z_scores = subtracted_scores.divide(stdev) #this divides the subtracted scores by the stdev to give the z score 
					ranked_gene_list['Z_SCORES'] = z_scores #this adds the z-scores as their own column in the dataframe
				
					significant_downregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] > 2)] #this pulls out the genes who's z-scores are greater than 2
					names_significant_downregulated_data = significant_downregulated_data.loc[:,('NAME')]
					significant_downregulated_data.loc[:,('NAME')] = 'DN_' + names_significant_downregulated_data.astype(str) #this appends the prefix "DN" to the gene names
					significant_downregulated_genes = significant_downregulated_data['NAME'].tolist() #this pulls out the downregulated classified gene names and turns them into a list
					significant_upregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] < -2)] #this pulls out the genes who's z-scores are less than -2
					names_significant_upregulated_data = significant_upregulated_data.loc[:,('NAME')]
					significant_upregulated_data.loc[:,('NAME')] = 'UP_' + names_significant_upregulated_data.astype(str) #this appends the prefix "UP" to the gene names
					significant_upregulated_genes = significant_upregulated_data['NAME'].tolist() #this pulls out the upregulated classified genes and turns them into a list
				
				elif os.path.isfile("%s/ranked_gene_list_Sensitive_versus_Resistant_%s.xls" %(folder_name, trial_num)): #this takes into account the fact that gsea might have named the file Sensitive vs. Resistant instead of S vs R
					ranked_gene_list = pd.read_csv("%s/ranked_gene_list_Sensitive_versus_Resistant_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the ranked gene list from the gsea report folder
					average = ranked_gene_list['SCORE'].mean(axis=0) #this computes the average of the enrichment scores
					stdev = ranked_gene_list['SCORE'].std(axis=0) #this computes the stdev of the enrichment scores
					subtracted_scores = ranked_gene_list['SCORE'].sub(average) #this subtracts the average from the scores
					z_scores = subtracted_scores.divide(stdev) #this divides the subtracted scores by the stdev to give the z score 
					ranked_gene_list['Z_SCORES'] = z_scores #this adds the z-scores as their own column in the dataframe
				
					significant_downregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] > 2)] #this pulls out the genes who's z-scores are greater than 2
					names_significant_downregulated_data = significant_downregulated_data.loc[:,('NAME')]
					significant_downregulated_data.loc[:,('NAME')] = 'DN_' + names_significant_downregulated_data.astype(str) #this appends the prefix "DN" to the gene names
					significant_downregulated_genes = significant_downregulated_data['NAME'].tolist() #this pulls out the downregulated classified gene names and turns them into a list
					significant_upregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] < -2)] #this pulls out the genes who's z-scores are less than -2
					names_significant_upregulated_data = significant_upregulated_data.loc[:,('NAME')]
					significant_upregulated_data.loc[:,('NAME')] = 'UP_' + names_significant_upregulated_data.astype(str) #this appends the prefix "UP" to the gene names
					significant_upregulated_genes = significant_upregulated_data['NAME'].tolist() #this pulls out the upregulated classified genes and turns them into a list
				
				
				elif os.path.isfile("%s/ranked_gene_list_R_versus_S_%s.xls" %(folder_name, trial_num)):  #this classifies the genes if the enrichment report is for R vs S
					ranked_gene_list = pd.read_csv("%s/ranked_gene_list_R_versus_S_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the ranked gene list from the gsea report folder
					average = ranked_gene_list['SCORE'].mean(axis=0) #this computes the average of the enrichment scores
					stdev = ranked_gene_list['SCORE'].std(axis=0) #this computes the stdev of the enrichment scores
					subtracted_scores = ranked_gene_list['SCORE'].sub(average) #this subtracts the average from the scores
					z_scores = subtracted_scores.divide(stdev) #this divides the subtracted scores by the stdev to give the z score 
					ranked_gene_list['Z_SCORES'] = z_scores #this adds the z-scores as their own column in the dataframe
				
					significant_downregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] < -2)] #this pulls out the genes who's z-scores are less than -2
					names_significant_downregulated_data = significant_downregulated_data.loc[:,('NAME')]
					significant_downregulated_data.loc[:,('NAME')] = 'DN_' + names_significant_downregulated_data.astype(str) #this appends the prefix "DN" to the gene names
					significant_downregulated_genes = significant_downregulated_data['NAME'].tolist() #this pulls out the downregulated classified gene names and turns them into a list
					significant_upregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] > 2)] #this pulls out the genes who's z-scores are greater than 2
					names_significant_upregulated_data = significant_upregulated_data.loc[:,('NAME')]
					significant_upregulated_data.loc[:,('NAME')] = 'UP_' + names_significant_upregulated_data.astype(str) #this appends the prefix "UP" to the gene names
					significant_upregulated_genes = significant_upregulated_data['NAME'].tolist() #this pulls out the upregulated classified genes and turns them into a list
					
				else: #this takes into account the fact that gsea might have named it Resistant versus Sensitive instead of R versus S
					ranked_gene_list = pd.read_csv("%s/ranked_gene_list_Resistant_versus_Sensitive_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the ranked gene list from the gsea report folder
					average = ranked_gene_list['SCORE'].mean(axis=0) #this computes the average of the enrichment scores
					stdev = ranked_gene_list['SCORE'].std(axis=0) #this computes the stdev of the enrichment scores
					subtracted_scores = ranked_gene_list['SCORE'].sub(average) #this subtracts the average from the scores
					z_scores = subtracted_scores.divide(stdev) #this divides the subtracted scores by the stdev to give the z score 
					ranked_gene_list['Z_SCORES'] = z_scores #this adds the z-scores as their own column in the dataframe
				
					significant_downregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] < -2)] #this pulls out the genes who's z-scores are less than -2
					names_significant_downregulated_data = significant_downregulated_data.loc[:,('NAME')]
					significant_downregulated_data.loc[:,('NAME')] = 'DN_' + names_significant_downregulated_data.astype(str) #this appends the prefix "DN" to the gene names
					significant_downregulated_genes = significant_downregulated_data['NAME'].tolist() #this pulls out the downregulated classified gene names and turns them into a list
					significant_upregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] > 2)] #this pulls out the genes who's z-scores are greater than 2
					names_significant_upregulated_data = significant_upregulated_data.loc[:,('NAME')]
					significant_upregulated_data.loc[:,('NAME')] = 'UP_' + names_significant_upregulated_data.astype(str) #this appends the prefix "UP" to the gene names
					significant_upregulated_genes = significant_upregulated_data['NAME'].tolist() #this pulls out the upregulated classified genes and turns them into a list
				
			
				#add the two lists of genes together
				significant_dysregulated_genes = significant_downregulated_genes + significant_upregulated_genes
			
				#then we save the list of dysregulated genes in the dictionary for each drug
				list_significant_genes_rna[drug_name] = significant_dysregulated_genes





	#takes the excel files produced by GSEA for CNV expression data
	#takes the enrichment score files
	#calculates z-scores
	#generates the list of genes with up and down instances of CNVs

	drugs_run = pd.read_excel('%s_Drugs_Run_Through_Gsea.xlsx' % (ic50_basename)) #creates the list of drugs that were run

	for folder_name in gsea_folder_names_cnv: #this iterates through the list of folders 
	
		for drug_name in drugs_run.keys(): 
			if drug_name in folder_name: #this lets us identify the drug name of the folder and use the variable drug_name
	
				trial_num = folder_name.split('.')[2] #retrieves the trial number of the gsea run

			
				#now we generate the list of significantly enriched genes
				if os.path.isfile("%s/ranked_gene_list_S_versus_R_%s.xls" %(folder_name, trial_num)): #this checks to see if the list is S vs R
					ranked_gene_list = pd.read_csv("%s/ranked_gene_list_S_versus_R_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the ranked gene list from the gsea report folder
					average = ranked_gene_list['SCORE'].mean(axis=0) #this computes the average of the enrichment scores
					stdev = ranked_gene_list['SCORE'].std(axis=0) #this computes the stdev of the enrichment scores
					subtracted_scores = ranked_gene_list['SCORE'].sub(average) #this subtracts the average from the scores
					z_scores = subtracted_scores.divide(stdev) #this divides the subtracted scores by the stdev to give the z score 
					ranked_gene_list['Z_SCORES'] = z_scores #this adds the z-scores as their own column in the dataframe
				
					significant_downregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] > 2)] #this pulls out the genes who's z-scores are greater than 2
					names_significant_downregulated_data = significant_downregulated_data.loc[:,('NAME')]
					significant_downregulated_data.loc[:,('NAME')] = 'DN_' + names_significant_downregulated_data.astype(str) #this appends the prefix "DN" to the gene names
					significant_downregulated_genes = significant_downregulated_data['NAME'].tolist() #this pulls out the downregulated classified gene names and turns them into a list
					significant_upregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] < -2)] #this pulls out the genes who's z-scores are less than -2
					names_significant_upregulated_data = significant_upregulated_data.loc[:,('NAME')]
					significant_upregulated_data.loc[:,('NAME')] = 'UP_' + names_significant_upregulated_data.astype(str) #this appends the prefix "UP" to the gene names
					significant_upregulated_genes = significant_upregulated_data['NAME'].tolist() #this pulls out the upregulated classified genes and turns them into a list
				
				elif os.path.isfile("%s/ranked_gene_list_Sensitive_versus_Resistant_%s.xls" %(folder_name, trial_num)): #this takes into account the fact that gsea might have named the file Sensitive vs. Resistant instead of S vs R
					ranked_gene_list = pd.read_csv("%s/ranked_gene_list_Sensitive_versus_Resistant_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the ranked gene list from the gsea report folder
					average = ranked_gene_list['SCORE'].mean(axis=0) #this computes the average of the enrichment scores
					stdev = ranked_gene_list['SCORE'].std(axis=0) #this computes the stdev of the enrichment scores
					subtracted_scores = ranked_gene_list['SCORE'].sub(average) #this subtracts the average from the scores
					z_scores = subtracted_scores.divide(stdev) #this divides the subtracted scores by the stdev to give the z score 
					ranked_gene_list['Z_SCORES'] = z_scores #this adds the z-scores as their own column in the dataframe
				
					significant_downregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] > 2)] #this pulls out the genes who's z-scores are greater than 2
					names_significant_downregulated_data = significant_downregulated_data.loc[:,('NAME')]
					significant_downregulated_data.loc[:,('NAME')] = 'DN_' + names_significant_downregulated_data.astype(str) #this appends the prefix "DN" to the gene names
					significant_downregulated_genes = significant_downregulated_data['NAME'].tolist() #this pulls out the downregulated classified gene names and turns them into a list
					significant_upregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] < -2)] #this pulls out the genes who's z-scores are less than -2
					names_significant_upregulated_data = significant_upregulated_data.loc[:,('NAME')]
					significant_upregulated_data.loc[:,('NAME')] = 'UP_' + names_significant_upregulated_data.astype(str) #this appends the prefix "UP" to the gene names
					significant_upregulated_genes = significant_upregulated_data['NAME'].tolist() #this pulls out the upregulated classified genes and turns them into a list
				
				
				elif os.path.isfile("%s/ranked_gene_list_R_versus_S_%s.xls" %(folder_name, trial_num)):  #this classifies the genes if the enrichment report is for R vs S
					ranked_gene_list = pd.read_csv("%s/ranked_gene_list_R_versus_S_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the ranked gene list from the gsea report folder
					average = ranked_gene_list['SCORE'].mean(axis=0) #this computes the average of the enrichment scores
					stdev = ranked_gene_list['SCORE'].std(axis=0) #this computes the stdev of the enrichment scores
					subtracted_scores = ranked_gene_list['SCORE'].sub(average) #this subtracts the average from the scores
					z_scores = subtracted_scores.divide(stdev) #this divides the subtracted scores by the stdev to give the z score 
					ranked_gene_list['Z_SCORES'] = z_scores #this adds the z-scores as their own column in the dataframe
				
					significant_downregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] < -2)] #this pulls out the genes who's z-scores are less than -2
					names_significant_downregulated_data = significant_downregulated_data.loc[:,('NAME')]
					significant_downregulated_data.loc[:,('NAME')] = 'DN_' + names_significant_downregulated_data.astype(str) #this appends the prefix "DN" to the gene names
					significant_downregulated_genes = significant_downregulated_data['NAME'].tolist() #this pulls out the downregulated classified gene names and turns them into a list
					significant_upregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] > 2)] #this pulls out the genes who's z-scores are greater than 2
					names_significant_upregulated_data = significant_upregulated_data.loc[:,('NAME')]
					significant_upregulated_data.loc[:,('NAME')] = 'UP_' + names_significant_upregulated_data.astype(str) #this appends the prefix "UP" to the gene names
					significant_upregulated_genes = significant_upregulated_data['NAME'].tolist() #this pulls out the upregulated classified genes and turns them into a list
					
				else: #this takes into account the fact that gsea might have named it Resistant versus Sensitive instead of R versus S
					ranked_gene_list = pd.read_csv("%s/ranked_gene_list_Resistant_versus_Sensitive_%s.xls" %(folder_name, trial_num), delimiter='\t') #this retrieves the ranked gene list from the gsea report folder
					average = ranked_gene_list['SCORE'].mean(axis=0) #this computes the average of the enrichment scores
					stdev = ranked_gene_list['SCORE'].std(axis=0) #this computes the stdev of the enrichment scores
					subtracted_scores = ranked_gene_list['SCORE'].sub(average) #this subtracts the average from the scores
					z_scores = subtracted_scores.divide(stdev) #this divides the subtracted scores by the stdev to give the z score 
					ranked_gene_list['Z_SCORES'] = z_scores #this adds the z-scores as their own column in the dataframe
				
					significant_downregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] < -2)] #this pulls out the genes who's z-scores are less than -2
					names_significant_downregulated_data = significant_downregulated_data.loc[:,('NAME')]
					significant_downregulated_data.loc[:,('NAME')] = 'DN_' + names_significant_downregulated_data.astype(str) #this appends the prefix "DN" to the gene names
					significant_downregulated_genes = significant_downregulated_data['NAME'].tolist() #this pulls out the downregulated classified gene names and turns them into a list
					significant_upregulated_data = ranked_gene_list[(ranked_gene_list['Z_SCORES'] > 2)] #this pulls out the genes who's z-scores are greater than 2
					names_significant_upregulated_data = significant_upregulated_data.loc[:,('NAME')]
					significant_upregulated_data.loc[:,('NAME')] = 'UP_' + names_significant_upregulated_data.astype(str) #this appends the prefix "UP" to the gene names
					significant_upregulated_genes = significant_upregulated_data['NAME'].tolist() #this pulls out the upregulated classified genes and turns them into a list
					
				
				
			
				#add the two lists of genes together
				significant_dysregulated_genes = significant_downregulated_genes + significant_upregulated_genes
			
				#then we save the list of dysregulated genes in the dictionary for each drug
				list_significant_genes_cnv[drug_name] = significant_dysregulated_genes



#converts the dictionaries of the gene and pathway lists into dataframes and saves them as files
#it names the files based on the input files, if more than 3 databases are combined it's simply named combined databases

rna_dataframe_significant_pathways = pd.DataFrame.from_dict(list_significant_pathways, orient='index') #this converts the dictionary into a dataframe with the indexes as keys - which means that shorter columns have "NaN" added on to extend their length
transpose_rna_dataframe = rna_dataframe_significant_pathways.transpose() #this transposes the dataframe so that the keys are drug names again 
pathways_dataframe_nan = transpose_rna_dataframe.fillna(value='NaN') #this replaces empty cells with Nan
if num_databases == 1:
	pathways_dataframe_nan.to_csv('%s_Pathways_Results.txt' % (rna_basename), index=False, sep='\t')
elif num_databases == 2:
	pathways_dataframe_nan.to_csv('Combined_%s_%s_Pathways_Results.txt' %(list_rna_basenames[0], list_rna_basenames[1]), index=False, sep='\t')
elif num_databases == 3:
	pathways_dataframe_nan.to_csv('Combined_%s_%s_%s_Pathways_Results.txt' %(list_rna_basenames[0], list_rna_basenames[1], list_rna_basenames[2]), index=False, sep='\t')
else:
	pathways_dataframe_nan.to_csv('Combined_Databases_Pathways_Results.txt', index=False, sep='\t')
	
rna_dataframe_significant_genes = pd.DataFrame.from_dict(list_significant_genes_rna, orient='index') #this converts the dictionary into a dataframe with the indexes as keys 
transpose_rna_genes_dataframe = rna_dataframe_significant_genes.transpose() #this transposes the dataframe so that the keys are drug names again 
rna_genes_dataframe_nan = transpose_rna_genes_dataframe.fillna(value='NaN') #this replaces empty cells with Nan
if num_databases == 1:
	rna_genes_dataframe_nan.to_csv('%s_Gene_Results.txt' % (rna_basename), index=False, sep='\t')
elif num_databases == 2:
	rna_genes_dataframe_nan.to_csv('Combined_%s_%s_Gene_Results.txt' %(list_rna_basenames[0], list_rna_basenames[1]), index=False, sep='\t')
elif num_databases == 3:
	rna_genes_dataframe_nan.to_csv('Combined_%s_%s_%s_Gene_Results.txt' %(list_rna_basenames[0], list_rna_basenames[1], list_rna_basenames[2]), index=False, sep='\t')
else:
	rna_genes_dataframe_nan.to_csv('Combined_Databases_Gene_Results_rna.txt', index=False, sep='\t')

cnv_dataframe_significant_genes = pd.DataFrame.from_dict(list_significant_genes_cnv, orient='index') #this converts the dictionary into a dataframe with the indexes as keys 
transpose_cnv_genes_dataframe = cnv_dataframe_significant_genes.transpose() #this transposes the dataframe so that the keys are drug names again 
cnv_genes_dataframe_nan = transpose_cnv_genes_dataframe.fillna(value='NaN') #this replaces the empty cells with NaN
if num_databases == 1:
	cnv_genes_dataframe_nan.to_csv('%s_Gene_Results.txt' % (cnv_basename), index=False, sep='\t')
elif num_databases == 2:
	cnv_genes_dataframe_nan.to_csv('Combined_%s_%s_Gene_Results.txt' %(list_cnv_basenames[0], list_cnv_basenames[1]), index=False, sep='\t')
elif num_databases == 3:
	cnv_genes_dataframe_nan.to_csv('Combined_%s_%s_%s_Gene_Results.txt' %(list_cnv_basenames[0], list_cnv_basenames[1], list_cnv_basenames[2]), index=False, sep='\t')
else:
	cnv_genes_dataframe_nan.to_csv('Combined_Databases_Gene_Results_cnv.txt', index=False, sep='\t')



#This program makes a network of cancer drugs by comparing how many pathway elements are shared by any two drugs
def compare(drug1, drug2):
    return set(drug1).intersection(set(drug2))

def get_pathways(f, sep = "\t"):
    fh = open(f, "rU")
    s = '   '
    names = {}
    dataset = {}
    for line_num, line in enumerate(fh):
        line = line.strip().split(sep)
        #If this is the first line, get the drug names.
        if line_num == 0:
            for col, name in enumerate(line):
                names[col] = name
            print len(names)
        #Or else, collect the data
        else:
            for col, data in enumerate(line):
                if data != "NaN":
                    if names[col] in dataset.keys():
                        dataset[names[col]].append(data)
                        #print col
                        #print names[col] 
                        #print data
                    else:
                        dataset[names[col]] = [data] 
                        #print col
                        #print names[col]
                        #print data
    return dataset

# the next script uses the programs defined above to generate the edge weights for result text files (pathway, RNA_gene, and CNV_gene)

#for the pathway CADREs
if num_databases == 1:
	readfile = '%s_Pathways_Results.txt' % (rna_basename)
elif num_databases == 2:
	readfile = 'Combined_%s_%s_Pathways_Results.txt' %(list_rna_basenames[0], list_rna_basenames[1])
elif num_databases == 3:
	readfile = 'Combined_%s_%s_%s_Pathways_Results.txt' %(list_rna_basenames[0], list_rna_basenames[1], list_rna_basenames[2])
else:
	readfile = 'Combined_Databases_Pathways_Results.txt' 
network = open(readfile.split(".")[0]+"_network.txt", "w")
index = open(readfile.split(".")[0]+"_index.txt", "w")
network.write("Drug1\tWeight\tDrug2\n") #Write header
data = get_pathways(readfile, "\t")
names = data.keys()
for d1 in range(len(names)):
    for d2 in range(d1+1 ,len(names)):
        network.write(names[d1]+"\t"+str(len(compare(data[names[d1]], data[names[d2]])))+"\t"+names[d2]+"\n")
        index.write(names[d1]+"\t"+names[d2]+"\t"+str(", ".join(list(compare(data[names[d1]], data[names[d2]]))))+"\n")
network.close()
index.close()

#for the differential expression CADRE
if num_databases == 1:
	readfile = '%s_Gene_Results.txt' % (rna_basename)
elif num_databases == 2:
	readfile = 'Combined_%s_%s_Gene_Results.txt' %(list_rna_basenames[0], list_rna_basenames[1])
elif num_databases == 3:
	readfile = 'Combined_%s_%s_%s_Gene_Results.txt' %(list_rna_basenames[0], list_rna_basenames[1], list_rna_basenames[2])
else:
	readfile = 'Combined_Databases_Gene_Results_rna.txt'
network = open(readfile.split(".")[0]+"_network.txt", "w")
index = open(readfile.split(".")[0]+"_index.txt", "w")
network.write("Drug1\tWeight\tDrug2\n") #Write header
data = get_pathways(readfile, "\t")
names = data.keys()
for d1 in range(len(names)):
    for d2 in range(d1+1 ,len(names)):
        network.write(names[d1]+"\t"+str(len(compare(data[names[d1]], data[names[d2]])))+"\t"+names[d2]+"\n")
        index.write(names[d1]+"\t"+names[d2]+"\t"+str(", ".join(list(compare(data[names[d1]], data[names[d2]]))))+"\n")
		#this returns all drug pair combinations, even if edge weight equals zero
network.close()
index.close()

#for the cnv CADRE
if num_databases == 1:
	readfile = '%s_Gene_Results.txt' % (cnv_basename)
elif num_databases == 2:
	readfile = 'Combined_%s_%s_Gene_Results.txt' %(list_cnv_basenames[0], list_cnv_basenames[1])
elif num_databases == 3:
	readfile = 'Combined_%s_%s_%s_Gene_Results.txt' %(list_cnv_basenames[0], list_cnv_basenames[1], list_cnv_basenames[2])
else:
	readfile = 'Combined_Databases_Gene_Results_cnv.txt'
network = open(readfile.split(".")[0]+"_network.txt", "w")
index = open(readfile.split(".")[0]+"_index.txt", "w")
network.write("Drug1\tWeight\tDrug2\n") #Write header
data = get_pathways(readfile, "\t")
names = data.keys()
for d1 in range(len(names)):
    for d2 in range(d1+1 ,len(names)):
        network.write(names[d1]+"\t"+str(len(compare(data[names[d1]], data[names[d2]])))+"\t"+names[d2]+"\n")
        index.write(names[d1]+"\t"+names[d2]+"\t"+str(", ".join(list(compare(data[names[d1]], data[names[d2]]))))+"\n")
        #this returns all drug pair combinations, even if edge weight equals zero
network.close()
index.close()




#the following removes edges whose weights are three or less
#it combines the rna and cnv data
#and it removes edges whose weight is less than two standard deviations away from the mean

if num_databases == 1:
	pathways_network = pd.read_csv('%s_Pathways_Results_network.txt' %(rna_basename), delimiter='\t')
	if pathwaycutoff == 'TRUE': #this applies the threshold if the user wants it
		pathways_network_morethan3 = pathways_network[(pathways_network['Weight'] > 3)] #this removes edges whose weights are three or less from the pathways network
		pathways_network_morethan3.to_csv('Cytoscape_%s_PathwayCADRE.txt' %(rna_basename), index=False, sep='\t')
	else: #otherwise no threshold is used, just removes the edge weights that equal zero
		pathways_network_nozeros = pathways_network[(pathways_network['Weight'] > 0)]
		pathways_network_nozeros.to_csv('Cytoscape_%s_PathwayCADRE.txt' %(rna_basename), index=False, sep='\t')
elif num_databases == 2:
	pathways_network=pd.read_csv('Combined_%s_%s_Pathways_Results_network.txt' %(list_rna_basenames[0], list_rna_basenames[1]), delimiter='\t')
	if pathwaycutoff == 'TRUE': #applies the threshold if indicated by user
		pathways_network_morethan3 = pathways_network[(pathways_network['Weight'] > 3)] #this removes edges whose weights are three or less from the pathways network
		pathways_network_morethan3.to_csv('Cytoscape_Combined_%s_%s_PathwayCADRE.txt' %(list_rna_basenames[0], list_rna_basenames[1]), index=False, sep='\t')
	else: #otherwise no threshold is used
		pathways_network_nozeros = pathways_network[(pathways_network['Weight'] > 0)]
		pathways_network_nozeros.to_csv('Cytoscape_Combined_%s_%s_PathwayCADRE.txt' %(list_rna_basenames[0], list_rna_basenames[1]), index=False, sep='\t')
elif num_databases == 3:
	pathways_network = pd.read_csv('Combined_%s_%s_%s_Pathways_Results_network.txt' %(list_rna_basenames[0], list_rna_basenames[1], list_rna_basenames[2]), delimiter='\t')
	if pathwaycutoff == 'TRUE': #applies the threshold if indicated by user
		pathways_network_morethan3 = pathways_network[(pathways_network['Weight'] > 3)] #this removes edges whose weights are three or less from the pathways network
		pathways_network_morethan3.to_csv('Cytoscape_Combined_%s_%s_%s_PathwayCADRE.txt' %(list_rna_basenames[0], list_rna_basenames[1], list_rna_basenames[2]), index=False, sep='\t')
	else: #otherwise no threshold is used
		pathways_network_nozeros = pathways_network[(pathways_network['Weight'] > 0)]
		pathways_network_nozeros.to_csv('Cytoscape_Combined_%s_%s_PathwayCADRE.txt' %(list_rna_basenames[0], list_rna_basenames[1]), index=False, sep='\t')
else:
	pathways_network = pd.read_csv('Combined_Databases_Pathways_Results_network.txt' % (rna_basename), index=False, sep='\t')
	if pathwaycutoff == 'TRUE': #applies the threshold if indicated by user
		pathways_network_morethan3 = pathways_network[(pathways_network['Weight'] > 3)] #this removes edges whose weights are three or less from the pathways network
		pathways_network_morethan3.to_csv('Cytoscape_Combined_Databases_PathwayCADRE.txt', index=False, sep='\t')
	else: #otherwise no threshold is used 
		pathways_network_nozeros = pathways_network[(pathways_network['Weight'] > 0)]
		pathways_network_nozeros.to_csv('Cytoscape_Combined_Databases_PathwayCADRE.txt', index=False, sep='\t')



#opens the rna network file from ComparePathways
if num_databases == 1:
	rna_gene_network = pd.read_csv('%s_Gene_Results_network.txt' % (rna_basename), delimiter='\t')
elif num_databases == 2:
	rna_gene_network = pd.read_csv('Combined_%s_%s_Gene_Results_network.txt' %(list_rna_basenames[0], list_rna_basenames[1]), delimiter='\t')
elif num_databases == 3:
	rna_gene_network = pd.read_csv('Combined_%s_%s_%s_Gene_Results_network.txt' %(list_rna_basenames[0], list_rna_basenames[1], list_rna_basenames[2]), delimiter='\t')
else:
	rna_gene_network = pd.read_csv('Combined_Databases_Gene_Results_rna_network.txt', delimiter='\t')

#opens the cnv network file from ComparePathways
if num_databases == 1:
	cnv_gene_network = pd.read_csv('%s_Gene_Results_network.txt' % (cnv_basename), delimiter='\t')
elif num_databases == 2:
	cnv_gene_network = pd.read_csv('Combined_%s_%s_Gene_Results_network.txt' %(list_cnv_basenames[0], list_cnv_basenames[1]), delimiter='\t')
elif num_databases == 3:
	cnv_gene_network = pd.read_csv('Combined_%s_%s_%s_Gene_Results_network.txt' %(list_cnv_basenames[0], list_cnv_basenames[1], list_cnv_basenames[2]), delimiter='\t')
else:
	cnv_gene_network = pd.read_csv('Combined_Databases_Gene_Results_cnv_network.txt', delimiter='\t')


#this creates a new dataframe to store the rna drug pairs and the cnv drug pairs in
combined_gene_network = pd.merge(rna_gene_network, cnv_gene_network, how='outer', on=['Drug1', 'Drug2'], suffixes=('_rna', '_cnv')) #merges the rna and cnv data in a dataframe
combined_gene_network = combined_gene_network.fillna(value=0) #replaces any NaN's created by the merge with zeros

#this adds the weights of the rna and cnv drug pairs together, removes the edge weights of zero, and applies the edge weight cut off 
combined_gene_network['Combined_Weight'] = combined_gene_network['Weight_rna'] + combined_gene_network['Weight_cnv']
combined_gene_network_no_zeros = combined_gene_network[(combined_gene_network['Combined_Weight'] > 0)]
edited_combined_gene_network = pd.DataFrame()
edited_combined_gene_network['Drug1'] = combined_gene_network_no_zeros.loc[:,('Drug1')]
edited_combined_gene_network['Weight'] = combined_gene_network_no_zeros.loc[:,('Combined_Weight')]
edited_combined_gene_network['Drug2'] = combined_gene_network_no_zeros.loc[:,('Drug2')]
average = combined_gene_network_no_zeros['Combined_Weight'].mean(axis=0) #this computes the average of the weights
stdev = combined_gene_network_no_zeros['Combined_Weight'].std(axis=0) #this computes the stdev of the weights
cutoff = average + (2 * stdev) #this calculates the cutoff that we're going to apply
if featurecutoff == 'TRUE': #applies threshold if indicated by user
	average = combined_gene_network_no_zeros['Combined_Weight'].mean(axis=0) #this computes the average of the weights
	stdev = combined_gene_network_no_zeros['Combined_Weight'].std(axis=0) #this computes the stdev of the weights
	cutoff = average + (2 * stdev) #this calculates the cutoff that we're going to apply
	filtered_combined_gene_network = edited_combined_gene_network[(edited_combined_gene_network['Weight'] > cutoff)]
else: #otherwise doesn't apply threshold 
	filtered_combined_gene_network = edited_combined_gene_network #changes the name of the dataframe to the name that is later used to save to a file for cytoscape 


if num_databases == 1:
	filtered_combined_gene_network.to_csv('Cytoscape_%s_FeatureCADRE.txt' % (rna_basename), index=False, sep='\t')
elif num_databases == 2:
	filtered_combined_gene_network.to_csv('Cytoscape_Combined_%s_%s_FeatureCADRE.txt' %(list_rna_basenames[0], list_rna_basenames[1]), index=False, sep='\t')
elif num_databases == 3:
	filtered_combined_gene_network.to_csv('Cytoscape_Combined_%s_%s_%s_FeatureCADRE.txt' %(list_rna_basenames[0], list_rna_basenames[1], list_rna_basenames[2]), index=False, sep='\t')
else:
	filtered_combined_gene_network.to_csv('Cytoscape_Combined_Databases_FeatureCADRE.txt', index=False, sep='\t')
#this saves the gene network to a file that can be run by cytoscape



#if they don't give a folder for the gsea output, just use the current working directory
#figure out a way so that they don't have to give  a full name for the gsea downloads
#still can't call the gsea trial runs from the separate folder gsea_trials that I store them in, if this doesn't work just have them stored in the current working directory so the user doesn't have to pass an argument
