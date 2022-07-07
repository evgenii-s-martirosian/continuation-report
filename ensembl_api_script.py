# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:59:47 2022

@author: 44784
"""

from datetime import datetime
import os, sys, getopt, requests
import argparse
import pandas as pd
import numpy as np
import re

chromatin_related_genes = pd.read_table("C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/Final-chromatin-related-genes-list.txt")

# Here we create gene class which will store information for each gene object of this class
class Gene():
	# Similarly as for previous class we create an array where positional information in memory
	# for gene object will be stored
	_generegistry = []

	# Here we define function to store gene object as 'ENSEMBL gene ID. Other attributes are
	# None at this point as they will be fetched from the ENSEMBL Rest API.
	def __init__(self, geneid):
		self._generegistry.append(self)
		self.geneid = geneid
		self._name = None
		self._ensemblid = None
		self._chr = None
		self._start = None
		self._end = None
		self._strand = None
		self._description = None
		self._canonicaltranscript = None

	# Here we define a function which will be used to fetch the gene data from the ENSEMBL Rest API.
	def fetch_ensembl_info(self):
		# Firstly, we store server address and extension for the database as variables.
		# We use lookup search instead of overlap as it will only show the information for the specific
		# ENSEMBL gene ID. However, overlap may show information for genes which positionally overlap with the
		# requested ENSEMBL gene ID, which is not suitable.
		server = "https://rest.ensembl.org"
		ext1 = "/lookup/symbol/homo_sapiens/" + self.geneid + "?"
		# Here, we send the request to get gene information from the API server in the JSON format
		# which is returned as HTML, which we then convert to human-readable JSON.
		r1 = requests.get(server+ext1, headers={ "Content-Type" : "application/json"})
		decoded = r1.json()

		# Here we create an if-elif conditional statement, which sends instruction on how to handle the requests.
		# If there is an error appears during the request, it might be due to the absence of gene ID in
		# in the database (might be a gene encoding uncharacterised protein). If this is the case, then
		# the gene ID is passed an the function looks at the next one.
		# The elif statment is set to only extract genes which have protein_coding biotype, which means
		# information from protein coding gene will only be stored. The else statement (other than protein_coding
		# genes) passess the gene ID and function checks the next gene ID.
		if 'error' in decoded:
			pass
		else: #decoded['biotype'] == 'protein_coding':
			# Here, we first try and extract gene information from the API request for each gene, which are stored
			# object attribtes.
			try:
				self._name = decoded['display_name']
				self._ensemblid = decoded['id']
				self._chr = decoded['seq_region_name']
				self._start = decoded['start']
				self._end = decoded['end']
				self._strand = decoded['strand']
				self._description = decoded['description']
				self._canonicaltranscript = decoded['canonical_transcript']

			# KeyError may arise when some information is not available or absent. Most likely, gene name information
			# or gene description maybe absent, as strand, start, end and chromosome number are usually known from NGS
			# data. This exception statement has 2 if-else conditional statements which check if gene name and/or description
			# are missing and sets their respective object attributes as None type.
			except KeyError:
				self._ensemblid = decoded['id']
				self._chr = decoded['seq_region_name']
				self._start = decoded['start']
				self._end =  decoded['end']
				self._strand = decoded['strand']

				if decoded.get('display_name') is None:
					self._name = 'None'
				else:
					self._name = decoded['display_name']

				if decoded.get('description') is None:
					self._description = 'None'
				else:
					self._description = decoded['description']
				if decoded.get('canonical_transcript') is None:
					self._canonicaltranscript = 'None'
				else:
					self._canonicaltranscript = decoded['canonical_transcript']
		#else:
		#	pass

# This function is used to create class objects by looking at the list of
# object which are needed to be the obects of the particular class.
# For each list object we parse it into the class __init__ function
# to make it a class object. Furthermore, we create a dictionary with these objects.
def create_instances(class_name, class_type, list_type):
	name = {class_type: class_name(class_type) for class_type in list_type}
	return name

# This function builds the array with the gene information from which dataframe is
# created needed for Task3. Firstly, we make empty array which will store gene information.
# Then, we have a for loop which iterates over each object of a gene class
# and adds gene information in to the array which are stored as the attributed of the Gene
# object. After all gene objects were iterated, the array is converted into the oandas dataframe.
def create_ensembl_array():
	resultArr = []
	#counter = 0
	for geneobject in Gene._generegistry:
		print(geneobject.geneid)
		geneobject.fetch_ensembl_info()
		if geneobject._chr is not None:
			resultArr.append({'Geneid': geneobject.geneid,
					  'Gene Name': geneobject._name,
                      'Ensembl ID': geneobject._ensemblid,
					  'Chromosome': geneobject._chr,
					  'Start': geneobject._start,
					  'End': geneobject._end,
					  'Strand': geneobject._strand,
					  'Description': geneobject._description,
                      'CanonicalTranscript': geneobject._canonicaltranscript})
		else:
			resultArr.append({'Geneid': geneobject.geneid,
    					  'Gene Name': geneobject.geneid,
                          'Ensembl ID': 'N/A',
    					  'Chromosome': 'N/A',
    					  'Start': 'N/A',
    					  'End': 'N/A',
    					  'Strand': 'N/A',
    					  'Description': 'N/A',
                          'CanonicalTranscript': 'N/A'})
		#print(resultArr)

	ensembl_df = pd.DataFrame(resultArr)
	return ensembl_df

geneid_list = chromatin_related_genes['OFFICIAL_GENE_SYMBOL'].tolist()
geneids = create_instances(Gene, 'geneid', geneid_list)
ensembl_data_array = create_ensembl_array()
#combined_list = pd.merge(ensemb_data_array, chromatin_related_genes, on='Geneid')

ensembl_data_array.to_csv('C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/chromatin-related-genes-ensembl.tsv', sep="\t", index=False)