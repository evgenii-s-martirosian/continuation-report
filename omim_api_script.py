# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:55:43 2022

@author: 44784
"""

from datetime import datetime
import os, sys, getopt, requests
import argparse
import pandas as pd
import numpy as np
import re

chromatin_related_genes = pd.read_table("C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/gene-list-22-04-22-filtered-terms.tsv")

class Gene():
    # Here we create two arrays:
    # 1. One will store the hexadecimal code for each object of the class which
    # is stored in the memory.
    # 2. Another will store the memory code for genecount attribute of each 
    # object (genecount df).
    # This allows to make the class iterable for future code.
    _generegistry = []
    
    # Here we define function to store gene object as Gene Symbol. Other
    # Other attribute will be MIM number of associated phenotype/s.
    def __init__(self, geneid):
        self._generegistry.append(self)
        self.geneid = geneid
        self._phenotypeMIMnumber = list()
        
    # Here we define a function which will be used to fetch the data from OMIM
    def fetch_omim_data(self):
        server = "https://api.omim.org/"
        ext1 = "api/geneMap/search?search=" + self.geneid + "&"
        RequestFormat = "format=json" + "&"
        ApiKey = "apiKey=nWOiOxxtSfyqNneLIbZ5FQ"
        
        r1 = requests.get(server+ext1+RequestFormat+ApiKey, headers={ "Content-Type" : "application/json"})
        decoded = r1.json()
        gene_map_list = decoded['omim']['searchResponse']['geneMapList']
        gene_map_list_len = len(gene_map_list)
        for genemaplist in range(0, gene_map_list_len):
            try:
                if decoded['omim']['searchResponse']['geneMapList'][genemaplist]['geneMap']['approvedGeneSymbols'] == self.geneid:
                    try:
                        phenotype_map_list = decoded['omim']['searchResponse']['geneMapList'][genemaplist]['geneMap']['phenotypeMapList']#[0]['phenotypeMap']
                        phenotype_map_list_len = len(phenotype_map_list)
                        for number in range(0, phenotype_map_list_len):
                            self._phenotypeMIMnumber.append(phenotype_map_list[number]['phenotypeMap']['phenotypeMimNumber'])
                    except KeyError:
                        self._phenotypeMIMnumber = None
                    except IndexError:
                        self._phenotypeMIMnumber = None
                else:
                    pass
            except KeyError:
                pass

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
def create_omim_array():
	resultArr = []
	#counter = 0
	for geneobject in Gene._generegistry:
		print(geneobject.geneid)
		geneobject.fetch_omim_data()
		resultArr.append({'Geneid': geneobject.geneid,
                    'phenotypeMimNumber': geneobject._phenotypeMIMnumber})
		#print(resultArr)

	ensembl_df = pd.DataFrame(resultArr)
	return ensembl_df

geneid_list = chromatin_related_genes['Geneid'].tolist()

geneids = create_instances(Gene, 'geneid', geneid_list)

gene_omim_phenotypes = create_omim_array()

combined_df = pd.merge(chromatin_related_genes, gene_omim_phenotypes, on='Geneid')
#combined_df.to_csv('C:/Users/44784/OneDrive - The University of Manchester/PhD_Genomics/Chromatin-related genes/GO-terms-new-22-04-22/chromatin_related_genes_ensembl_omim.tsv', sep="\t", index=None)
#server = "https://api.omim.org/"
#ext1 = "api/geneMap/search?search=" + "MYC" + "&"
#RequestFormat = "format=json" + "&"
#ApiKey = "apiKey=nWOiOxxtSfyqNneLIbZ5FQ"
        
#r1 = requests.get(server+ext1+RequestFormat+ApiKey, headers={ "Content-Type" : "application/json"})
#decoded = r1.json()
#print(decoded['omim'])

#print(decoded['omim']['searchResponse']['geneMapList'][0]['geneMap']['approvedGeneSymbols'])

#phenotype_map = decoded['omim']['searchResponse']['geneMapList'][0]['geneMap']['phenotypeMapList']#[0]['phenotypeMap']
#phenotype_map_len = len(phenotype_map)
#print(phenotype_map[0:4])
#list_mim_n = list()
#for number in range(0, phenotype_map_len):
#    print(phenotype_map[number]['phenotypeMap']['phenotypeMimNumber'])
#    list_mim_n.append(phenotype_map[number]['phenotypeMap']['phenotypeMimNumber'])