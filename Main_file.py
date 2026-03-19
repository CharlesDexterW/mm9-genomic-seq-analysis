#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 25 16:56:13 2025

@author: Benjamin Garcés Cifuentes
email='agarces2381@gmail.com'

"""
#%%
'''
# Sequence Analysis
Importing necessary packages & Functions
'''

import csv # This will provide functionality for reading and writing tabular data in CSV (Comma Separated Values) format.
import gc # This will provide an interface to the garbage collector, allowing manual control over memory management.
import io # This will provide Python's main facilities for handling various types of I/O (input/output), including file I/O, string I/O, and in-memory streams.
import zipfile # This  will allow you to create, read, write, append, and list files in a ZIP archive.

from LoadFASTA_Function import LoadFastaFile, LoadGene, TSSChroms

#%% Loading genome into a variable
gene_file= 'mm9_sel_chroms_knownGene.txt' 
file_open = open(gene_file).readlines() 

gene_inf=LoadGene(gene_file) # creating a dictionary
#%%
'''
To access any needed information such as the chromosome and position of 
the gene on the chromosome, just use brackets like this:
'''
gene_inf['uc009auw.1']['chr']

#%% Evaluating the Data
chroms=[] # First we have to create a list

for k in gene_inf.keys(): # For each gene, we look to see if the chromosome is in list.  If not, add it:
  chr=gene_inf[k]['chr']
  if chr not in chroms:
    chroms=chroms+[chr]

#%% Counting Genes
gene_counts={} # this codeline creates a dictionary to store the relevant gene information:
for chr in chroms: # Now this new dictionary must be filled using two for loops: one for each  chromosome and one for each gene:
  chrom_count=0 
  for k in gene_inf.keys():
    if gene_inf[k]['chr']==chr:
      chrom_count+=1
  gene_counts[chr]=chrom_count

#%% Locating identified Genes in Genomic Data
fasta_file='selChroms_mm9.fa.zip' #Now the raw gene sequences can be read from the FASTA file.
seq_dict=LoadFastaFile(fasta_file) # Now to load the sequences into a dictionary. This might take a minute...or two depending on your computer’s RAM and processor.
cntn4="uc009dcr.2" # Now, you can start using the gene information within the chromosomes. 

#%%
gene_inf[cntn4]['chr'] # we can find out the chromosomal location of Cntn4 using our gene_info dictionary

#%% Dictionary of Data
len(seq_dict['chr6']) # Use len() to see the length of any chromosome, like chromosome 6.

#%% Extracting Genetic Information
'''
You can peek into the entire chromosomal sequence where we know Cntn4 to reside:
'''
hchr=gene_inf[cntn4]['chr']
hst=gene_inf[cntn4]['start']
hen=gene_inf[cntn4]['end']
cntn4_seq=seq_dict[hchr][hst:hen] # This variable pulls out the sequence of interest.

cntn4_seq[5:200] # Each sequence is a string of characters now, you can identify specific parts of that gene.  
#%% We can identify where translation begins in the gene sequence using the index command, built into the Python string library:

cntn4_seq.index('ATG')
#%% Genomic Statistics
"""
Here I create data structures. So far it works properly with python 3.12.7
Using the gene_inf dictionary, find out the length of every gene on the four chromosomes collected. 
First, import the NumPy module to include more complicated math calculations:
if by any chance you find yourself with any errors, try running "pip install numpy" In the terminal
"""

import numpy as np
#%%
gene_lengths={} # Let's find the start and end of each gene, and store the difference:
for g in gene_inf.keys():
        st=gene_inf[g]['start']
        en=gene_inf[g]['end']
        gene_lengths[g]=np.absolute(en-st)

import matplotlib.pyplot as plt # By using pyplot's basic plotting tools in matplotlib we create a histogram:
"""If there's any problem importinb matplotlib try running "pip install matplotlib" in the terminal"""
plt.hist(gene_lengths.values(), bins=50, log=True, facecolor='green')

#%% basic mathematic operations with genes.
"""
To compare lengths of genes we can call out a gene by its identifier in the dictionary and add or substract it from 
another gene sequence such as cntn4.
"""
gene_lengths["uc012enb.1"]-len(cntn4_seq)

#%% Beyond the Coding Region
"""
Let’s use the gene_inf data structure to count the base pairs in genes compared to those that aren’t in coding regions.
By using a set data structure, which stores unique items, you can remove duplicates. To make sure there are no double countings, let’s track the indices, of chromosome x so as not to add in a position in the list of coding regions if they’re already considered in another gene this way.
"""
# This can be done using a boolean array.

chr6_len=len(seq_dict['chr6'])

ingene_numpy=np.zeros(chr6_len,dtype=bool)

for gene in gene_inf.keys():
    if gene_inf[gene]['chr']=='chr6':
        start_in=gene_inf[gene]['start']
        end_in=gene_inf[gene]['end']
        ingene_numpy[start_in:end_in]=True

ingene_numpy
#%% Counting the coding base pairs

sum_gene=ingene_numpy.sum() # Now by summing this array, we will find out the number of coding index sites (lenght of coding sequence in chr6)
print(sum_gene)

len(ingene_numpy)-sum_gene # Amount of non coding base-pairs of the chromosome 6
(sum_gene/len(ingene_numpy))*100 # Fractions of non-coding DNA

#%%
chr6_starts=TSSChroms(gene_inf,'chr6') # Finding TATA motif which work as regulatory regions
'''
This dictionary contains the starting positions of every gene in the chromosome.
By default, the TATA binding motif relevant to transcriptions is located in a small window of the starting sites.
Let's give it a try with a 40 base pair window upstream.
Now with this dictionary, we can identify which genes contain a TATA binding motif near the starting site, which means
these genes can be regulated with a TATA binding protein.
'''
tata_dis = {}
for g in chr6_starts.keys():
    e = chr6_starts[g]
    strand = gene_inf[g]['strand']
    if strand == '+':
        s = e - 40
        if 'TATA' in seq_dict['chr6'][s:e].upper():
            tata_dis[g] = seq_dict['chr6'][s:e].upper().rindex('TATA')
    else:
        s = e
        e = s + 40
        if 'TATA' in seq_dict['chr6'][s:e].upper():
            tata_dis[g] = seq_dict['chr6'][s:e].upper().index('TATA')

len(tata_dis) # Now how many genes have a motif within 40 bps of each transcription starting site?

#%%  Calculate the mean distance of the transcription start sites and the TATA motif
if tata_dis:
    mean_tata_dist = sum(tata_dis.values()) / len(tata_dis)
    print(f'Mean TATA distance from TSS: {mean_tata_dist:.2f} bp')
else:
    print('No TATA motifs found within the search window — mean distance cannot be calculated.')
 