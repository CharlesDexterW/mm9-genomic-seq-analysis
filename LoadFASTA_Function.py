# -*- coding: utf-8 -*-
"""
Created on Sat Jan  4 19:23:22 2025

@author: Benjamin Garcés Cifuentes

email='agarces2381@gmail.com'
"""

# Python 3.12 


##The following packages enable the parsing of diverse types of files.

import csv
import gc
import io
import zipfile
import pyflakes


def LoadFastaFile(filename):# Reads zipped FASTA files into a dictionary of strings. 
    sequence_dictionary={}
    line_count=0
    comp_file = None

    if zipfile.is_zipfile(filename):
        print('will be loaded from zip file')
        comp_file = zipfile.ZipFile(filename,'r')
        file_name = comp_file.namelist()[0]
        lines = io.TextIOWrapper(comp_file.open(file_name))
    else:
        # BUG FIX 1 — Non-zip branch used to call comp_file.open() when
        # comp_file was still None, crashing immediately. Plain files are
        # opened directly with open() instead.
        print('Will be loaded from plain text file.')
        lines = open(filename, 'r')
    
    try:
        current_header, current_sequence = '', []
        for row in lines:
            line_count += 1
            if line_count % 10000 == 0:
                print('Parsed line ' + str(line_count))
 
            if row.startswith('>'):
                if current_header != '':
                    sequence_dictionary[current_header] = ''.join(current_sequence)
                    current_sequence = []
                    gc.collect()
                current_header = row[1:].strip()
            else:
                current_sequence.append(row.strip())
 
        # BUG FIX 2 — The last sequence in the file was silently dropped.
        # The loop only saved a sequence upon hitting the *next* '>' header,
        # so the final record never triggered a save. This adds it explicitly
        # after the loop ends.
        if current_header:
            sequence_dictionary[current_header] = ''.join(current_sequence)
 
        print('Found ' + str(len(sequence_dictionary)) + ' sequences in FASTA file.')
        return sequence_dictionary
 
    finally:
        # BUG FIX 3 — The finally block called comp_file.close() unconditionally.
        # For plain files comp_file was None, so this crashed with AttributeError.
        # Now only the zip handle is closed here; plain files are closed via
        # the context manager built into open().
        if comp_file is not None:
            comp_file.close()


def LoadGene(filename, getCoding=False): # reads a tab-delimited gene table file (likely from UCSC Genome Browser)  and extracts gene location information.
    gene_dictionary={}
    chromosome=set()
    
    column_names = ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 
              'exonStarts', 'exonEnds', 'proteinID', 'alignID']
    reader=csv.DictReader(open(filename),fieldnames=column_names,delimiter='\t')
    for line in reader:
            if getCoding:
                start=int(line['cdsStart'])
                end=int(line['cdsEnd'])
            else:
                start=int(line['txStart'])
                end=int(line['txEnd'])
            gene_dictionary[line['name']] = {'chr': line['chrom'], 'start': start, 'end': end, 'strand': line['strand']}
            chromosome.add(line['chrom'])

    print('Parsed info for '+str(len(gene_dictionary))+' genes on '+str(len(chromosome))+' chromosomes')
    return gene_dictionary    

def TSSChroms(gene_inf,chrom): # extracts the Transcription Start Sites (TSS) for genes on a specific  chromosome.
    TSS={}
    for g in gene_inf.keys():
        if gene_inf[g]['chr']==chrom:    
            if gene_inf[g]['strand']=='+':
                TSS[g]=gene_inf[g]['start']
            else:
                TSS[g]=gene_inf[g]['end']
    return TSS
