#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 25 16:08:20 2019

@author: ismini
"""

import argparse
import contextlib
import io
#import requests
#import re
#import sys
#import matplotlib.pyplot as plt
#import matplotlib.patches as mpatch
#from collections import defaultdict
#import matplotlib.ticker as ticker
#from classes import*
from functions import*


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='The best python project ever by Ismini')
       
    parser.add_argument('-q', '--question', help='Your question in english')
    
    parser.add_argument('-t', '--task', help = 'The task in english please')

    args = parser.parse_args()
    question = args.question
    task = args.task
    
    if question:
        print ('The question asked is:', question)

        #find the target in question (gene,trancript etc)
                          
        try:
            
            target = check_input(question, hgnc)
            
            if type(target) == Transcript:
                
                #as basic transcript I'll take the transcript itself
                bs_trans = target.name
                
                #find the entries for that basic transcript
                transcript_components = target.db_entries(db=gencode)
                    
            if type(target) == Gene:
                
                #find the basic transcript
                bs_trans = target.basic_transcript()
                
                #find the entries for that basic transcript
                transcript_components = target.db_entries_basic_transcript(db=gencode)
                
            if type(target) == DNA_site:

                    #Silence print message of basic_transcript method, as it was printed at the class assignment
                    with contextlib.redirect_stdout(io.StringIO()):
                        
                        #find the basic transcript
                        bs_trans = target.basic_transcript()
                    
                    #find entries for that basic transcript
                    transcript_components = target.db_entries_basic_transcript(db=gencode)
              
                
            if type(target) == Polymorphism:
                
                if type(target.site) == Gene:
                    
                    #find the basic transcript
                    bs_trans = target.site.basic_transcript()
                    
                    #find the entries for that basic transcript
                    transcript_components = target.site.db_entries_basic_transcript(db=gencode)
                    
                    #position of the mutation
                    CDSs = [range(int(y[3]), int(y[4])+1) for y in [x.split('\t') for x in transcript_components] if y[2]== 'CDS']
                    
                    if transcript_components[0].split('\t')[6] == '+':
                        
                        mut_pos = [v for k,v in enumerate([y for x in map(list,CDSs) for y in x]) if k == int(re.findall('\d+',target.mutation[0])[0])+1 ][0]
                    
                    if transcript_components[0].split('\t')[6] == '-':
                        
                        mut_pos = [v for k,v in enumerate(sorted([y for x in map(list,CDSs) for y in x], reverse = True)) if k == int(re.findall('\d+',target.mutation[0])[0])+1 ][0]
                    
                    
                if type(target.site) == Transcript:
                    
                    bs_trans = target.site.name
    
                    transcript_components= target.site.db_entries(db=gencode)
    #####here is the problem: transcript_components = []. I must match refseq id to ensemble id, but I can't match them
    #####exactly because if map refseq to ensg_id first and then take the principal transcript, there is a chance to choose
    #####a transcript that isn't a perfect match to that refseq id (as being told in http://www.ensembl.info/2019/04/16/we-are-making-mane-changes/)                
    #### This choice will result in showing a variant outside the range of the chosen transcript, if not annotated as MANE in gencode or ensembl (it is mentioned that almost 
    #### 50% of transcripts are annotated as MANE and thus ensuring for a perfect match between ids)               
    
    #### solution: 
    #### target.site as Gene() in check_input function and handle its ensg_id as another 
    #### type of Class.Gene input. Ask ensembl to map it on a gene. And then continue as if it was a gene input                
                                        
                    CDSs = [range(int(y[3]), int(y[4])+1) for y in [x.split('\t') for x in transcript_components] if y[2]== 'CDS']
                    if transcript_components[0].split('\t')[6] == '+':
                        
                        mut_pos = [v for k,v in enumerate([y for x in map(list,CDSs) for y in x]) if k == int(re.findall('\d+',target.mutation[0])[0])+1 ][0]
                    
                    if transcript_components[0].split('\t')[6] == '-':
                        
                        mut_pos = [v for k,v in enumerate(sorted([y for x in map(list,CDSs) for y in x], reverse = True)) if k == int(re.findall('\d+',target.mutation[0])[0])+1 ][0]
                    
                if type(target.site) == DNA_site: 
                    
                    #Silence print message of basic_transcript method, as it was printed at the class assignment
                    with contextlib.redirect_stdout(io.StringIO()):
                        
                        #find the basic transcript
                        bs_trans = target.site.basic_transcript()
                    
                    #find entries for that basic transcript
                    transcript_components = target.site.db_entries_basic_transcript(db=gencode)
                    
                    #position of the mutation
                    mut_pos = int(target.site.pos)
                    
                if type(target.site) == Dbsnp_var:
                    
                    #find the basic transcript
                    bs_trans = target.site.basic_transcript()
                    
                    #find entries for that basic transcript
                    transcript_components = target.site.db_entries_basic_transcript(db=gencode)
                    
                    #position of the variant
                    
                    mut_pos = int(target.site.pos)
                    
            #find gene name
            gene = re.findall('gene_name=(.+);transcript_type', transcript_components[0])[0]
        
            #find strand
            if transcript_components[0].split('\t')[6] == '+':
                strand = 'forward strand'
            if transcript_components[0].split('\t')[6] == '-':
                strand = 'reverse strand'
            
            #find chromosome
            if re.findall('chr(\d+)', transcript_components[0].split('\t')[0]) != []:
                
                chromosome = re.findall('chr(\d+)', transcript_components[0].split('\t')[0])[0]
            else:
                chromosome = re.findall('chr(.)', transcript_components[0].split('\t')[0])[0]
            
            #define mutation if any    
            if hasattr(target, 'mutation'):
                mut = target.mutation
                
            else:
                mut = None
                mut_pos = None
                
            
            #merge fields of each entry from transcript_components
            aa = Gencode_List(transcript_components)
    
            #dictionary with transcript regions
            map_regions= aa.map_regions(8,8)               
                
            #exon/intron map of basic trancript    
            map_plot(map_regions, bs_trans, gene, strand, chromosome, mut, mut_pos)
            
        except Exception as e:
            print(e)
            
    if task:
        print("The task assigned is:", task)
    
    
