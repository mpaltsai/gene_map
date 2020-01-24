#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 16:13:20 2019

classes used in the project

@author: ismini
"""
import re
from collections import defaultdict
import requests
import sys
import os

#download hgnc file with gene names, ids etc if not in directory
if not os.path.exists('protein-coding_gene.txt'):
       
    print('downloading gene list file from HGNC')
    os.system("wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/protein-coding_gene.txt")

#read HGNC list of genes
with open('protein-coding_gene.txt') as fi:
    
    hgnc = fi.read()

hgnc = hgnc.strip().split('\n')

hgnc = [x.split('\t') for x in hgnc]

#there are some withdrawn entries, so we must exclude them
hgnc = [x for x in hgnc if x[5] != 'Entry Withdrawn']

#download gff3 file from gencode
if not os.path.exists('gencode.v29.annotation.gff3'):
       
    print('downloading gene list file from HGNC')
    os.system('wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gff3.gz')
    
    os.system('gunzip gencode.v29.annotation.gff3.gz')
#read gencode file   
with open("gencode.v29.annotation.gff3") as f:
    
    gencode = f.read()

#read lines of gencode file
gencode = gencode.strip()
gencode = gencode.split('\n')

class Gene:    
    def __init__(self, name, data):
        self.name = name
        
        try:
            self.ensg_id = list(set([x[19] for x in data if x[1] == self.name]))[0]
        except Exception:
            
            parameters = {
              'fields': 'genomic_pos',
              'species': 'human',
              'q' : 'refseq:' + self.name
            }
            url = 'http://mygene.info/v3/query'
    
            response = requests.get(url, params=parameters)

   
            if response.ok:
                mygene_data = response.json()
                if mygene_data['hits'] == []:
                    raise Exception("Can't find a gene this refseq transcript matches")
                mygene_data = response.json()
                ensg_id = mygene_data['hits'][0]['genomic_pos']['ensemblgene']
                self.ensg_id = ensg_id 
            if not response.ok:
                response.raise_for_status()
                
        
        
    def db_entries(self, db):
        db=gencode
        return  [line for line in db if self.ensg_id in line]
    
    def basic_transcript(self,):
        #read tab separated entries
        res = [x.split('\t') for x in self.db_entries(db=gencode)]
    
        #expand field 8, with information of each entry
        #res_expand = [x[8].split(';') for x in res]
        
        #find all transcripts, with a tag=basic (appris system) and the smaller 
        #value in transcript support level, for a given gene
        gene_transcripts = [(re.findall(r"transcript_id=(ENST.+);gene_type", x[8]), 
                             re.findall(r"transcript_support_level=(\d)", x[8]) ) 
                for x in res if re.match(r"(.+)basic,appris_principal_(\d),CCDS",x[8])]
   
        if gene_transcripts != []:   
            #take the first transcript for this site
            return sorted(gene_transcripts, key=lambda x: x[1])[0][0][0]
        else:
            raise Exception("This region doesn't fall into a gene region or the region of the basic transcript")  
    
    
    def db_entries_basic_transcript(self, db):
        
        basic_transcript = self.basic_transcript() 
    
        return [x for x in self.db_entries(db=gencode)  if basic_transcript in x]

class Transcript:
    def __init__(self, name):
        self.name = name

    def db_entries(self, db):
        
        db = gencode
        return [line for line in db if self.name in line]
       
      
class Polymorphism:
    def __init__(self, mut, site):
        self.mutation = mut
        self.site = site
        
class Dbsnp_var:
    def __init__(self, gene, mut, rs_id):
        data = hgnc
        self.gene = gene
        self.mut = mut
        self.rs_id = rs_id

        try:
            self.ensg_id = list(set([x[19] for x in data if x[1] == self.gene]))[0]
        except Exception:
            print("Can't find an ensembl gene_id")
            print("This variant doesn't match to a gene region")
            sys.exit('Choose a variant that matces a gene')
        
         
        if self.ensg_id:
            
            print("Awaiting a response from ensembl db")
            #retrieve coordinates of variant from ensembl using the ensembl gene id
         
            server = "http://rest.ensembl.org"
            ext = "/overlap/id/" + self.ensg_id + "?feature=variation"
            r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
             
            if not r.ok:
              r.raise_for_status()
              raise Exception("This variant doesn't fall within a gene region")
              
            else:
                response = r.json()        
            
            self.pos = [x['start'] for x in response if x['id'] == self.rs_id][0]
        
           
    def db_entries(self, db):
        db=gencode
        return  [line for line in db if self.ensg_id in line]
        
    def basic_transcript(self,):
        #read tab separated entries
        res = [x.split('\t') for x in self.db_entries(db=gencode)]
    
        #find all transcripts, with a tag=basic (appris system) and the smaller 
        #value in transcript support level, for a given gene
        gene_transcripts = [(re.findall(r"transcript_id=(ENST.+);gene_type", x[8]), 
                             re.findall(r"transcript_support_level=(\d)", x[8]) ) 
                for x in res if re.match(r"(.+)basic,appris_principal_(\d),CCDS",x[8])]
   
        if gene_transcripts != []:   
            #take the first transcript for this site
            return sorted(gene_transcripts, key=lambda x: x[1])[0][0][0]
        else:
            raise Exception("This region doesn't fall into a gene region or the region of the basic transcript")  
       
            
    def db_entries_basic_transcript(self, db):
        
        basic_transcript = self.basic_transcript() 
    
        return [x for x in self.db_entries(db=gencode)  if basic_transcript in x]
        
class DNA_site:
    def __init__(self, site):
        self.site = site
        self.chr = re.findall('(.+):', self.site)[0]
        self.pos = re.findall(':(\d+)', self.site)[0]
                
        
    def basic_transcript(self,):
        db = gencode
        ##take entries for the specific chromosome
        chrom = self.chr
        res = [line for line in db if re.findall('chr'+(chrom+'\t'),line) and not re.search("^#",line)]
        
        #find those entries that our site falls within
        res_untab = [x.split('\t') for x in res]
        
        res_entries_site = [x for x in res_untab if int(self.pos) >=int(x[3]) and int(self.pos) <= int(x[4])]
         
        #find all transcripts, with a tag=basic (appris system) and the smaller 
        #value in transcript support level, for a given gene
        transcripts = [(re.findall(r"transcript_id=(ENST.+);gene_type", x[8]), 
                             re.findall(r"transcript_support_level=(\d)", x[8]) ) 
                for x in res_entries_site if re.match(r"(.+)basic,appris_principal_(\d),CCDS",x[8])]
        if transcripts != []:   
            
            if len(transcripts) > 1:
                print("An exon/intron map for transcript {} will be shown".format(sorted(transcripts, key=lambda x: x[1])[0][0][0]))
                print("If you wish to see a map for the rest transcripts, you should provide them as a question input, one at a time\n These transcripts are:")
                rest_transcripts = [x[0][0] for x in sorted(transcripts, key=lambda x: x[1])[1:] ]
                print(rest_transcripts)
                
            if len(transcripts) <= 1:
                pass
                
            #take the first transcript for this site
            return sorted(transcripts, key=lambda x: x[1])[0][0][0]
        
            
        else:
            raise Exception("This region doesn't fall into a gene region or the region of the basic transcript")
            
            
                
    def db_entries_basic_transcript(self, db):
        
        basic_transcript = self.basic_transcript() 
    
        return  [line for line in db if basic_transcript in line]
      
class Gencode_List: 
        
        '''a class to handle lists derived from gencode'''
        
        def __init__(self, a_list):
            self.list = a_list
            
        def tab_sep(self,):
            return  [x.split('\t') for x in self.list]
        
        def semicolon_sep(self,index):
            return  [x[index].split(';') for x in self.tab_sep()]
        
        def merge_fields(self,fields,index):
            return [x+y for x,y in zip([x[0:fields] for x in self.tab_sep() ], self.semicolon_sep(index)) ]
            
        def map_regions(self,fields,index):
            
            '''make a dictionary with regions in a transcript so as
                to use them when plotting the intron/exon map'''
                
            map_regions = defaultdict(list)
    
            for x in self.merge_fields(fields,index):
                region_id = re.findall('ID=(.+):ENST', x[8])
                if region_id == []:
                    
                    map_regions[x[2]].append({'start' : x[3], 'end' : x[4]})
                
                else:
                    exon = re.findall('(.+)_.+=(\d+)', x[16])[0]
                    exon = ''.join(exon)
                    region_id = re.findall('ID=(.+):ENST', x[8])[0] 
                    map_regions[exon].append({'region_id' : region_id,  'start' : x[3], 'end' : x[4]})
            
            return map_regions 