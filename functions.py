#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 16:08:57 2019

functions called by the project

@author: ismini
"""

#import re

import matplotlib.pyplot as plt
import matplotlib.patches as mpatch
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
from math import ceil
from classes import*

def check_input(question, hgnc):
    
    '''chech words for biological meaning in a question'''
    
    #keep only letters, numbers and special characters needed
    rx = re.compile('[^a-zA-Z0-9 ._:>]')
    
    #Replace all characters except those in regex with '' (nothing)
    question = rx.sub('', question)
    
    #All words in question in uppercase
    question = question.upper()
    
    #match gene names with their ensembl id
    hgnc_gene_names = [(x[1],x[19]) for x in hgnc[1:]]
      
    result = ''
    
    while result == '':
        
        #catch an input like A>G at chromosome 12 position 12345
        if re.findall(".>.+\d+",question):
            result = re.findall(".>.+\d+",question)[0]
            print("The polymorphism you entered is {}".format(result))
            site = DNA_site(':'.join(re.findall("\d+", question)))
            mut = re.findall(".>.", question)
            target = Polymorphism(mut, site)
            break
        
        #catch an input like chromosome 12 position 12345 A>G or chr12 position 12345 A>G
        if re.findall("CHR.+>.",question):
            result = re.findall("CHR.+>.",question)[0]
            print("The polymorphism you entered is {}".format(result))
            site = DNA_site(':'.join(re.findall("\d+", question)))
            mut = re.findall(".>.", question)
            target = Polymorphism(mut, site)
            break
        
        for x in question.split(' '):
            
            if re.findall("(.+):(.+)>.",x):
                result = x
                print("The polymorphism you entered is {}".format(x))
           
                #catch an input like <Γονίδιο>:c.100A>G 
                if re.findall("(.+):", x)[0] in (gene for gene_names in hgnc_gene_names for gene in gene_names if gene != ''):
                    site = Gene(name = re.findall("(.+):", x)[0], data = hgnc)
                    mut = re.findall(".(\d+.>.)", x)
                    target = Polymorphism(mut, site)
                    break
                    
                #catch an input like <ENST Μετάγραφο>:c.100A>G
                if re.findall("(ENST.+):", x):
                    print('found a polymorphism at an ensembl transcript')
                    site = Transcript(re.findall("(.+):", x)[0])
                    mut = re.findall(".(\d+.>.)", x)
                    target = Polymorphism(mut, site)
                    break
                
                #catch an input like <refseq Μετάγραφο>:c.100A>G
                if re.findall("(NM_.+):", x):
                    print("found a polymorphism at a refseq transcript")
                    site = Gene(re.findall("(.+):", x)[0], data = hgnc)
                    mut = re.findall(".(\d+.>.)", x)
                    target = Polymorphism(mut, site)
                    break
                
                #catch an input like chr1:1234567A>G or 1:1234567A>G                  
                if re.findall("(.+):(.+)>.",x):
                    print('found a polymorphism at a chromosome site')
                    
                    site = DNA_site(':'.join(re.findall("\d+", x)))
                    mut = re.findall(".>.", x)
                    target = Polymorphism(mut, site)
                    break
                
            #catch a dbsnp variant
            if re.findall("(RS\d+)",x):
               
                result = x 
                print('The snp variant you entered is {}'.format(x.lower()))
                print("Awaiting a response from myvariant db")
                parameters = {
                            'q': x.lower(),
                            'fields': {'dbsnp' : 'gene'} 
                            }                

                url = 'http://myvariant.info/v1/query'
                
                response = requests.get(url, params=parameters)
                
                if response.ok:
                    myvariant_data = response.json()
                    
                if not response.ok:
                    response.raise_for_status()
                
                rs_id = x.lower()
                mut = re.findall("\d+(.>.)", myvariant_data['hits'][0]['_id'])
                try:
                    gene = myvariant_data['hits'][0]['dbsnp']['gene']['symbol']
                except Exception as e:
                    print('There is no {} matching this variant'.format(e))
                    sys.exit()

                site = Dbsnp_var(gene, mut, rs_id)
                
                target = Polymorphism(mut, site)
                break
            
            #catch a transcript with an ensembl id                            
            if re.findall("ENST(\d+)", x):
                result = x
                print('The transcript you entered is {}'.format(x))
                target = Transcript(result)
                break 
            #catch a refseq transcript
            if re.findall("NM_\d+", x):
                result = x
                print('The transcript you entered is {}'.format(x))
                target = Gene(result, data = hgnc)
                break
                        
            if re.findall("(.+):(.+)",x):
                
                rx = re.compile('[^0-9:]')
                result = rx.sub('', x)

                print("The site you entered is {}".format(x))
                target = DNA_site(result)
                break
                
            if x in (gene for gene_names in hgnc_gene_names for gene in gene_names if gene != ''):
                result = x
                print("The gene you entered is {}".format(x))
                target = Gene(result, data = hgnc)
                break
            
            #catch an input like chromosome 12 position 12345
            if re.findall("CHR.+",question):
                result = re.findall("CHR.+",question)[0]
                print("The site you entered is {}".format(result))
                
                target = DNA_site(':'.join(re.findall("\d+", question)))
                break
        
        
        if result == '':
            raise Exception("Wrong input. Check your question.")
    return(target)


def map_plot(d, bs_trans, gene, strand, chromosome, mut = None, mut_pos = None):
    
    '''plot an intron/exon map of a given dictionary containing genomic regions
       bs_trans is the basic transcript'''

    #start and end of transcript
    start_trans = float( d['transcript'][0]['start'])
    end_trans = float( d['transcript'][0]['end'])
    
    #figure properties
    fig=plt.figure(frameon = False)
    fig.size=(60,4)
    fig.dpi = 200
    
    ax=fig.add_subplot(111)
     
    #plt.rc('font', size=0.1 * fig.dpi) 
    #fontsize = 0.4 * fig.dpi
   
    #set axes range
    x_min = start_trans
    x_max = end_trans
    y_min = 0
    y_max =1
    ax.set_xbound(lower = x_min-2, upper = x_max+2)
    ax.set_ybound(lower = y_min, upper = y_max)
    ax.get_xaxis().get_major_formatter().set_scientific(False)
   

    #draw canvas to set axes
    fig.canvas.draw()
    
    #set maximum 3 xticks
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax.tick_params(direction='out', length=6, width=2, colors='k',
               grid_color='k', grid_alpha=0.5, labelsize=7)
   
    limx = ax.get_xlim()
    limy = ax.get_ylim()
    
    #add a xtick for the mutation
    if mut:    
        ax.set_xticks(list(ax.get_xticks()) + [int(start_trans),int(end_trans), mut_pos])
    else:
        ax.set_xticks(list(ax.get_xticks()) + [int(start_trans),int(end_trans)])
        
    plt.xticks(rotation=45, ha = 'right')
   
    #ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('{:,}'))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter('{:,.0f}'.format))
    ax.set_xlim(limx)
    ax.set_ylim(limy)
    
    
    #hide axes
    
    #ax.xaxis.set_visible(True)
      
    ax.yaxis.set_visible(False)
   
    
    #initial y position on y axis where rectangles for all components will be drawn
    Y = .6
    
    #add strand notation
    plt.text(start_trans, .02, 'chromosome ' + chromosome + ' ' + strand, ha = 'left', fontsize = 10)
    
    
    
    
    #initialize a list to add rectangles when plotted , so as to use them at mouse hovering   
    rectangles = []
    #d = map_regions 
    first_line = []
    second_line = []
    third_line = []
    fourth_line = []
    fifth_line = []
    
    if strand == 'reverse strand':
        dd = [v for k,v in sorted(list(enumerate(d.items())), reverse = True)]
        
        
    else:
        dd = d.items()
        
    
    for component, sites in dd:

        if component =='transcript':
            
           
            #plot a line indicating the transcript's length
            plt.plot((start_trans, end_trans), (Y + 0.1, Y + 0.1), 'k--', 
                     label = "transcript", lw = 0.4, zorder=1)
            
            #plot transcript's name
            plt.text(start_trans + (end_trans-start_trans)/2, 0.85 , 'transcript: ' + bs_trans, 
                     ha = 'center', fontsize = 15)
            
            #plot gene's name
            plt.text(start_trans + (end_trans-start_trans)/2, 0.90 , 'gene: ' + gene, 
                     ha = 'center', fontsize = 15)
            
            #add background horizontal line indicating transcript's range at UTR plot region
            plt.plot((start_trans, end_trans), 
                                     (Y -.1- .15, Y -.1 - .15), 'k--', lw = 0.4, zorder=1)
            
            #add background horizontal line indicating transcript's range at start/stop codon region
            plt.plot((start_trans, end_trans), 
                                     (Y -.1- 2 * .15, Y -.1 - 2 * .15), 'k--', lw = 0.4, zorder=1)
            
            #add background horizontal line indicating transcript's range to CDS' plot region
            plt.plot((start_trans, end_trans), 
                                     (Y - .1, Y - .1), 'k--', lw = 0.4, zorder=1)

                        
                
        else:
            
            #length of each drawing box (text drawing map) for each region in sites
            l_dbox = 7

            n_reg = len(sites)-1 # I don't want to plot a region for exons, only CDS, UTRs and start/stop codons
            first_line.append(u'\u250c' + u'\u2500'*(n_reg*l_dbox + (n_reg-1)*2) + u'\u2510')
            
            span = n_reg*l_dbox + (n_reg-1)*2
            s = ceil((span -len(component))/2)
            if len(component)==5:
                second_line.append(u'\u2502' +(' '*s) + component + (' '*s) + u'\u2502')
            if len(component) == 6:
                second_line.append(u'\u2502' +(' '*s) + component + (' '*(s-1)) + u'\u2502')
            
            #a list to store drawing boxes for each site and then append them in third_line
            sdb = []
            
            #the same for fourth & fifth line
            sdb4 = []
            sdb5 = []
            
            ddd = sorted(sites, key = lambda x: x['end'])

            for regions in ddd:

                #start and end of each exon, CDS, UTR, start/stop codon
                start = float(regions['start'])
                end = float(regions['end'])
                
                
                #height of rectangles 
                rect_height = 0.1
                
                if regions['region_id'] == 'exon':
                    #add rectangle to transcript's plot region
                    bb = mpatch.Rectangle((start, Y), end-start, 0.2, 
                                fc = 'none', alpha = 1, ec= 'k', linewidth=0.01, zorder = 2)
                    ax.add_patch(bb)
                        
                     
                    #add text
                    #plt.text(start + (end-start)/2, Y-0.05, regions['region_id'], ha = 'center', wrap=True)
                    
                    rectangles += bb,
                    
                    
                if regions['region_id'] == 'UTR5' or regions['region_id'] == 'UTR3':
                    #add rectangle to transcript's plot region
                    bb = mpatch.Rectangle((start, Y), end-start, 0.2, 
                                fc = 'olive', alpha = 1, ec= 'none', linewidth=0.01, zorder = 2)
                    ax.add_patch(bb)
                    
                    #add rectangle to UTR plot region
                    bb1 = mpatch.Rectangle((start, Y - 2*0.15), end-start, rect_height, 
                                fc = 'olive', alpha = 1, ec= 'olive', linewidth=0.5, zorder = 2)
                    ax.add_patch(bb1)
                    
                                        
                    #add text
                    plt.text(start + (end-start)/2, Y + 0.005 - 2*0.1, regions['region_id'], ha = 'center', 
                             fontsize=10)
                    
                    rectangles += bb,bb1
                    
                    if mut:
                    
                        if (mut_pos>= start and mut_pos<= end):
                            
                            sdb4.append(u'\u2502' +(' '*3) + u'\u2506' + (' '*3) + u'\u2502')
                            sdb5.append(u'\u2514' +(u'\u2500'*2) + re.findall('(.>.)',mut[0])[0] + (u'\u2500'*2) + u'\u2518')
                            
                        else:
                            sdb4.append(u'\u2502' + (' '*7) + u'\u2502')
                            sdb5.append(u'\u2514' +  (u'\u2500'*7) + u'\u2518')
                            
                    if not mut:
                        sdb4.append(u'\u2514' + (u'\u2500'*7) + u'\u2518') 
                    
                    sdb.append(u'\u2502' +'  ' + regions['region_id'] + ' ' + u'\u2502')
                    
                if regions['region_id'] == 'start_codon' or regions['region_id'] == 'stop_codon':
                    #add rectangle to transcript's plot region
                    bb = mpatch.Rectangle((start, Y), end-start, 0.2, 
                                      fc = 'red', alpha = 1, ec = 'red', linewidth=0.8, zorder=3)
                    ax.add_patch(bb)
                    
                    #add rectangle to start/stop codon plot region
                    bb1= mpatch.Rectangle((start, Y - 3*0.15), end-start, rect_height, 
                                fc = 'red', alpha = 1, ec= 'red', linewidth=0.8, zorder = 2)
                    ax.add_patch(bb1)
                    
                    
                    #add text
                    plt.text(start + (end-start)/2, Y- 2*0.1 - .145, regions['region_id'].replace('_',' '), 
                             ha = 'center', fontsize=10)
      
                    rectangles += bb,bb1
                    
                    if re.findall('(.+)_',regions['region_id'])[0] == 'start':
                        sdb.append(u'\u2502' +' ' + re.findall('(.+)_',regions['region_id'])[0] + ' ' + u'\u2502')
                        
                    if re.findall('(.+)_',regions['region_id'])[0] == 'stop':
                        sdb.append(u'\u2502' +'  ' + re.findall('(.+)_',regions['region_id'])[0] + ' ' + u'\u2502')
                        
                    if mut:
                    
                        if (mut_pos>= start and mut_pos<= end):
                            
                            sdb4.append(u'\u2502' +(' '*3) + u'\u2506' + (' '*3) + u'\u2502')
                            sdb5.append(u'\u2514' +(u'\u2500'*2) + re.findall('(.>.)',mut[0])[0] + (u'\u2500'*2) + u'\u2518')
                            
                        else:
                            sdb4.append(u'\u2502' + (' '*7) + u'\u2502')
                            sdb5.append(u'\u2514' +  (u'\u2500'*7)  + u'\u2518')
                            
                    if not mut:
                        sdb4.append(u'\u2514' + (u'\u2500'*7) + u'\u2518')  
                        
                if regions['region_id'] == 'CDS':
                    
                    #add rectangle to transcript's plot region
                    bb = mpatch.Rectangle((start, Y), end-start, 0.2, 
                                      fc = 'grey', alpha = 1, ec = 'none', linewidth=0.01, zorder=2)
                    ax.add_patch(bb)
                    
                    #add rectangle to CDS' plot region
                    bb1 = mpatch.Rectangle((start, Y - 0.15), end-start, rect_height, 
                                fc = 'grey', alpha = 1, ec= 'none', linewidth=0.8, zorder = 2)
                    ax.add_patch(bb1)
                    
                    rectangles += bb,bb1

                    sdb.append(u'\u2502' +'  ' + regions['region_id'] + '  ' + u'\u2502')
                    
                    if mut:
                    
                        if (mut_pos>= start and mut_pos<= end):
                            
                            sdb4.append(u'\u2502' +(' '*3) + u'\u2506' + (' '*3) + u'\u2502')
                            sdb5.append(u'\u2514' +(u'\u2500'*2) + re.findall('(.>.)',mut[0])[0] + (u'\u2500'*2) + u'\u2518')
                            
                        else:
                            sdb4.append(u'\u2502' + (' '*7) + u'\u2502')
                            sdb5.append(u'\u2514' + (u'\u2500'*7) + u'\u2518')
                            
                    if not mut:
                        sdb4.append(u'\u2514' + (u'\u2500'*7) + u'\u2518')                         

                #generate annotation to be displayed on mouse hovering over rectangles
                annot = ax.annotate("", xy=(0,0), xytext=(0,15),textcoords="offset points",
                                    fontsize = 10, ha ='center')
                   
                annot.set_visible(False)
            if sdb !=[]:
                third_line.append(''.join(sdb))
            if sdb4 != []:
                fourth_line.append(''.join(sdb4))
            if sdb5 != []:
                fifth_line.append(''.join(sdb5))
            
    #plot mutation if any
    if mut:
        
        if mut_pos<start_trans or mut_pos>end_trans:
            print("Variant falls outside basic transcript's region, at position: {} ".format(mut_pos))
       
        plt.text(mut_pos +1, .1, re.findall('(.>.)',mut[0])[0], ha = 'left', fontsize = 10, color = 'b')

        bb = mpatch.Rectangle((mut_pos, 0.0), 1, 0.8, 
                                    fc = 'b', alpha = 1, ec= 'b', linewidth=0.4, zorder = 4)
        ax.add_patch(bb)     
    
#mouse hovering following https://stackoverflow.com/questions/50560525/how-to-annotate-the-values-of-x-and-y-while-hovering-mouse-over-the-bar-graph        
    #update annotation for each rectangle          
    def update_annot(bb):
        #x,y coordinates of annotation
        x = bb.get_x()+bb.get_width()/2.
        y = bb.get_y()+bb.get_height()
        annot.xy = (x,y)
        #formating text annotation
        start = bb.get_x()
        end =bb.get_x() + bb.get_width()
        text = "{:,} - {:,}".format( int(start), int(end) )
        annot.set_text(text)
            
    #control display of annotations when hovering over each rectangle
    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            for bb in rectangles:
                cont, ind = bb.contains(event)
                if cont:
                    update_annot(bb)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                    return
        if vis:
            annot.set_visible(False)
            fig.canvas.draw_idle()
            
    #use "motion_notify_event" to control annotation display
    fig.canvas.mpl_connect("motion_notify_event", hover)
    
    #hide plot frame
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    ###add legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)

    utr_patch = mpatch.Patch(fc = 'olive', alpha = 1, ec= 'olive', linewidth=0.01, label='UTRs')
    codon_patch = mpatch.Patch(fc = 'red', alpha = 1, ec= 'red', linewidth=0.01, label = 'start/stop codons')
    cds_patch = mpatch.Patch(fc = 'grey', alpha = 1, ec= 'grey', linewidth=0.01, label = 'CDSs')
    introns_line = mlines.Line2D([], [], color='k', linewidth=0.6, linestyle='--', label='non-coding regions')

    handles=[utr_patch, codon_patch, cds_patch, introns_line]
    if mut:
        mut_line = mlines.Line2D([], [], color='b', linewidth=0.8, linestyle='-', label='polymorphism')
        handles += [mut_line]
    ll = len(handles)    
    
    plt.legend(handles= handles, bbox_to_anchor=(0.1, 0.02, 0.87, .102), loc='lower center',
           ncol=ll, mode="expand", borderaxespad=0., bbox_transform=plt.gcf().transFigure)
    
    #make space for x axis
    plt.tight_layout()
    
    #print text map

    nums = [i for i in range(0,len(first_line),9)]  
    
    if len(nums)==1:
        
         print(*first_line, sep = '  ')
         print(*second_line, sep = '  ')
         print(*third_line, sep = '--')
         print(*fourth_line, sep = '  ')
         if fifth_line != []:
            print(*fifth_line, sep = '  ')
    
    if len(nums)>1:
        y = nums[0]
        print(*first_line[y:y+9],'  ', sep = '  ')
        print(*second_line[y:y+9], '  ', sep = '  ')
        print(*third_line[y:y+9], '--', sep = '--')
        print(*fourth_line[y:y+9], '  ', sep = '  ')
        if fifth_line != []:
            print(*fifth_line[y:y+9], sep = '  ')
                
        for i in nums[1:-1]:
            print('  ', *first_line[i:i+9],'  ', sep = '  ')
            print('  ', *second_line[i:i+9], '  ', sep = '  ')
            print('--',*third_line[i:i+9], '--', sep = '--')
            print('  ',*fourth_line[i:i+9], '  ', sep = '  ')
            if fifth_line != []:
                print('  ', *fifth_line[i:i+9], sep = '  ')
        y = nums[-1]
        print('  ', *first_line[y:y+9], sep = '  ')
        print('  ', *second_line[y:y+9], sep = '  ')
        print('--', *third_line[y:y+9], sep = '--')
        print('  ', *fourth_line[y:y+9], sep = '  ')
        if fifth_line != []:
            print('  ',*fifth_line[y:y+9], sep = '  ')
    
        
                
        
    #print(*first_line, sep = '  ')
    #print(*second_line, sep = '  ')
    #print(*third_line, sep = '--')
    #print(*fourth_line, sep = '  ')
    #if fifth_line != []:
        #print(*fifth_line, sep = '  ')
       
                      
    #plot figure    
    plt.show()
    

if __name__ == '__main__':
   question = 'Where is gene SF3B1:c.100C>G?'
   check_input(question,hgnc)
   map_plot(map_regions, bs_trans, gene, strand, chromosome, mut, mut_pos)
   
