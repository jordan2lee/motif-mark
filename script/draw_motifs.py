#!/usr/bin/env python3

import cairo
import math
import re
import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description='draw seq motifs to show motif location in intron or exon. handles ambiguous characters in motif sequences and uracil')
    parser.add_argument("-fa", "--fasta_path", help ="path to fasta file that want to scan seq for motifs", required=True, type=str)
    parser.add_argument("-m", "--motif_file_path", help ="path to file of motifs of interest. must have one motif seq per line", required=True, type=str)
    parser.add_argument("-t", "--temp_fa", help ="must include name of output file. path where want to create temp fasta. is reformatted so entire seq on 1 line", required=True, type=str)
    parser.add_argument("-svg", "--svg_path", help ="must include name. path to save output svg file", required=True, type=str)
    return parser.parse_args()
    
args = get_arguments()
fa_path = args.fasta_path
motif_path = args.motif_file_path
temp_path = args.temp_fa
svgfile = args.svg_path


#####################################
############# Functions #############
#####################################

def format_fa(path_to_fasta, path_to_outputfile):
    '''Condenses seq info to a single line. Outputs a copy of input file, but with all seq info on one line. Each read consists of only 2 lines (header and seq). Output file does not necessarily have to be already created'''
    
    # Open file and loop through file (one line at a time)
    with open(path_to_fasta, "r") as in_fh, open(path_to_outputfile, "w") as temp_fh:

        file_line1 = 'true' # first line in entire file
        
        for line in in_fh:
            
            # Write very first line in file to output
            if file_line1 == 'true':
                temp_fh.write(line)
                # update that no longer looking at first line in entire file
                file_line1 = 'false'
            
            # Looking at any line that isn't first line in file
            else:
                # If seq line --> strip whitespace and write to output
                if line[0] != ">":
                    seq_line = line.strip()
                    temp_fh.write(seq_line)

                # If header line --> add newline char and write to output 
                if line[0] == ">":
                    temp_fh.write("\n" + line) 
                    
        # Proper format: end of file - has last seq line end with \n 
        temp_fh.write("\n")

    return "Fasta reformatting completed. Temp file located at:", temp_path

## Test above function works ##
# fa_path="../data/testfile_INSR.fasta"
# temp_path="../data/temp.fa"
# format_fa(fa_path, temp_path) 

def get_nt_type(nt):
    '''Given a single nt, outputs string if exon (capital letter) or intron (lowercase letter)'''
    if nt.isupper() == True:
        current_type = 'exon'
        
    if nt.islower() == True:
        current_type = 'intron'
        
    return str(current_type)

## Test above function works #
# get_nt_type("z")

def ambig_found(target_motif_list):
    '''Input list of motif(s) that want to find in a read seq later on.
    Checks if only ATCG or atcg in motif seq. Returns False if only ATCG 
    or 'true' if any ambiguous characters found 
    (e.g. Y, etc).'''
    import re
    # Code below will exit function at first sight of 
        # any motif with ambi char
        # even if not finished looping through motif list
        
    i = 0
    for motif in target_motif_list:
        
        match = re.search(r"[^ATCGatcg]", motif) # match = match object in re

        # If motif has ambig character(s)
        if match:
            return 'true'

        # If only AGTC
        else:
            # look at next motif in list for ambig char
            i+=1
    return 'false'


## Test above function works ##
# test_list = ['agag', 'aaaq']
# ambig_found(test_list)

# must already created a dict containing ambig and their meanings = called ambig_dna_dict

def replace_ambig_char(motif_list):
    '''Input motif list and return mod_motif_list with these converted ambig characters, checks if ambig char present,
    outputs motif without ambig char as per regex.
    Motif seq can be capital or lowercase
    motif = yaa outputs to motif = [AT]aa. Assumes you have
    already created ambig_dna_dict where keys = ambig char
    and values = all possible nt that ambig char could be'''
    import re
    
    mod_motif_list = []

    i = 0
    for motif in motif_list:
        for key in ambig_dna_dict:
            if key in motif:
            # flags = allows to work for lower or capital letters, e.g. input yaa, sees key=y value=A or T, output [AT]aa
                meaning = ambig_dna_dict[key]
                
                mod_motif = re.sub(key, meaning, motif, flags = re.IGNORECASE)
                motif = mod_motif
                
                

        mod_motif_list.append(motif)
        i+=1

    return mod_motif_list

## Test above function works ##
# test_list = ['aaa', 'yat', 'ggg']
# replace_ambig_char(test_list)

def mot_coords(motif, nt_seq):
    '''scan if a motif in seq, then return start/stop coords, draw box, continue scanning for motif'''
    import re
    
    # list to store all start/stop coords for the given motif
    coords = []
    
    # look for motif in nt_seq, flags makes it work for both capital and lowercase letters
    match = re.finditer(motif, nt_seq, flags=re.IGNORECASE)
        
    if match:
        i=0
        for item in match:
            start = item.start()+1
            stop = item.end()
            coords.append(start)
            coords.append(stop)
            i+=1
        return coords
    
    else:
        return # returns empty list if motif not found in seq
    
    
## Test if above function works ##
# test_motif = 'gggg'
# test_nt_seq = 'atataggggtatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatataggggtatatatccccatatatatatatatatatatatatatatatAAAAACATACCCCGGCACTGGTGTCGAGGGGATAGatatataatatatatatatatatatatatatatatatatatatatatatatatatatatatatatccccatatatatatatatatatatatatatatatatatatatatatatatatataggggtatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatatccttccttccttctttctttttttccttctttctttccttctttcttttc'
# mot_coords(test_motif, test_nt_seq)

################# Arg parse #################
# user needs to tell how many total reads in file, num_rds
num_rds = 2



################# Obj set up or comment out for arg parse #################
#fa_path="../data/testfile_INSR.fasta"
# fa_path="../data/testfile_INSR.1read.fasta"
#temp_path="../data/temp.fa"
# temp_path="../data/temp.1read.fa"
#motif_path="../data/testfile_motifs.txt" # contains y
# motif_path="../data/motifs.txt" # contiains y


################# Reformat fa: all seq info on 1 line  #################        
# Create temp file that is copy of original fa
#but reformatted for prot seq on one line only 
# Will use temp file as fasta for all downstream work
format_fa(fa_path, temp_path) 


################# Save motifs of interest in motif_list #################
with open(motif_path, "r") as motif_fh:
    # Initialize list
    motif_list = []
    # Add motifs to list
    for line in motif_fh:
        line = line.strip()
        motif_list.append(line)


################# Create_ambig_dict,  key = nt code, value = tuple of nts it represents #################
# ambig_dna_dict contains for ambiguous nt meanings, 
    #e.g. Y represents a nt that can be a cytosine or thymine.
    #Only works for Y, R, W, S, K, M, D, , H, B, X, and N.
    #Will populate dict containing lower and uppercase versions
    
# Create dict ambig_dna, 
ambig_dna_dict = {}
nt_key = ["U", "u", "Y", "y", "R", "r", "W", "w", "S", "s", "K", "k", \
            "M", "m", "D", "d", "V", "v", "H", "h", "B", "b", \
            "X", "x", "N", "n"]
repre = ["[T]", "[t]", '[CT]', '[ct]', '[AG]', '[ag]', '[AT]', '[at]', '[GC]', '[gc]', '[TG]', '[tg]', \
         '[CA]', '[ca]', '[AGT]', '[agt]', '[ACG]', '[acg]', '[ACT]', '[act]', '[CGT]', '[cgt]', \
         '[ATCG]', '[atcg]', '[ATCG]', '[atcg]']
# Populate ambig_dna_dict
i = 0
for symbol in nt_key:
    meaning = repre[i]
    # Populate dict
    ambig_dna_dict[symbol] = meaning
    i+=1 
# print(ambig_dna_dict)


################# Set up for using pycairo #################
# Outline for all remaining code:
#Notes: 1)loop through one read and plot exons/introns, 2)loop thorugh same read and plot motifs, 3)move to next read and repeat

# Set up 'canvas' where drawing image
width, height = 800, 500
# Create coordaintes where to display image 
    # produces SVG file
#svgfile = "../output/motif_location_image.svg"
surface = cairo.SVGSurface(svgfile, width, height)
# Create create the coordinates where I will be drawing on
context = cairo.Context(surface)
# Width of line
context.set_line_width(1)
# All set to draw image downstream

################# Read file once: Draw intron line #################

with open(temp_path, 'r') as mod_fa:
    # set starting pt
    x_init = 5
    y_init = 160 
    rd_counter = 1
    
    while True:
        header = mod_fa.readline().strip()
        seq = mod_fa.readline().strip()
        if not header: # if end of file
            break  
        

        if seq:
            ## 1. Set up Canvas ##
            # Scale canvas, width = full width of canvas
            draw_width = width - 10
            # Determine length of seq for a given read and incorporate into scaling
            seq_length = len(seq)
            adj_seq_w = draw_width / seq_length # how many canvas units represent 1 nt
            
            ## 2. Draw horizontal line that is as long as total nts in seq (then scaled) ##
            # how far down to start, arbitrarily picked 4
            y = y_init + (40 *(rd_counter -1))
#             y = y_init + (rd_counter * 4)
            # Move to starting spot of line
            context.move_to(x_init, y)
            # Move to ending spot of line (aka how forare right to move)
            x = (x_init + seq_length) * adj_seq_w # x = (total nts) * canvas_units_per_1_nt
            context.line_to(x, y)
            # Draw it
            context.stroke()
            # Update rd counter
            rd_counter +=1
            
################# Read file Second Time: Draw exon boxes #################
with open(temp_path, 'r') as mod_fa:
    # set starting pt
    x_init = 5
    y_init = 160
    rd_counter = 1
    font_x_init = 5
    font_y_init = 130
    
    while True:
        header = mod_fa.readline().strip()
        seq = mod_fa.readline().strip()
        if not header: # if end of file
            break
        
        if header:
            ## Write out legend ##
            context.select_font_face("Times New Roman", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            context.set_font_size(16)
            x_leg = 5
            y_leg = 25
            context.move_to(x_leg, y_leg)
            context.show_text('Legend') #add legend name
            # Exon legend
            context.move_to(x_leg, y_leg+20)
            context.rectangle(x_leg, y_leg+20, 5, 20)
            context.fill()
            context.move_to(x_leg+10, y_leg+35)
            context.show_text('Exon')
            # Intron legend
            context.move_to(x_leg+80, y_leg+30)
            context.line_to(x_leg+110, y_leg+30)
            context.stroke()
            context.move_to(x_leg+120, y_leg+35)
            context.show_text('Intron')
            # Motif legend
            col = 1
            move = 0
            for motif in motif_list:
                # Draw colored box
                context.move_to(x_leg+move, y_leg+55)
                context.rectangle(x_leg+move, y_leg+55, 5, 20)
                context.set_source_rgb(0.45*col, 0.25*col, 0.15)
                context.fill()
                # write text
                context.move_to(x_leg+10+move, y_leg+70)
                context.set_source_rgb(0,0,0)
                context.show_text(motif)
                col +=0.75
                move+=100

            ## Type out header info to output ##
            context.set_source_rgb(0,0,0)
            context.select_font_face("Times New Roman", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            context.set_font_size(13)
            font_y = font_y_init + (40 *(rd_counter -1) + 15)
            context.move_to(font_x_init, font_y)
            context.show_text(header)
        
        if seq:
            ## 1. Set up Canvas ##
            # Scale canvas, width = full width of canvas
            draw_width = width - 10
            # Determine length of seq for a given read and incorporate into scaling
            seq_length = len(seq)
            adj_seq_w = draw_width / seq_length # how many canvas units represent 1 nt
#             print('rd_counter is', rd_counter)
            
    
            ## 2. Draw exons ##
            i = 0
            start = 'no exon yet'
            stop = 'null'
            for nt in seq:
                if start == 'no exon yet':
                    # If intron, then move on to next pos
                    if get_nt_type(nt) == 'intron':
                        i+=1
                    # If exon
                    if get_nt_type(nt) == 'exon':
                        start = i
                        i+=1
                elif stop == 'null':
                    # if still see exon
                    if get_nt_type(nt) == 'exon':
                        i+=1
                    # if now see introns
                    if get_nt_type(nt) == 'intron':
                        stop = i
#                         print(start, stop)
                        
                        # Now have start and stop coordinates of exon..let's draw
                        context.move_to(x_init, y_init)
                        # Draw rect cords. 
                        #x0 --> x coord of top left rect corner
                        #y0 --> y coord of top left rect corner
                        #x1 --> rectangle lenght (how far right to go)
                        #y1 --> rectangle height (how far down to go)   
#                         print('exon start and stop coords:', start, stop)
                        x0 = start * adj_seq_w # multip by canvas units per 1 nt
                        x1 = (stop-start) * adj_seq_w # multip by canvas units per 1 nt
                        y0 = (y_init + (40 *(rd_counter -1))) - 10
                        y1 = 20
                        context.rectangle(x0, y0, x1, y1)
                        context.fill()
                        # Reset for next nt in seq
                        start = 'no exon yet'
                        stop = 'null'
                        i+=1
                  
                
            # Update things before moving to next read in fadd
            rd_counter +=1
                           
                
################# Read file Third time: Draw motif boxes #################

with open(temp_path, 'r') as mod_fa:
    # set starting pt
    x_init = 5
    y_init = 160
    rd_counter = 1
    
    while True:
        header = mod_fa.readline().strip()
        seq = mod_fa.readline().strip()
        if not header: # if end of file
            break
        
        if seq:
            
            ## 4. Replace any ambig char that exist in motifs of interest ##
            if ambig_found(motif_list) == 'true':
                # Create mod_mofit_list that has ambigs as regex, y --> [ct]
                motif_list = replace_ambig_char(motif_list)
#                 print(motif_list)
                
                
                
            ## 5. Draw motif ##
            col=1
            for motif in motif_list:
#                 create function: scan if a motif in seq, then return start/stop coords, draw box, continue scanning for motif
                start_stop = mot_coords(motif, seq)
#                 print(start_stop)
                
                
                
                even = 0
                odd = even + 1
        
        
                for cood in start_stop:
                    if odd < len(start_stop):
#                         print('even and odd is', even, odd)
                        start = start_stop[even]
                        stop = start_stop[odd]
#                         print(start, stop)

                        x0 = start * adj_seq_w # multip by canvas units per 1 nt
                        x1 = (stop-start) * adj_seq_w # multip by canvas units per 1 nt
                        y0 = (y_init + (40 *(rd_counter -1))) - 10
                        y1 = 20
                        context.rectangle(x0, y0, x1, y1)
                        context.set_source_rgb(0.45*col, 0.25*col, 0.15)
                        context.fill()

                        #update 
                        even +=2
                        odd+=2
                col+=0.75
                        
            rd_counter+=1


                


                
# Image complete
surface.finish()