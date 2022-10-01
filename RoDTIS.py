#!/usr/bin/env python

"""
Documentation: https://github.com/nielsend/RoDTIS 
Note: Assumes a space separates the accession from other NCBI values 
(e.g.: >AATDWI010000001.1 Escherichia coli strain 197080 SAMN05750861-rid6247213.denovo.001, whole genome shotgun sequence)

"""
import re
import sys
import csv
import math
import pandas as pd
#new 10.01.2022
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Check to see if a value is a float
def isfloat(value):
  try:
    float(value)
    return True
  except:
    return False


#Read in file
parameter = sys.argv[1]
if parameter=='-s':
#    print(sys.argv[2])
    if isfloat(sys.argv[2]) and -0.000000000000000000001 < float(sys.argv[2]) < 1.0000000000000000000000001:
        threshold = sys.argv[2]
        inFiles = sys.argv[3:]
        start=3
    else:
        raise TypeError('Please enter a vailid number between 0.0-1.0')

else: 
    parameter = 0.1
    threshold = parameter
    inFiles = sys.argv[1:]
    start=1

sequence = ''
genome= ''


def abs_list_search(list1,list2,min_allow,max_allow):
    #check that lists are not empty 
    if list1==-1 or list2==-1:
        return False
    #if both lists are singleton lists, job is easy
    if len(list1)==1 and len(list2)==1:
        list1=list1[0];
        list2=list2[0];
        return bool(min_allow<=abs(list1-list2)<=max_allow)
    else: #else, sad day:must iterate through all matches
        for item1 in list1:
            for item2 in list2:
                if bool(min_allow<=abs(item1-item2)<=max_allow):
                    list1=item1;
                    list2=item2;
                    return True
    return False

#differences between two strings; hamming distance normalized to length of string.
def hamming_distance(str1, str2):
    assert len(str1) == len(str2)
    dist=0
    for i in range(0,len(str1)):
        if str1[i] != str2[i]:
            dist=dist+1
    dist=float(dist)/len(str1)
    return dist

def primer_search(sequence, primer):
    #threshold for acceptance; similarity metric between primer and sequence found will 
    #be between 0 (totally similar) and 1 (totally dissimilar) (lower values --> more simlar).
    #For the given set up, the number should not be above 0.5, since the first half must match
    #by the way we have decided to pull out "primer-like" sequences.
    # #testing 
    #print primer
    half_primer = primer[:int(math.ceil(.5*len(primer)))] #first half of the primer (index zero to floor(len/2)).
    
    #Regular expression for matching a sequence of length(primer)
    # with the first bit matching EXACTLY, and the second bit being anything.
    regex_primer = half_primer+'.{' + str(len(primer)-len(half_primer))+'}'
    #assert regex_primer='GGCATATTGC.{10}'
    #print regex_primer 
    
    #find instances of the regex primer in the full genome; gives lists of potential matches 
    regex_match_list = re.findall(regex_primer,sequence)
    #print "PRIMER: " + primer
    #print regex_match_list
    #for item in regex_match_list:
    #    print sequence.find(item)
    
    #optional print statements;... probably remove them
    #print('Primer: ' + primer)
    #print('Match list:')
    #print(regex_match_list);
    
    #if length of list is 0, then no matches are found, 
    # exit function bit early; no point in wasting CPU time when there is no match.
    if len(regex_match_list) < 1:
         return -1
    
    #at least one match has been found if we've gotten here. 
    best_match_ind = -1;#index of best match; defaults to -1.
    best_match_score = 1;#best match score (defaults to 1, worst score possible)
    i=0; # current index in list
    for item in regex_match_list: 
        #print primer
        #print item
        distance = hamming_distance(item,primer) #get distance between purported match and primer
        if distance < best_match_score: #check if new distance in better than best stored score
    		#store new best distance and index of occurrence if better
    		    best_match_score=distance
    		    best_match_ind=i
        i=i+1
    #print best_match_score
    
    #check if the best match found is below the threshold:
    if best_match_score < threshold:
    	#match found! :D
        #print 'Match found'

        #Note: best match may occur more than once; we shall return a list of indices instead of just a single index.
        genome_start_ind=[]
        ind=0
        end=len(sequence)
        while ind!= -1 and ind < end - len(primer) :
            curr = sequence.find(regex_match_list[best_match_ind],ind+len(primer))
            if curr!=-1:
                genome_start_ind.append(curr)
            else:
                break
            ind=curr
        #print genome_start_ind
    	#do whatever
        return genome_start_ind
    else:
    	#closest match found was not good enough. :((
     	return -1

def get_element(my_list, position):
    return my_list[position]


#Read in file
df= pd.read_excel(sys.argv[start])    
count = len(sys.argv)
rows = len(df.index)
start = start+1

for x in range(start, count):

    inFile=sys.argv[x]
#    x = x + 1

    with open(inFile, 'r') as myfile:
        inFile=myfile.read().replace('\n', '')
            
    sequence = inFile
    accession = sequence.partition(' ')[0]    
    sequence=sequence.upper()
    
    #Loop over rows in dataframe
    
    for y in range(len(df)):
        #get the name of the gene
        geneName = df.iloc[y,0]
        print(">"+str(geneName)+"_"+accession[1:])
        # get F primer for said gene
    #    print 'primerF'
        primerF = df.iloc[y,1]
#        print primerF
        #get R primer for said gene
    #    print primerF
        primerR = df.iloc[y, 2]
#        print 'PRIMER R'
#        print primerR
    
    
    
    
        minLength = int(str(df.iloc[y, 3]))
#        print minLength
        maxLength = int(str(df.iloc[y, 4]))
#        print maxLength


            
            
        
        
        # determine acceptable pcr fragment size
#        if bool(str(minLength)=='' or str(maxLength)=='')==True:
#            minLength = int(200)
#            maxLength = int(2500)
#        else:
#            minLength=int(df.iloc[y, 3])
#            maxLength=int(df.iloc[y, 4])
        

    
    #    print minLength
    #    print 'maxLength:'
    #    print maxLength
        
        #Format for Biopython
        
        #primerF_seq=Seq('"'+primerF+'"', generic_dna)
        #primerR_seq=Seq('"'+primerR+'"', generic_dna)
        primerF_seq=Seq('"'+primerF+'"')
        primerR_seq=Seq('"'+primerR+'"')

    #    print primerF_seq
    #    print primerR_seq
        #Reverse complement primers
        primerFRT = (str(primerF_seq.reverse_complement()))
        primerRRT = (str(primerR_seq.reverse_complement()))
    #    print 'primerFRT'
    #    print primerFRT
    #    print 'primerRRT'
    #    print primerRRT
        
        #Remove biopython format
        primerF = (str(primerF.replace('"',''))).upper()
        primerR = (str(primerR.replace('"',''))).upper()
        primerFRT = (str(primerFRT.replace('"',''))).upper()
        primerRRT = (str(primerRRT.replace('"',''))).upper()
        
        #Search for all primers
        temp1 = primer_search(sequence, primerF)
    #    print 'temp1'
        #print temp1  ####
        temp2 = primer_search(sequence, primerR)
    #    print 'temp2'
       # print temp2 ####
        temp3 = primer_search(sequence, primerFRT)
    #    print 'temp3'
      #  print temp3 ####
        temp4 = primer_search(sequence, primerRRT)
    #    print 'temp4 '
      #  print temp4 ####
        
        #Find primers with the band size specified
        if abs_list_search(temp1,temp4,minLength,maxLength):
            temp1=temp1[0]
            temp4=temp4[0]
#            print temp1
#            print temp2
            #(geneSeq=sequence[temp1:temp2] )if temp1 < temp2 else geneSeq=sequence[temp2:temp1];
            if temp1 < temp4:
                #geneSeq = sequence.split(sequence[temp1])[1].split(sequence[temp2])[0]
                geneSeq=sequence[temp1:temp4]
            else: 
                #geneSeq = sequence.split(sequence[temp2])[1].split(sequence[temp1])[0]
                geneSeq=sequence[temp4:temp1]
            print(geneSeq)
        elif abs_list_search(temp3,temp2,minLength,maxLength):
            temp3=temp3[0]
            temp2=temp2[0]
#            print temp3
#            print temp4
            #geneSeq=sequence[temp3:temp4] if temp3 < temp4 else geneSeq=sequence[temp4:temp3];
            if temp3 < temp2:
                #geneSeq = sequence.split(sequence[temp3])[1].split(sequence[temp4])[0]
                geneSeq=sequence[temp3:temp2]
            else: 
                #geneSeq = sequence.split(sequence[temp4])[1].split(sequence[temp3])[0]
                geneSeq=sequence[temp2:temp3]
            print(geneSeq)
    
        print('\n')



    
