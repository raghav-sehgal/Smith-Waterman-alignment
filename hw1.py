#!/usr/bin/python
__author__ = "Raghav Sehgal"
__email__ = "raghav.sehgal@yale.edu"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt

import argparse
import numpy as np
import math
import pandas as pd

#### Passing arguments
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False,
default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False,
default=-1)
args = parser.parse_args()

def readfile(file_name):
   """ Function to read input text file. Expects two sequences in the first two lines"""
   with open(file_name, "r") as fin:
       data =fin.read().splitlines()
   seq1 = data[0]
   seq2 = data[1]
   return seq1, seq2

def Score(char1,char2,simMatrix):
   """Scoring function to get mismatch and match scores from scoring matrix"""
   return simMatrix.loc[char1][char2]
   
def Backtrace(matrix, seq1, seq2, openGap, extGap, mi, mj,simMatrix):
   """Backtrace algorithm"""
   align1 = ''
   align2 = ''
   current_score = matrix[mi][mj]
   while mi >= 1 and mj >= 1 and current_score>0:
       """Loop till scores greater than zero and all paths have been travered"""
       current_score = matrix[mi][mj]
       matchScore = simMatrix.loc[seq1[mi-1]][seq2[mj-1]]
       if current_score == matrix[mi-1][mj-1] + matchScore:
           align1 = align1+ seq1[mi-1]
           align2 = align2+ seq2[mj-1]
           mi -=1
           mj -=1
       elif (current_score == matrix[mi][mj-1] + extGap) or (current_score == matrix[mi][mj-1] + openGap):
           align2 = align2 + seq2[mj-1]
           align1 = align1 + "-"
           mj -=1
       elif (current_score == matrix[mi-1][mj] + extGap) or (current_score == matrix[mi-1][mj] + openGap):
           align1 = align1 + seq1[mi-1]
           align2 = align2 + "-"
           mi -=1
   align1 = align1[::-1]
   align2 = align2[::-1]
   return align1, align2

def Alignment(seq1,seq2,openGap, extGap,simMatrix):
   """Function which does scoring followed by calling backtrace and finally cleans up data for file writing"""
   row = len(seq1)+1
   col = len(seq2)+1
   #### Making Scoring matrices
   X = np.zeros((row,col), dtype=int)
   Y = np.zeros((row,col), dtype=int)
   m = np.zeros((row,col), dtype=int)
   M = np.zeros((row,col), dtype=int)
   score = 0
   mi, mj = (0,0)
   #### Scoring the matrices
   for i in range(1,row):
       for j in range(1,col):
           X[i][j] = max(X[i][j-1]+extGap,M[i][j-1]+openGap)
           Y[i][j] = max(Y[i-1][j]+extGap,M[i-1][j]+openGap)
           m[i][j] = M[i-1][j-1]+ Score(seq1[i-1],seq2[j-1],simMatrix)
           M[i][j] = max(0,m[i][j],X[i][j],Y[i][j])
           if score <= M[i][j]:
               score = M[i][j]
               mi,mj = i,j
   print(score)
   
   #### Calling traceback
   align = Backtrace(M, seq1, seq2, openGap, extGap, mi, mj, simMatrix)
   
   #### Traceback returns the exactly aligned sequence, following steps are done to make it similar to the results in the file
   #### Step 1 Add bracket at end and any excess sequence at the end of the original string missing in the aligned string
   align1 = align[0]+")"+seq1[mi:row]
   align2 = align[1]+")"+seq2[mj:col]
   
   #### Step 2 Add Bracket to start of alignment sequence and add the starting parts of the sequence which were not aligned
   #### Dum align is used to count all the characters which were aligned
   dum_align1 = align1.replace("-","")
   dum_align2 = align2.replace("-","")
   align1 = seq1[0:(row-len(dum_align1))] + "(" + align1
   align2 = seq2[0:(col-len(dum_align2))] + "(" + align2
   
   #### Step3 Add Spaces before and after the shorter string so that its visually similar in our data
   if(len(align1)>len(align2)):
       start_space = align1.find('(') - align2.find('(')
       end_space = len(align1) - align1.find(')')
       front_space = " "*start_space
       back_space =" "*(end_space-1)
       align2 = front_space + align2 + back_space
   elif(len(align2)>len(align1)):
       start_space = align2.find('(') - align1.find('(')
       end_space = len(align2) - align2.find(')') 
       front_space = " "*start_space
       back_space =" "*(end_space-1)
       align1 = front_space + align1 + back_space
   
   #### Step4 The exact match seuence is printed to show all the letters which had exact matches
   match = ""
   if(len(align1)==len(align2)):
       for k in range(len(align1)):
           if(align1[k]==align2[k] and align1[k].isalpha()):
               match += "|"
           else:
               match+=" "
   
   print("align1:{}".format(align1))
   print("match:{}".format(match))
   print("align2:{}".format(align2))
   
   #### Step 5 Making the scoring matrix easily readable
   M = np.vstack([np.array(list(" "+seq2)),M])
   M = np.vstack([np.array(list("  "+seq1)),M.T])
   
   return seq1,seq2, M,score, align1, match, align2

def writefile(file_name,result):
   """ Function to write to an external file"""
   #### Open file
   f = open(file_name, "w")

   #### Add seuqnces
   f.write("-----------\n")
   f.write("|Sequences|\n")
   f.write("-----------\n")
   f.write("sequence1\n")
   f.write(result[0]+"\n")
   f.write("sequence2\n")
   f.write(result[1]+"\n")
   
   #### Add Scoring matrix to file by looping over it
   f.write("--------------\n")
   f.write("|Score Matrix|\n")
   f.write("--------------\n")
   for i in range(len(result[2])):
       for j in range(len(result[2][i])):
           f.write(result[2][i,j]+"\t")
       f.write("\n")
    
   #### Add alignment scores and final alignment
   f.write("----------------------\n")
   f.write("|Best Local Alignment|\n")
   f.write("----------------------\n")
   f.write("Alignment Score:" + result[3].astype(str)+"\n")
   f.write("Alignment Results:\n")
   f.write(result[4]+"\n")
   f.write(result[5]+"\n")
   f.write(result[6]+"\n")
   f.close()
   return


def runSW(inputFile, scoreFile, openGap, extGap):
    #### Read input file
    pattern = readfile(inputFile)
    
    #### Read scoring file
    simMatrix = pd.read_fwf(scoreFile, index_col =0)
    
    #### Call Allignment
    result = Alignment(pattern[0],pattern[1],openGap, extGap, simMatrix)
    
    #### Write output file
    outputFile = inputFile.replace(".txt","")
    outputFile = outputFile +"_output.txt"
    writefile(outputFile, result)
    return result

runSW(args.input, args.score, args.opengap, args.extgap)
