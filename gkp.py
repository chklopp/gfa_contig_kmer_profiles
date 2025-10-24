#Copyright (c) 2022-2025, INRAE - MIAT
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
#
#

import argparse
import sys
import subprocess
import textwrap

intervals = [[1,2],[2,3],[3,4],[4,11],[11,21],[21,101],[101,1000],[1001,10001],[10001,100001],[100001,100000000000]]

def load_gfa_file(gfa):
    """open a gfa file and produces the fasta and length coverage file, plus returns 
    length and coverage dictionnaries

    Parameters:
    gfa file name

    Returns:
    contig length and coverage dictionnary extracted from the gfa file
    """
    fout = open(gfa+".fasta","w")
    lcout = open(gfa+".length_coverage","w")
    coverage = {}
    length = {}
    
    with open(gfa,"r") as f :
        for l in f :
            if l.startswith("S") :
                b = l[:-1].replace("LN:i:","").replace("rd:i:","").split("\t")
                lcout.write(b[1]+"\t"+b[3]+"\t"+b[4]+"\n")
                length[b[1]] = int(b[3])
                coverage[b[1]] = int(b[4])
                fout.write(">"+b[1]+"\n")
                fout.write(textwrap.fill(b[2], 80)+"\n")
            
    return length, coverage


def check_fasta_length_coverage(fasta, length_coverage):
    """open a fasta file and a length_coverage file and checks the file lengths and 
    contig names found in the file as well and the contig lengths from both sources

    Parameters:
    fasta file
    length_coverage file 

    Returns:
    contig length and coverage dictionnary extracted from the length_coverage file
    """    
    fasta_d = {}
    length = {}
    coverage = {}
    name = ""
    # load fasta contigs
    with open(fasta,"r") as f :
       for line in f :
           if line[:1] == ">" :
               name = line[1:-1].split(" ")[0]
               #print(name)
               fasta_d[name]=""
           else :
               fasta_d[name]+=line[:-1]
    f.close()
    # load length and coverage
    with open(length_coverage,"r") as f :
       for line in f :
           length[line[:-1].split("\t")[0]] = int(line[:-1].split("\t")[1])
           coverage[line[:-1].split("\t")[0]] = int(line[:-1].split("\t")[2])
    f.close()
    # check contig names 
    #print(len(fasta_d) ,len(length))
    if len(fasta_d) != len(length) :
        sys.exit("Fasta and length_coverage file have not the same length")
    for k,v in fasta_d.items() :
        if len(v) != length[k] :
            print(len(v), length[k])
            sys.exit("Sequence "+k+"has not the same length in fasta and length_coverage file")       
           
    return length,coverage
            
def run_jellyfish(fasta, threads):
    """run jellyfish commands : kmer dictionnary building and query_by_sequence to count kmer from 
    the dictionnary in the contig fasta file 

    Parameters:
    fasta file
    number of theads used for jellyfish kmer database construction  

    Returns:
    name of the kmer count file .gps
    """
    try:
        print("generating kmer database...this can take some time")
        command = "jellyfish count -s 10000000 -m 21 -C -t"+str(threads)+" -o "+fasta+".jf "+fasta
        print(command)
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if result.returncode == 0 :
            print("database generated")
        else:
            print("database not generated")
    except FileNotFoundError:
        print("Could not be generate kmer database")
    try:
        print("producing kmer counts...this can take some time")
        command = "query_per_sequence "+fasta+".jf "+fasta+" > "+fasta+".gps_out"
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        if result.returncode == 0 :
            print("kmer counts produced")
        else:
            print("kmer counts not produced")
    except FileNotFoundError:
        print("Could not be generate kmer counts")
        
    return fasta+".gps_out"
                      
# Function to check if all specified software can be run on the system
def check_softwares(softwares):
    """verifies that the expected software packages are present in the environment 

    Parameters:
    list of software packages names
    
    Returns:
    None
    """    
    for software in softwares:
        #print(software)
        try:
            # Run the command with no arguments
            result = subprocess.run([software], stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
            # print(result.returncode)       
            # Check the return code
            if result.returncode == 0 or result.returncode == 1 :
                print("The "+software+" command is accessible.")
            else:
                print("The "+software+" command is not accessible.")
                sys.exit(1) 
                return False
        except FileNotFoundError:
            print("The "+software+" command is not found in the environment.")
            return False
        except Exception as e:
            print(f"An error occurred: {e}")
            return False

def format_line(output, l, coverage):
    """formats the count line  

    Parameters:
    vector of count values
    contig length
    dictionnary of contig coverages 
    
    Returns:
    formated count vector 
    """    
    b = [0,0,0,0,0,0,0,0,0,0]
    for ll in l :
        for j,i in enumerate(intervals) :
            if int(ll) >= i[0] and int(ll) < i[1] :
                b[j] +=1
        som = 0
        for i in b :
            som+=i
        for j,i in enumerate(b):
            b[j]=float(i*1000)/som
    return output+"\t"+str(coverage[output])+"\t"+str(len(l))+"\t"+"\t".join(map(str,b))+"\n"

output = ""
cmpt = 0

def process_counts(counts, coverage):
    """processes count file to extract vector of count values   

    Parameters:
    Jellyfish count file names
    dictionnary of contig coverages 
    
    Returns:
    None 
    """    
    output = ""
    lin = ""
    fout = open(counts+".permils", "w")
    with open(counts,"r") as f :
        for cmpt,line in enumerate(f) :
            block = line[:-1].split(" ")
            name = ""
            if line[0:1] == ">" :
                if cmpt != 0 :
                    fout.write(format_line(output,lin, coverage))
                output = block[0][1:]
            else :
                lin = block
            #print(block)
    fout.write(format_line(output,lin,coverage))

            
def main():
    parser = argparse.ArgumentParser(description="generate kmer profile.")
    parser.add_argument("--assembly_fasta", type=str, required=False, help="input assembly fasta file")
    parser.add_argument("--assembly_length_coverage", type=str, required=False, help="input assembly length and coverage file")
    parser.add_argument("--assembly_gfa", type=str, required=False, help="input assembly gfa file")
    parser.add_argument("--output", type=str, default="gkp_output.tsv", required=False, help="output file name")
    parser.add_argument("--mpthreads", type=int, default=4, required=False, help="threads to run jellyfish")

    args = parser.parse_args()
    fasta_file = ""
    
    # check software 
    softwares = ["jellyfish -h"]
    check_softwares(softwares)
    
    # check configuration 
    run = 0
    if args.assembly_fasta and args.assembly_length_coverage :
        run = 1
    if args.assembly_gfa :
        run = 2

    if run == 0 :
        print("Check your inputs, you have to provide fasta + length_coverage or gfa file")
    else :    
        if run == 1 :
            print("Checking fasta and length coverage file content...this can take a while")
            # load fasta file  
            (flength, fcoverage) = check_fasta_length_coverage(args.assembly_fasta, args.assembly_length_coverage)
            fasta_file = args.assembly_fasta
            
        if run == 2 :
            print("Producing fasta from gfa file...this can take a while")
            (flength, fcoverage) = load_gfa_file(args.assembly_gfa)
            fasta_file = args.assembly_gfa+".fasta"

        print("Running Jellyfish count...this can take a while")
        counts = run_jellyfish(fasta_file, args.mpthreads)
        print("Running query_per_sequence...this can take a while")
        fcounts = process_counts(counts, fcoverage)

if __name__ == "__main__":
    main()
