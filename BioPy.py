
'''
Justin Levins
jlevs800@gmail.com
https://github.com/JustinLevins


BioPy is a short script designed to be used as a command line tool for obtaining more information
about proteins. Currently, it is capable of:

Displaying residues
Creating a Ramachandran plot using matplotlib
Displaying BLAST results by sending a GET request to NCBI

2/3: Program and first three options created     
2/5: include hydrophobicity plot using index values from Kyte and Doolittle, exit option in menu

It will download any .pdb file to your current working directory. After it is downloaded, use 
"pdb<code>.ent" if you are running the script again with the PDB file already downloaded.
'''

from Bio.PDB import *
import re
import os
import sys
import urllib.request
import matplotlib.pyplot as plt
import webbrowser
import requests
import Nyx #Nyx.py
import numpy as np
import pandas as pd


def download(code): #fetch a .pdb file and place in current directory
    bio = PDBList()
    return bio.retrieve_pdb_file(code, pdir=os.getcwd(), file_format='pdb')

def prompt(): 

    while True:
        
        query = input("Welcome to BioPy. Do you have a pdb file? Y/n: ")
        if query == 'y' or query == 'Y':
            return 0
        elif query == 'n' or query == 'N':
            return 1
        else:
            print("Invalid option! Please enter Y/n:\n")


def check(code): #regex for code validity
    code = str(code)
    if re.match(r"[0-9A-Za-z]{4}", code): 
        return 0
    else:
        return 1


def menu(code):
    print("\n------HOMEPAGE------\n")
    while True:
        query = str(input("Enter an option:\nA.List residues\nB.Ramachandran plot\nC.BLAST search\nD.Hydropathy plot\nE.Exit\nEnter an option: "))
        if query == 'a' or query == 'A':
            residues(code)
            
        elif query == 'b' or query == 'B':
            Rama(code)
            
        elif query == 'c' or query == 'C':
            blast()
            
        elif query == 'd' or query == 'D':
            hydropathy(code)

        elif query == 'e' or query == 'E':
            print("Now exiting...\n")
            sys.exit(1)
        else:
            print("Invalid response")


def conversion(value):
    return value * 57.2958
    #function to convert radians to degrees...180/pi = 57.2958
    #As of right now, this isn't used because radians produce the same image for the ramachandran plot. However, a Ramachandran is typically constructed in degrees.


def blast(): #uses RESTful web service to obtain blast results

    while True:
        code = input("Enter the 4-digit alphanumeric code again: ")
        if check(code) == 0:
            break
        else:
            print("invalid code:")
            pass

    serviceLocation = "https://www.rcsb.org/pdb/rest/getBlastPDB2?structureId=" + str(code) + "&chainId=A&eCutOff=10.0&matrix=BLOSUM62&outputFormat=HTML"
    r = requests.get(serviceLocation)
    webbrowser.open(str(r.url)) 
    

def size(structure): #  # of residues in structure
    i = 0
    model = structure[0]
    for chain in model:
        for residue in chain:
            if residue.get_resname() != 'HOH':
                i += 1
    return i

def Rama(code):
    
    #PHI IS CO->CO, x-axis
    #PSI IS N->N, y-axis

    parser = PDBParser()
    structure = parser.get_structure('molecule', code) 
    x = []
    y = []
    
    model = structure[0]
    
    for chain in model:
        polypeptides = PPBuilder().build_peptides(chain)
        for peptide in polypeptides:
            for item in peptide.get_phi_psi_list():
                x.append(item[0]) #get_phi_psi_list() returns a tuple--these values must be separated
                y.append(item[1]) 
       
    '''
   Conceptually, what the above code is doing:
    model = structure[0]
    try:
        for chain in model:
            for i in range (1, size(structure)):
                
                #define torsion angle atoms
                a0 = structure[0][chain.get_id()][i]['N']#changed from 'A'
                ab = structure[0][chain.get_id()][i]['CA']
                a1 = structure[0][chain.get_id()][i]['C']
                a2 = structure[0][chain.get_id()][i+1]['N'] 
                a3 = structure[0][chain.get_id()][i+1]['CA']
                a4 = structure[0][chain.get_id()][i+1]['C']
                
                #calculate atomic vectors to be used in torsion angle
                vector0 = a0.get_vector()
                vectorb = ab.get_vector()
                vector1 = a1.get_vector()
                vector2 = a2.get_vector()
                vector3 = a3.get_vector()
                vector4 = a4.get_vector()

                #create a dihedral angle using vectors
                phi = calc_dihedral(vector1, vector2, vector3, vector4)
                x.append(conversion(phi)) # <- convert radians to degrees and append to x

                psi = calc_dihedral(vector0, vectorb, vector1, vector2)
                y.append(conversion(psi)) # <- convert radians to degrees and append to y
'''
    print("(Close the plot to continue)")
    plt.xlabel('Phi, Radians') 
    plt.ylabel('Psi, Radians')
    plt.title("Ramachandran Plot of " + os.path.relpath(code))
    plt.axhline(y=0)
    plt.axvline(x=0)
    plt.scatter(x, y, label='dihedral angle', color='k', s=10) #creation of the Ramachandran plot
    plt.legend()
    plt.show()
    

def residues(code): 
    print("(Omitting water molecules)")
    parser = PDBParser()
    structure = parser.get_structure('molecule', code)
    model = structure[0] 
    for residue in model.get_residues():
        if residue.get_resname() != 'HOH':
            print(residue) 
    count = size(structure)
    print("Total residues found:")
    print(count)


def movingAverage(working_array, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(working_array, window, 'same')

def hydropathy(id): #Kyte Doolittle scale--utilizes a moving average
    
    parser = PDBParser()
    structure = parser.get_structure('molecule', id) 
    i = 1
    x = []
    #moving average size=19(or 21) center=10(9w/indices) 
    model = structure[0]
    window = []


    #fetch data    
    for residue in model.get_residues():
        code = residue.get_resname()
        if code != 'HOH' and code in Nyx.hydropathy:
            x.append(i)
            i += 1
            window.append(Nyx.hydropathy[code])

    avg = movingAverage(window, 19)
    print("(Close the plot to continue)") #more aesthetic
    plt.axhline(y=0)
    plt.xlabel("Residue Position")
    plt.ylabel("Hydropathy Index")
    plt.title("Hydropathy Plot for User Entered Structure")
    plt.plot(x, avg)
    plt.show()
    

def main():

    if prompt() == 1:
        fetch = input('Would you like to download a pdb file? Y/n: ')
        if fetch == 'y' or fetch == 'Y':
            molecu = input("Enter the 4 charater code: ")
            if check(molecu) == 0:
                downloaded = download(molecu)
                menu(downloaded)
            else:
                print("Error in reading the code. Please try again.\n")
        elif fetch == 'n' or fetch == 'N':
            print("Program will now exit.\n")
            sys.exit(1)
    else: 

        old_file = input("Move the file to the current directory and enter its filename (pdb files are stored as 'pdb<code>.ent'): ")
        menu(old_file)
        sys.exit(2)


if __name__ == "__main__":
    main()
