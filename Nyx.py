'''---------------------HYDROPATHY INDICES---------------------
Kyte Doolittle scale
Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. J Mol Biol. (1982) 157:105–32. 10.1016/0022-2836(82)90515-0
 A higer hydropathy index indicates higher hydrophobicity.


This space can also be used for indices of other operations


For a degree of secondary structure prediction:
Koehl P, Levitt M. Structure-based conformational preferences of amino acids. Proc Natl Acad Sci U S A. 1999 Oct 26; 96(22): 12524–12529.
a-helix and b-sheet propensity values will also be used here.
'''

# uses a dynamic 'window' to calculate a moving average utilizing the defined values
hydropathy = {'ALA': 1.8, 'ARG': -4.5, 'ASN': -3.5, 'ASP': -3.5, 'CYS': 2.5, 'GLU': -3.5, 'GLN': -3.5, 'GLY': -0.4,
'HIS': -3.2, 'ILE': 4.5, 'LEU': 3.8, 'LYS': -3.9, 'MET': 1.9, 'PHE': 2.8, 'PRO': -1.6, 'SER': -0.8, 'THR': -0.7, 'TRP': -0.9, 'TYR': -1.3, 'VAL': 4.2}


