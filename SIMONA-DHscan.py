#!/usr/bin/env python
import os
import sys
import glob
import tarfile
import csv
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import subprocess
import yaml
import xml.etree.ElementTree as ET
import pandas as pd
#PYMOL
import __main__
__main__.pymol_argv = [ 'pymol', '-qc' ]
import pymol
pymol.finish_launching()
from pymol import cmd


def SMILE_Preprocessor(MOL, SmileCode): 
    """
    Process 
    -------------------------
    1. create PDB from smile code
    2. fix PDB file with right labels
    3. get Parameters with Acpype

    INPUT
    -------------------------
    1. MoleculeName = User input, name of molecule
    2. SmileCode = Code to create all-model 
    3. Charge = Net charge of system

    OUTPUT
    -------------------------
    1. Generates "Moleculename".acpype parameter directory
    2. Creates SMILES code into PDB and input for acpype
    """
    #Create Smilecode file to use later.
    with open("{}smile.txt".format(MOL), 'w') as outfile:
        outfile.write(SmileCode)

    #subprocess.run(["obabel", "-:{}".format(SmileCode), "-O", "{}tmp.pdb".format(MOL), "--gen3d", "--ff","GAFF"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    subprocess.run(["obabel", "-:{}".format(SmileCode), "-O", "{}tmp.pdb".format(MOL), "--gen3d", "--ff","GAFF"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #os.system("obabel -:" + SmileCode + " -O " + MOL + "tmp.pdb" + " --gen3d --ff GAFF")
 
    
    AtomsDictionary = {}

    with open('{}tmp.pdb'.format(MOL), 'r') as infile:
        LabelCoordinate = []
        Labels = []
        for line in infile:
            if 'HETATM' in line:
                column = line.split()
                Label = column[2]
                coordinate = line[20:80]
                Labels.append(Label)
                LabelCoordinate.append([Label, coordinate])

        Elements = list(sorted(set(Labels)))

        print('PDB has {} elements'.format(Elements))
        #print(LabelCoordinate)
        for element in Elements:
            ElementList = []
            for idx in range(0, len(LabelCoordinate)):
                    
                if element in LabelCoordinate[idx]:
                    #print(LabelCoordinate[idx])
                    ElementList.append(LabelCoordinate[idx])
            #print(ElementList)
            for id, info in enumerate(ElementList):
                #print(id, info)
                AtomsDictionary[str(info[1])] = '{}{}'.format(info[0], str(id+1))
    #print(AtomsDictionary)

    with open('{}.pdb'.format(MOL), 'w') as outfile:
        Coords = []
        with open('{}tmp.pdb'.format(MOL), 'r') as infile:        
            for line in infile:
                if 'HETATM' in line:
                    column = line.split()
                    coordinate = line[20:80]
                    Coords.append(coordinate)
                    
        #Write new PDB file
        outfile.write("REMARK temporal PDB for parameterization \n")
        outfile.write("AUTHOR SIMSTACK WANO SIMONA-DHSCAN " + str(date.today()) + " \n")
        for ai in range(0, len(Coords)):
            if ai < 9:
                if int(AtomsDictionary.get(Coords[ai])[1:]) < 10:
                    PDBline = "ATOM      " + str(ai+1) + "  " +  AtomsDictionary.get(Coords[ai]) +  "  " + MOL  + Coords[ai] + " \n"
                else:
                    PDBline = "ATOM      " + str(ai+1) + "  " +  AtomsDictionary.get(Coords[ai]) +  " " + MOL  + Coords[ai] + " \n"


            elif ai >= 99:
                if int(AtomsDictionary.get(Coords[ai])[1:]) < 10:
                    PDBline = "ATOM    " + str(ai+1) + "  " +  AtomsDictionary.get(Coords[ai]) +  "  "+ MOL + Coords[ai] + " \n"
                else:
                    PDBline = "ATOM    " + str(ai+1) + "  " +  AtomsDictionary.get(Coords[ai]) +  " "+ MOL + Coords[ai] + " \n"                

            else:
                if int(AtomsDictionary.get(Coords[ai])[1:]) < 10:
                    PDBline = "ATOM     " + str(ai+1) + "  " +  AtomsDictionary.get(Coords[ai]) +  "  "+ MOL + Coords[ai] + " \n"
                else:
                    PDBline = "ATOM     " + str(ai+1) + "  " +  AtomsDictionary.get(Coords[ai]) +  " "+ MOL + Coords[ai] + " \n"
               
            #print(PDBline)
            outfile.write(PDBline)
        outfile.write("END")
        os.system("rm {}tmp.pdb".format(MOL)) 



def AcPyPeProcess(MOL, charge=0):
    """
    Process 
    -------------------------
    Runs ACPYPE and it makes sure it finish!

    INPUT
    -------------------------
    1. MOl and charge

    OUTPUT
    -------------------------
    1. Creates MOL.acpype folder with all parameter files needed.
    """
       
    #ACPYPE Process TODO: solved!
    acpypeProcess = subprocess.Popen(["acpype", "-i", MOL + ".pdb", "-n", str(charge), "-o", "gmx"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")

    out, err = acpypeProcess.communicate()
    #print('output MONTSE ACPYPE {0}'.format(out))
    acpype_out = []
    acpype_out.append(out)
    for line in acpype_out:
        if "Total time of execution:" in line:
            print("ACPYPE IS FINISHED! :D")
            acpypeProcess.terminate()
        


def SIMONA_files(MOL, SIMONAPATHS): #TODO: workinglocally, not yet in simstack
    """
    Process 
    -------------------------
    This function depends on SMILE_Preprocessor. SIMONA Pre-processor.

    INPUT
    -------------------------
    1. It will look for Acpype Working Directory

    OUTPUT
    -------------------------
    1. Creates SML, custom_radii.itp pre-processor inputs for SIMONA pre-processor
    2. Creates epqr, epqr.dh, mol2, spf and xml files for SIMONA
    3. Check if atomtypes in SIMONA are correct respecting GMX atomtypes
    """
    
    AcpypeWorkingDir = glob.glob("*acpype")
    #print(AcpypeWorkingDir)
    os.mkdir("SIMONA_inputs")
    os.chdir('SIMONA_inputs')
    #print("Current working directory: {0}".format(os.getcwd()))
    #os.system("cp ../{}/*GMX.* .".format(AcpypeWorkingDir))
    os.system("cp ../{}/*_GMX.* .".format(AcpypeWorkingDir[0]))
    

    #Create the custom_radii.itp from library
    #dictionary for atom radii
    Radii_Lib = {
    'c': '0.170',
    'c1': '0.170',
    'c2': '0.170',
    'c3': '0.170',
    'ca': '0.170',
    'cp': '0.170',
    'cq': '0.170',
    'cc': '0.170',
    'cd': '0.170',
    'ce': '0.170',
    'cf': '0.170',
    'cg': '0.170',
    'ch': '0.170',
    'cx': '0.170',
    'cy': '0.170',
    'cu': '0.170',
    'cv': '0.170',
    'cz': '0.170',
    'h1': '0.120',
    'h2': '0.120',
    'h3': '0.120',
    'h4': '0.120',
    'h5': '0.120',
    'ha': '0.120',
    'hc': '0.120',
    'hn': '0.120',
    'ho': '0.120',
    'hp': '0.120',
    'hs': '0.120',
    'hw': '0.120',
    'hx': '0.120',
    'f': '0.147',
    'cl': '0.175',
    'br': '0.185',
    'i': '0.198',
    'n': '0.155',
    'n1': '0.155',
    'n2': '0.155',
    'n3': '0.155',
    'n4': '0.155',
    'na': '0.155',
    'nb': '0.155',
    'nc': '0.155',
    'nd': '0.155',
    'ne': '0.155',
    'nf': '0.155',
    'nh': '0.155',
    'no': '0.155',
    'ni': '0.155',
    'nj': '0.155',
    'nk': '0.155',
    'nl': '0.155',
    'nm': '0.155',
    'nn': '0.155',
    'np': '0.155',
    'nq': '0.155',
    'o': '0.152',
    'oh': '0.152',
    'os': '0.152',
    'op': '0.152',
    'oq': '0.152',
    'ow': '0.152',
    'p2': '0.180',
    'p3': '0.180',
    'p4': '0.180',
    'p5': '0.180',
    'pb': '0.180',
    'pc': '0.180',
    'pd': '0.180',
    'pe': '0.180',
    'pf': '0.180',
    'px': '0.180',
    'py': '0.180',
    's': '0.180',
    's2': '0.180',
    's4': '0.180',
    's6': '0.180',
    'sh': '0.180',
    'ss': '0.180',
    'sp': '0.180',
    'sq': '0.180',
    'sx': '0.180',
    'sy': '0.180' }

    Element_Lib = {
    'c': 'C',
    'c1': 'C',
    'c2': 'C',
    'c3': 'C',
    'ca': 'C',
    'cp': 'C',
    'cq': 'C',
    'cc': 'C',
    'cd': 'C',
    'ce': 'C',
    'cf': 'C',
    'cg': 'C',
    'ch': 'C',
    'cx': 'C',
    'cy': 'C',
    'cu': 'C',
    'cv': 'C',
    'cz': 'C',
    'h1': 'H',
    'h2': 'H',
    'h3': 'H',
    'h4': 'H',
    'h5': 'H',
    'ha': 'H',
    'hc': 'H',
    'hn': 'H',
    'ho': 'H',
    'hp': 'H',
    'hs': 'H',
    'hw': 'H',
    'hx': 'H',
    'f': 'F',
    'cl': 'Cl',
    'br': 'Br',
    'i': 'I',
    'n': 'N',
    'n1': 'N',
    'n2': 'N',
    'n3': 'N',
    'n4': 'N',
    'na': 'N',
    'nb': 'N',
    'nc': 'N',
    'nd': 'N',
    'ne': 'N',
    'nf': 'N',
    'nh': 'N',
    'no': 'N',
    'ni': 'N',
    'nj': 'N',
    'nk': 'N',
    'nl': 'N',
    'nm': 'N',
    'nn': 'N',
    'np': 'N',
    'nq': 'N',
    'o': 'O',
    'oh': 'O',
    'os': 'O',
    'op': 'O',
    'oq': 'O',
    'ow': 'O',
    'p2': 'P',
    'p3': 'P',
    'p4': 'P',
    'p5': 'P',
    'pb': 'P',
    'pc': 'P',
    'pd': 'P',
    'pe': 'P',
    'pf': 'P',
    'px': 'P',
    'py': 'P',
    's': 'S',
    's2': 'S',
    's4': 'S',
    's6': 'S',
    'sh': 'S',
    'ss': 'S',
    'sp': 'S',
    'sq': 'S',
    'sx': 'S',
    'sy': 'S' }

    AtomTypesDic = {} # atomlabel : GMX_atomtype

    with open("custom_radii.itp", "w") as outfile:
        outfile.write("[ implicit_genborn_params ] \n#This table is taken from gromacs/share/gromacs/top/amber99sb-ildn.ff/gbsa.itp -->  (if you are not satisfied with these values, look for a better source for the radii (e.g. c3,os,hn was missing and added by hand with dummy value for everything except the radius)) \n# the following header might not be correct (copied from other .itp file) e.g. gbr is the radius \n; atype      sar      st     pi       gbr       hct \n;\n")
        #GET the atomtypes from GMX.itp
        Parmfile = "{}_GMX.itp".format(MOL)
        AllLines = []
        with open(Parmfile, 'r') as input:
            for line in input:
                AllLines.append(line)
        LineLimits = []
        for index, line in enumerate(AllLines):
            if "atoms" in line:
                #print(index, line)
                LineLimits.append(index)
            if "bonds" in line:
                #print(index, line)
                LineLimits.append(index)
        AtomTypes = []
        for line in AllLines[LineLimits[0]+2:LineLimits[1]-1]:
            line = line.split()
            atomtype = line[1]
            AtomTypes.append(atomtype)
            atomlabel = line[4]
            AtomTypesDic[str(atomlabel)] = str(atomtype)
            #print(atomtype, atomlabel)
        Types = list(sorted(set(AtomTypes)))
        print(Types)

        #Add the lines in custom_radii.itp
        for atom in Types:
            #print(atom, Radii_Lib.get(atom))
            if len(atom) == 2:
                #br           0.1      1      1        0.185   0.85 ; H
                outfile.write(atom + "           0.0      0      0        " + Radii_Lib.get(atom) + "   0.00 ;" + Element_Lib.get(atom)+ "\n")
            else: #n            0.155    1      1.028    0.155   0.79 ; N
                outfile.write(atom + "          0.0      0      0        " + Radii_Lib.get(atom) + "   0.00 ;" + Element_Lib.get(atom)+ "\n")


    # HEAD: 1-2, MOVES: 3- ?, if you need to modify something, just update the value of the key.
    sml_LIB = {
    '1': 'colorize: false', 
    '2': 'forcefield_spec:', 
    '3': 'moves:', 
    '4': '  analysis_moves:', 
    '5': '  - - print_energy', 
    '6': '    - begin_step: 0', 
    '7': '      last_step: 0', 
    '8': '      step_mod: 1000', 
    '9': '  - - BestConfigurationOutput', 
    '10': '    - begin_step: 0', 
    '11': '      data_type: pdb', 
    '12': '      fname: best.pdb', 
    '13': '      last_step: 0',
    '14': '      step_mod: 1', 
    '15': '  - - trajectory', 
    '16': '    - begin_step: 0', 
    '17': '      fname: trajectory.pdb', 
    '18': '      last_step: 0', 
    '19': '      only_new: 0', 
    '20': '      step_mod: 1000', 
    '21': '  - - energy', 
    '22': '    - begin_step: 0', 
    '23': '      last_step: 0', 
    '24': '      step_mod: 1000', 
    '25': '  initial: []', 
    '26': '  list:', 
    '27': '  - - new_dihedrals', 
    '28': '    - allow:', 
    '29': '      - all', 
    '30': '      angles: ' + repr('*') , 
    '31': '      delta_phi_max: 1.5235987755982988', 
    '32': 'nsteps: 100000', 
    '33': 'peptide_spec: false', 
    '34': 'preprocessor:', 
    '35': '  algorithm:', 
    '36': '    name: metropolis', 
    '37': '    params:', 
    '38': '      kB: 0.0019858775', 
    '39': '  atom_params: topol.spf', 
    '40': '  name: nano', 
    '41': '  simonaparser_use_bonds: false', 
    '42': '  simonapdbparser_auto_rename: false', 
    '43': '  simonapdbparser_connects: false', 
    '44': '  simonapdbparser_occ_as_charge: false', 
    '45': '  treat_unknown: delete', 
    '46': '  use_simona_pdb_parser: true', 
    '47': '  analysis_moves:', 
    '48': 'print_level: 1', 
    '49': 'seed: random', 
    '50': 'sourceFormat: 5', 
    '51': 'tend: 300.0', 
    '52': 'tstart: 300.0', 
    '53': 'verboseDihedral: false', 
    '54': 'warn_level: 5', 
    '55': 'xml_indent: true'}
    
    # FF : Nonbonded2 : LJ, Coulomb and GB.
    #FF_LIB = {'1': '- - BornRadii3', '2': '  - all_atoms: true', '3': '    scale: 1.0', '4': '- - Nonbonded2', '5': '  - depth: 3', '6': '    e_in: 4.0', '7': '    e_out: 80.0', '8': '    scale: 1.0'}

    # FF : Nonbonded2 : LJ, Coulomb, pit potential and GB.
    FF_LIB = {
    '1': '- !!python/tuple', 
    '2': '  - NonbondedVacuum', 
    '3' : '  - depth: 3', 
    '4': '    e_in: 1.0', 
    '5': '- - PitPotential', 
    '6': '  - scale: 1.0', 
    '7': '    steepness: 2000.0', 
    '8': '    subunits: '+ repr('-1'), 
    '9': '    x_max: 900.0', 
    '10': '    x_min: -30.0', 
    '11': '    y_max: 900.0', 
    '12': '    y_min: -30.0', 
    '13': '    z_max: 900.0', 
    '14': '    z_min: -30.0'}

    FF_terms = len(list(FF_LIB))

    #Print conf.sml according to FF terms. Customizable.
    with open('conf.sml', 'w') as outfile:
        for index in range(1, 3):
            outfile.write(sml_LIB.get(str(index))+ "\n")
        for index in range(1, FF_terms +1):
            outfile.write(FF_LIB.get(str(index))+ "\n")
        for index in range(3, 56):
            outfile.write(sml_LIB.get(str(index))+ "\n")

    #subprocess with SIMONA scripts

    MoleculeName = glob.glob("*_GMX.top")[0][:3]
    #gromacs_to_epqr.py *_GMX.top *_GMX.gro topol.epqr custom_radii.itp
    GRO2epqr = '{}gromacs_to_epqr.py'.format(SIMONAPATHS.get('PythonDir'))
    print(GRO2epqr)
    #GRO2epqr = 'gromacs_to_epqr.py'
    subprocess.run(['python', GRO2epqr, MoleculeName + "_GMX.top", MoleculeName + "_GMX.gro", "topol.epqr", "custom_radii.itp"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    #obabel -i pqr topol.epqr -o mol2 -O topol.mol2      
    subprocess.run(["obabel", "-i", "pqr", "topol.epqr", "-o", "mol2", "-O", "tmp.mol2"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    subprocess.run(["obabel", "-i", "pqr", "topol.epqr", "-o", "pdb", "-O", "topol.pdb"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    #we need to double check the atomtypes in MOL2 (SYBYL atom types) file since obabel guess the atomtyper for SIMONA and with C2 types gets bad.
    #TODO: SOME atom types MUST be double checked.
    Mol2_Lib = {
    'c': 'C.2',      # Sp2 C carbonyl group        
    'c1': 'C.1',      # Sp C
    'c2': 'C.2',      # Sp2 C  
    'c3': 'C.3',      # Sp3 C
    'ca': 'C.ar',     # Sp2 C in pure aromatic systems
    'cp': 'C.ar',     # Head Sp2 C that connect two rings in biphenyl sys. 
    'cq': 'C.ar',     # Head Sp2 C that connect two rings in biphenyl sys. identical to cp 
    'cc': 'C.2',      # Sp2 carbons in non-pure aromatic systems
    'cd': 'C.2',      # Sp2 carbons in non-pure aromatic systems, identical to cc
    'ce': 'C.2',      # Inner Sp2 carbons in conjugated systems
    'cf': 'C.2',      # Inner Sp2 carbons in conjugated systems, identical to ce
    'cg': 'C.1',      # Inner Sp carbons in conjugated systems
    'ch': 'C.1',      # Inner Sp carbons in conjugated systems, identical to cg
    'cx': 'C.3',      # Sp3 carbons in triangle systems
    'cy': 'C.3',      # Sp3 carbons in square systems
    'cu': 'C.2',      # Sp2 carbons in triangle systems
    'cv': 'C.2',      # Sp2 carbons in square systems
    'cz': 'C.2',      # Sp2 carbon in guanidine group
    'h1': 'H',        # H bonded to aliphatic carbon with 1 electrwd. group  
    'h2': 'H',        # H bonded to aliphatic carbon with 2 electrwd. group 
    'h3': 'H',        # H bonded to aliphatic carbon with 3 electrwd. group 
    'h4': 'H',        # H bonded to non-sp3 carbon with 1 electrwd. group 
    'h5': 'H',        # H bonded to non-sp3 carbon with 2 electrwd. group 
    'ha': 'H',        # H bonded to aromatic carbon  
    'hc': 'H',        # H bonded to aliphatic carbon without electrwd. group 
    'hn': 'H',        # H bonded to nitrogen atoms
    'ho': 'H',        # Hydroxyl group
    'hp': 'H',        # H bonded to phosphate 
    'hs': 'H',        # Hydrogen bonded to sulphur 
    'hw': 'H',        # Hydrogen in water 
    'hx': 'H',        # H bonded to C next to positively charged group  
    'f': 'F',        # Fluorine
    'cl': 'Cl',       # Chlorine 
    'br': 'Br',       # Bromine 
    'i': 'I',        # Iodine 
    'n': 'N.am',        # Sp2 nitrogen in amide groups
    'n1': 'N.1',      # Sp N  
    'n2': 'N.2',      # aliphatic Sp2 N with two connected atoms 
    'n3': 'N.3',      # Sp3 N with three connected atoms
    'n4': 'N.3',      # Sp3 N with four connected atoms 
    'na': 'N.2',      # Sp2 N with three connected atoms 
    'nb': 'N.ar',      # Sp2 N in pure aromatic systems 
    'nc': 'N.2',      # Sp2 N in non-pure aromatic systems
    'nd': 'N.2',      # Sp2 N in non-pure aromatic systems, identical to nc
    'ne': 'N.2',      # Inner Sp2 N in conjugated systems
    'nf': 'N.2',      # Inner Sp2 N in conjugated systems, identical to ne
    'nh': 'N.2',      # Amine N connected one or more aromatic rings (CHECK  again)
    'no': 'N.2',      # Nitro N  
    'ni': 'N.2',        # n in 3-memberred rings
    'nj': 'N.2',        # n in 4-memberred rings
    'nk': 'N.2',        # n4 in 3-memberred rings
    'nl': 'N.2',        # n4 in 4-memberred rings
    'nm': 'N.2',        # nh in 3-memberred rings
    'nn': 'N.2',        # nh in 4-memberred rings
    'np': 'N.2',        # n3 in 3-memberred rings
    'nq': 'N.2',        # n3 in 4-memberred rings
    'o': 'O.2',      # Oxygen with one connected atom
    'oh': 'O.3',      # Oxygen in hydroxyl group
    'os': 'O.3',        # Ether and ester oxygen
    'op': 'O.3',        # os in 3-memberred rings
    'oq': 'O.3',        # os in 4-memberred rings
    'ow': 'O.3',        # Oxygen in water 
    'p2': 'P.3',        # Phosphate with two connected atoms 
    'p3': 'P.3',        # Phosphate with three connected atoms, such as PH3
    'p4': 'P.3',        # Phosphate with three connected atoms, such as O=P(CH3)2
    'p5': 'P.3',        # Phosphate with four connected atoms, such as O=P(OH)3
    'pb': 'P.3',        # Sp2 P in pure aromatic systems 
    'pc': 'P.3',        # Sp2 P in non-pure aromatic systems
    'pd': 'P.3',        # Sp2 P in non-pure aromatic systems, identical to pc
    'pe': 'P.3',        # Inner Sp2 P in conjugated systems
    'pf': 'P.3',        # Inner Sp2 P in conjugated systems, identical to pe
    'px': 'P.3',        # Special p4 in conjugated systems
    'py': 'P.3',        # Special p5 in conjugated systems
    's': 'S.2',        # S with one connected atom 
    's2': 'S.3',        # S with two connected atom, involved at least one double bond  
    's4': 'S.O',        # S with three connected atoms 
    's6': 'S.O2',        # S with four connected atoms 
    'sh': 'S.3',        # Sp3 S connected with hydrogen 
    'ss': 'S.3',        # Sp3 S in thio-ester and thio-ether
    'sp': 'S.2',        # ss in 3-memberred rings
    'sq': 'S.2',        # ss in 4-memberred rings
    'sx': 'S.2',        # Special s4 in conjugated systems
    'sy': 'S.2'}        # Special s6 in conjugated systems

    #GMX_check = {} = AtomTypesDic 
    Mol2_check = {}
    AtomsList = []

    with open('topol.mol2','w') as outfile:
        
        #Get the atomtypes from tmp.mol2 file
        with open('tmp.mol2','r') as infile:
            AllLines = []
            LimitsLines = []
            LinesToCheck = {}
            for line in infile:
                AllLines.append(line)

            for index, line in enumerate(AllLines):
                if "@<TRIPOS>ATOM" in line:
                    LimitsLines.append(index)
                if "@<TRIPOS>BOND" in line:
                    LimitsLines.append(index)

            for lineID in range(LimitsLines[0]+1, LimitsLines[1]):
                #      1  C1         1.8800   -3.1200   -0.9700 C.3     1  RN11        1.7000
                line = AllLines[lineID]
                #print(lineID, line)
                column = line.split()
                atomname = column[1]
                AtomsList.append(atomname)
                mol2type = column[5]
                #print(atomname, mol2type)
                LinesToCheck[atomname] = line
                Mol2_check[atomname] = mol2type

            #compare the GMX atomtypes with the mol2 atomtypes
            for lineID in range(0, LimitsLines[0]+1):
                line = AllLines[lineID]
                outfile.write(line)
            with open('atomtype_CHECK.dat', 'w') as checkfile:
                for name in AtomsList:
                    AtomtypeGMX2Mol2 = Mol2_Lib.get(AtomTypesDic.get(name)) # is the right type from Mol2_Lib

                    if AtomtypeGMX2Mol2 != Mol2_check.get(name):
                    #correct the atomtype
                        outfile.write(LinesToCheck.get(name)[0:46] + AtomtypeGMX2Mol2 + LinesToCheck.get(name)[49:])
                        print("name {} had wrong atomtype".format(name))
                        checkfile.write("[REPLACED] {} : {} --> {} \n".format(name, Mol2_check.get(name), AtomtypeGMX2Mol2))
                    
                    else:
                    #keep the line
                        outfile.write(LinesToCheck.get(name))
                        checkfile.write("[OK] {} : {} \n".format(name, Mol2_check.get(name)))

            for lineID in range(LimitsLines[1], len(AllLines)):
                line = AllLines[lineID]
                outfile.write(line)

    MOL2SPF = '{}mol2spf.py'.format(SIMONAPATHS.get('PythonDir'))
    IdentRun = '{}SIMGenerateIdenticalRun.py'.format(SIMONAPATHS.get('PythonDir'))
    #MOL2SPF = 'mol2spf.py'
    #IdentRun ='SIMGenerateIdenticalRun.py'

    #mol2spf.py topol.mol2 topol.epqr topol.spf   
    subprocess.run([MOL2SPF, "topol.mol2", "topol.epqr", "topol.spf"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    #SIMGenerateIdenticalRun.py conf.sml topol.pdb tmp.xml
    subprocess.run([IdentRun, "conf.sml", "topol.pdb", "tmp.xml"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.chdir('../')
    return AtomTypesDic

def DihedralInfo(AtomTypesDic):
    # how many DH are in tmp.xml? which are the DH_ids?
    os.chdir('SIMONA_inputs')
    ids = []
    atomnames = []
    DH_atoms = []
    #DihedralInfoList = []
    atoms_ids = []
    tree = ET.parse('tmp.xml') 
    root = tree.getroot()
    #<coord chain_id="0" residue_id="0" id="0" name="C1" type="C" residue_name="MON">
    for x in root.iter('coord'):
        #print(x.attrib)
        atom_id = x.get('id')
        atom_name = x.get('name')
        #print(atom_name, atom_id)
        atomnames.append(atom_name)
        ids.append(atom_id)
    
    atoms_ids = dict(zip(atomnames, ids)) #Dictionary of atomname : atoms_id
    #print(atoms_ids)
    #<atomnames>H3.C1.C2.H5</atomnames>
    for x in root[1][0]:
        for ai in x:
            if ai.tag == "unique_id":
                #print(ai.tag)
                item = x.find('atomnames').text
                #print(item)
                atoms = item.split('.')
                DH_atoms.append(atoms)  

    # which are the aomt_ids dihedral_id per dihedral?
    DH_index = {}
    CurrentDihedrals_SIMONA = []  #simona xml format form id start in 0
    CurrentDihedrals_PDB = []     #PDB or xyz format form id start in 1
    with open('dihedral_details.txt', 'w') as outfile:
        #print('Total of available dihedrals : {}'.format(len(DH_atoms)))
        outfile.write('Total of available dihedrals : {} \n'.format(len(DH_atoms)))
        outfile.write('{} atom id : Atoms : Atom Types\n'.format(MOL))
        for dnumber, dihedral in enumerate(DH_atoms) :
            tmp_index = {}

            outfile.write('{} : {} : {} - {} - {} - {} \n'.format(dnumber, dihedral, AtomTypesDic.get(dihedral[0]), AtomTypesDic.get(dihedral[1]), AtomTypesDic.get(dihedral[2]), AtomTypesDic.get(dihedral[3])))
            for ai in dihedral:
                LabelIdPair = atoms_ids.get(str(ai))
                #print(ai , ":", int(LabelIdPair) +1) 
                tmp_index.update({str(ai): str(int(LabelIdPair) +1)})
                
                outfile.write('    - {} : {} \n'.format(ai , int(LabelIdPair) +1))
            DH_index.update({'{}_SIMONAScanStep'.format(dnumber): tmp_index})
    #print(DH_index)


    os.chdir('../')
    #print('CHECK')
    #print(DH_atoms)
    return DH_atoms, DH_index    # DihedralInfo_RETURNS

def XMLScanModification(DH_atoms, angle=0.174533, step=36): #TODO:working
    """
    Process 
    -------------------------
    This function depends on SIMONA_files. Creates new XML input files to run fast FF based DH scan calculations.
    Scan algorithm is save in a dictionary called ScanAlgorithmInput

    INPUT
    -------------------------
    1. tmp.xml file created in SIMONA_files. 

    OUTPUT
    -------------------------
    1. check the list of dihedrals recognized by SIMONA pre-processor
    2. Creates different XML to run different scans.
    """
    #I made a dict for the whole Algorithm protocol. Then the keys can be updated if some parameter need to be modify.
    ScanAlgorithmInput = {
    '1': '<Algorithm>', 
    '2': '<TransformationSequence repeats="36" weight="1.0">', 
    '3': '<TransformationSequence weight="1.0" repeats="1">', 
    '4': '<RecalcEnergies>',
    '5': '<beginstep>0</beginstep>', 
    '6': '<laststep>0</laststep>', 
    '7': '<stepmod>1</stepmod>', 
    '8': '</RecalcEnergies>', 
    '9': '<PrintEnergy>', 
    '10': '<beginstep>0</beginstep>', 
    '11': '<laststep>0</laststep>', 
    '12': '<stepmod>1</stepmod>', 
    '13': '</PrintEnergy>', 
    '14': '<BestConfigurationOutput>', 
    '15': '<beginstep>0</beginstep>', 
    '16': '<filename>best.pdb</filename>', 
    '17': '<laststep>0</laststep>',
    '18': '<stepmod>1</stepmod>', 
    '19': '<type>pdb</type>', 
    '20': '</BestConfigurationOutput>', 
    '21': '<ConfigurationOutput>', 
    '22': '<beginstep>0</beginstep>', 
    '23': '<filename>trajectory.pdb</filename>', 
    '24': '<laststep>0</laststep>', 
    '25': '<only_new>0</only_new>', 
    '26': '<stepmod>1</stepmod>', 
    '27': '<type>pdb</type>', 
    '28': '</ConfigurationOutput>', 
    '29': '<MetadataOutput>', 
    '30': '<beginstep>0</beginstep>', 
    '31': '<laststep>0</laststep>', 
    '32': '<metadataname>Summary</metadataname>', 
    '33': '<stepmod>1</stepmod>', 
    '34': '</MetadataOutput>', 
    '35': '</TransformationSequence>', 
    '36': '<SetDihedralRelative weight="1.0">', 
    '37': '<dihedral_id>2</dihedral_id>', 
    '38': '<value>0.174533</value>', 
    '39': '</SetDihedralRelative>', 
    '40': '</TransformationSequence>', 
    '41': '</Algorithm>'}

    #update angle, step
    ScanAlgorithmInput.update({'2': '<TransformationSequence repeats="{}" weight="1.0">'.format(str(step)), '38' : '<value>{}</value>'.format(str(angle))})

    #create a folder for each DH calculation with simona

    for dh_id in range(0, len(DH_atoms)):
        os.mkdir('{}_SIMONAScanStep'.format(dh_id))
        os.chdir('{}_SIMONAScanStep'.format(dh_id))

        with open('{}_dihedralscan.xml'.format(dh_id), 'w') as outfile:
            XMLLines = []
            AlgorithmLimits = []
            with open('../SIMONA_inputs/tmp.xml','r') as infile:
                for line in infile:
                    XMLLines.append(line)

                for index, line in enumerate(XMLLines):
                    if "<Algorithm>" in line:
                        AlgorithmLimits.append(index)
                    if "</Algorithm>" in line:
                        AlgorithmLimits.append(index)

                #write all lines from 0 to <Algorithm>
                for idx in range(0,AlgorithmLimits[0]): # not including Algorithm
                    outfile.write(XMLLines[idx])

                #replace with DH arbitrary dh scan lines + right dihedral id
                #update values in the dictionary per dh_id
                ScanAlgorithmInput.update({'37' : '<dihedral_id>{}</dihedral_id>'.format(dh_id)})

                for idx in range(1, len(ScanAlgorithmInput)+1):
                    outfile.write(ScanAlgorithmInput.get(str(idx))+ "\n")
                
                #print(AlgorithmLimits)
                for idx in range(AlgorithmLimits[1]+1,len(XMLLines)): # not including Algorithm
                    outfile.write(XMLLines[idx])

        os.chdir('../')



def XMLScanModification2(DH_atoms, angle=0.174533, step=36, sigma_0=3.14159, MC_sigma=0.523599, repeat=10000, t_0=300, t_i=300): #TODO:Running, test settings needed
    """
    Process 
    -------------------------
    This function depends on SIMONA_files. Creates new XML input files to run fast FF based DH scan calculations.
    Scan algorithm is save in a dictionary called ScanAlgorithmInput. In this function it is considered to relax sidechains after arbitrary rotation.

    INPUT
    -------------------------
    1. tmp.xml file created in SIMONA_files. 

    OUTPUT
    -------------------------
    1. check the list of dihedrals recognized by SIMONA pre-processor
    2. Creates different XML to run different scans.

    Algorithm
    ---------

        1. set dh in 180 with AbsoluteDHLib = sigma_0=3.14159
        2. Loop with rotation in steps and angle : angle=0.174533, step=36
        3. MC steps for side chain relaxation. It considers: MC_sigma=0.523599, repeat=10000, t_0=300, t_i=300
        4.print conformation and energy
    """

    Algorithm = {
    '1': '  <Algorithm>', 
    '2': '    <TransformationSequence repeats="1" weight="1.0">', 
    '3': 'AbsoluteTransformation', 
    '4': '    </TransformationSequence>', 
    '5': '    <TransformationSequence repeats="36" weight="1.0">', 
    '6': 'RelativeTransformation', 
    '7': '    <TransformationSequence repeats="1" weight="1.0">', 
    '8': '      <RepeatedMove>', 
    '9': '        <repeats>1000</repeats>', 
    '10': '        <tend>300.0</tend>', 
    '11': '        <tscaling>geometric</tscaling>', 
    '12': '        <tstart>300.0</tstart>', 
    '13': '        <TransformationSequence weight="1.0" repeats="1">', 
    '14': '          <TransformationSequence repeats="1" weight="1.0">', 
    '15': '            <RecalcEnergies>', 
    '16': '              <beginstep>0</beginstep>', 
    '17': '              <laststep>0</laststep>', 
    '18': '              <stepmod>1000</stepmod>', 
    '19': '            </RecalcEnergies>', 
    '20': '            <PrintEnergy>', 
    '21': '              <beginstep>0</beginstep>', 
    '22': '              <laststep>0</laststep>', 
    '23': '              <stepmod>1000</stepmod>', 
    '24': '            </PrintEnergy>', 
    '25': '            <BestConfigurationOutput>', 
    '26': '              <beginstep>0</beginstep>', 
    '27': '              <filename>best.pdb</filename>', 
    '28': '              <laststep>0</laststep>', 
    '29': '              <stepmod>1</stepmod>', 
    '30': '              <type>pdb</type>', 
    '31': '            </BestConfigurationOutput>', 
    '32': '            <ConfigurationOutput>', 
    '33': '                <beginstep>0</beginstep>', 
    '34': '                <filename>trajectory.pdb</filename>', 
    '35': '                <laststep>0</laststep>', 
    '36': '                <only_new>0</only_new>', 
    '37': '                <stepmod>1000</stepmod>', 
    '38': '                <type>pdb</type>', 
    '39': '            </ConfigurationOutput>', 
    '40': '            <MetadataOutput>', 
    '41': '                <beginstep>0</beginstep>', 
    '42': '                <laststep>0</laststep>', 
    '43': '                <metadataname>Summary</metadataname>', 
    '44': '                <stepmod>1000</stepmod>', 
    '45': '            </MetadataOutput>', 
    '46': '          </TransformationSequence>', 
    '47': '          <ConditionalTransformation weight="1.0">', 
    '48': '            <MetropolisAcceptanceCriterion>', 
    '49': '              <energymodel_nr>0</energymodel_nr>', 
    '50': '              <kB>0.0019858775</kB>', 
    '51': '            </MetropolisAcceptanceCriterion>', 
    '52': '            <TransformationChoice weight="1.0" repeats="1">', 
    '53': 'RelativeRandomTransformation', 
    '54': '            </TransformationChoice>', 
    '55': '          </ConditionalTransformation>', 
    '56': '        </TransformationSequence>', 
    '57': '      </RepeatedMove>', 
    '58': '      </TransformationSequence>', 
    '59': '    </TransformationSequence>', 
    '60': '  </Algorithm>'}

    #UPDATE step=36, repeat=10000, t_0=300, t_i=300
    Algorithm.update({'5': '    <TransformationSequence repeats="{}" weight="1.0">'.format(step)})
    Algorithm.update({'9': '        <repeats>{}</repeats>'.format(repeat)})
    Algorithm.update({'10': '        <tend>{}.0</tend>'.format(t_i)})
    Algorithm.update({'12': '        <tstart>{}.0</tstart>'.format(t_0)})
    # '50': '              <kB>0.0019858775</kB>'

    #dh_id = 2, sigma = 3
    AbsoluteDHLib = {
    '1': '      <SetDihedralAbsolute weight="1.0">', 
    '2': '        <dihedral_id>0</dihedral_id>', 
    '3': '        <value>3.14159</value>', 
    '4': '      </SetDihedralAbsolute>'}

    #dh_id = 2, sigma = 3
    RelativeDHLib = {
    '1': '      <SetDihedralRelative weight="1.0">', 
    '2': '        <dihedral_id>0</dihedral_id>', 
    '3': '        <value>0.174533</value>', 
    '4': '      </SetDihedralRelative>'}

    #dh_id = 2, sigma = 4, distribution_type = 3
    RelativeRandomDHLib = {
    '1': '              <SetDihedralRelativeRandom weight="1.0">', 
    '2': '                  <dihedral_id>1</dihedral_id>', 
    '3': '                  <distribution type="gaussian">', 
    '4': '                  <sigma>0.523599</sigma>', 
    '5': '                  <mean>0.0</mean>', 
    '6': '                  </distribution>', 
    '7': '              </SetDihedralRelativeRandom>'}

    for dh_id in range(0, len(DH_atoms)):
        os.mkdir('{}_SIMONAScanStep'.format(dh_id))
        os.chdir('{}_SIMONAScanStep'.format(dh_id))

        with open('{}_dihedralscan.xml'.format(dh_id), 'w') as outfile:
            XMLLines = []
            AlgorithmLimits = []
            with open('../SIMONA_inputs/tmp.xml','r') as infile: #'../SIMONA_inputs/tmp.xml'
                for line in infile:
                    XMLLines.append(line)

                for index, line in enumerate(XMLLines):
                    if "<Algorithm>" in line:
                        AlgorithmLimits.append(index)
                    if "</Algorithm>" in line:
                        AlgorithmLimits.append(index)

                #write all lines from 0 to <Algorithm>
                for idx in range(0,AlgorithmLimits[0]): # not including Algorithm
                    outfile.write(XMLLines[idx])
                #print(Algorithm.get(str(6)))
                #ALGORITHM section
                for AlgorithmLine in range(1, len(list(Algorithm))+1):

                    if Algorithm.get(str(AlgorithmLine)) == "AbsoluteTransformation":
                        AbsoluteDHLib.update({'2': '        <dihedral_id>{}</dihedral_id>'.format(dh_id)})
                        AbsoluteDHLib.update({'3': '        <value>{}</value>'.format(sigma_0)})
                        for idx in range(1, len(list(AbsoluteDHLib))+1):
                            outfile.write(AbsoluteDHLib.get(str(idx)) + "\n")
                    elif Algorithm.get(str(AlgorithmLine)) == "RelativeTransformation":
                        RelativeDHLib.update({'2': '        <dihedral_id>{}</dihedral_id>'.format(dh_id)})
                        RelativeDHLib.update({'3': '        <value>{}</value>'.format(angle)})
                        for idx in range(1, len(list(RelativeDHLib))+1):
                            outfile.write(RelativeDHLib.get(str(idx)) + "\n")
                    elif Algorithm.get(str(AlgorithmLine)) == "RelativeRandomTransformation":
                        RelativeRandomDHLib.update({'4': '                  <sigma>{}</sigma>'.format(MC_sigma)})
                        for dh_id2 in range(0, len(DH_atoms)):
                            if dh_id2 == dh_id:
                                continue
                            else:
                                RelativeRandomDHLib.update({'2': '                  <dihedral_id>{}</dihedral_id>'.format(dh_id2)})
                                for idx in range(1, len(list(RelativeRandomDHLib))+1):
                                    outfile.write(RelativeRandomDHLib.get(str(idx)) + "\n")
                    else:
                        outfile.write(Algorithm.get(str(AlgorithmLine))+ "\n")


                #Last lines of XML file
                #(AlgorithmLimits)
                for idx in range(AlgorithmLimits[1]+1,len(XMLLines)): # not including Algorithm
                    outfile.write(XMLLines[idx])

        os.chdir('../')



def SIMONAPoempp(SIMONAPATHS): #TODO: working
    """
    Process 
    -------------------------
    It runs the SIMONA calculations with Poempp

    INPUT
    -------------------------
    1. XML inputs files geneate with XMLScanModification

    OUTPUT
    -------------------------
    1. Simona outputs. 
    """

    #get the name of all the folders with name id_SIMONAScanStep
    SimonaFolders = sorted(glob.glob('*SIMONAScanStep')) #list with all the directories with name *SIMONAScanStep (as strings)
    #print(SimonaFolders)
    for idx in range(0, len(SimonaFolders)):
        #print("I am in working in folder {}".format(idx))
        os.chdir("{}_SIMONAScanStep".format(idx))
        
        POEMPP = '{}poempp'.format(SIMONAPATHS.get('PoemppDir'))
        poemppProcess = subprocess.Popen([POEMPP, "{}_dihedralscan.xml".format(idx)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding="utf-8")
        out, err = poemppProcess.communicate()
        poempp_out = []
        poempp_out.append(out)
        with open('run.out', 'w') as outfile:
            for line in poempp_out:
                outfile.write(line)
                if "SIMONA run finished." in poempp_out:
                    #print("Simulation in {} is finished".format(SimonaFolders[idx]))
                    poemppProcess.terminate()
            
        os.chdir('../')


def PostProcessSIMONA1(DHList): #TODO: Post process for XMLModi1
    """
    Process 
    -------------------------
    Reads all the output files in SIMONAScanStep folders and create a Results directory summarizing all the data obtained and plots. 

    INPUT
    -------------------------
    1. Outputs from poempp SIMONA calculation storage in  SIMONAScanStep folders

    OUTPUT
    -------------------------
    1. Generate the folder Results
    2. Get the raw data for each SIMONAScanStep folder
    3. Generate summary file of energies
    4. Create a plot per DH calculation (subplot).
    5. Create a summary plot with all DH calculations to determine the DH with bigger contribution to energy.
    
    """
    SimonaFolders = sorted(glob.glob('*SIMONAScanStep')) #list with all the directories with name *SIMONAScanStep (as strings)
    #print(sorted(SimonaFolders))
    os.mkdir('Results')
    os.chdir('Results')
    ScoreList = []
    Data = []
    for folderID in  range(0, len(SimonaFolders)):
        #print(folder)
        DataInFolder = []
        with open('../{}/trajectory.pdb'.format(SimonaFolders[folderID]), 'r') as infile:
            for line in infile:
                if '<E>' in line:
                    #print(line)
                    line = line.split()
                    DataInFolder.append(np.format_float_scientific(float(line[-1]), precision=4)) # last value is the total energy of the calculation 
            #Get min and max
            tmpDATA = [float(x) for x in DataInFolder]
            tmpArray = np.array(tmpDATA)
            tmpMin = np.format_float_scientific(np.min(tmpArray), precision=3)
            tmpMax = np.format_float_scientific(np.max(tmpArray), precision=3)
            tmpDiff = np.format_float_scientific(np.max(tmpArray) - np.min(tmpArray), precision=3)
            #print("{} --> min: {}, max : {}, diff : {}".format(folder, tmpMin, tmpMax, tmpDiff))
            ScoreList.append([SimonaFolders[folderID], tmpMin, tmpMax, tmpDiff])

        Data.append(DataInFolder)

    #After getting all the data, get the DataFrame file
    MyData = pd.DataFrame(Data) #pandas convert the Data nested list into a DataFrame
    MyData = MyData.transpose() #pandas do the transpose
    #print(MyData)
    MyData.to_csv(r'AllEnergies_outfile.txt', header=None, index=None, sep=' ', mode='a') #save the values as columns per simulation.
    def SortedScoreList(sublist):            
        sublist.sort(key = lambda x: x[1])            
        return sublist 


    print("\n Ranking of Dihedrals: \n")
    with open('FinalScore.dat', 'w') as outfile:
        for ranking, candidate in enumerate(SortedScoreList(ScoreList)[::-1]):
            
            print("{}. {} --> min: {}, max : {}, delta E : {} kcal/mol".format(ranking+1, candidate[0], candidate[1], candidate[2], candidate[3]))
            outfile.write("{}. {} --> min: {}, max : {}, delta E : {} kcal/mol \n ".format(ranking+1, candidate[0], candidate[1], candidate[2], candidate[3]))
    os.chdir('../')

def PlotDHProfiles(FileName, AngleList, EnergyList, TitlePlot):
    fig = plt.figure()
    x = np.array(AngleList)
    y = np.array(EnergyList)
    plt.scatter(x, y, color='r',zorder=1)
    plt.plot(x, y,color='b',zorder=2)
    plt.title(TitlePlot)
    plt.xlim([-190, 190]) #[np.min(x), np.max(x)]
    plt.ylim([np.min(y)-1, np.max(y)+1])
    plt.xlabel(r'Angle (Degree)', fontsize=14)
    plt.ylabel(r'Energy (kcal/mol)', fontsize=14)
    fig.tight_layout()
    plt.savefig('{}.png'.format(FileName), dpi=300)
    plt.clf()

def PostProcessSIMONA2(DHList, DHindex):
    """
    Process 
    -------------------------
    Reads all the output files in SIMONAScanStep folders and create a Results directory summarizing all the data obtained,  plots all profiles and extract the coordinates to feed Turbomole WANO. 

    INPUT
    -------------------------
    1. Outputs from poempp SIMONA calculation storage in  SIMONAScanStep folders
    2. DH_atoms list and DH_index list from DihedralInfo()

    OUTPUT
    -------------------------
    1. Generate the folder Results
    2. Get the raw data for each SIMONAScanStep folder : AllEnergies_outfile.txt
    3. Create a summary : FinalScore.dat with the ranking of dihedrals
    4. Create a plot per DH calculation (Profiles dir), depends on: PlotDHProfiles(FileName, AngleList, EnergyList, TitlePlot).
    5. Extract all the coordinates in PDB and xyz format to use in DFT-Turbomole WANO. Depends on: ExtractConf(SelectedFrames, SimonaFolder)
    
    """

    SimonaFolders = sorted(glob.glob('*SIMONAScanStep')) #list with all the directories with name*SIMONAScanStep (as strings)
    #print(sorted(SimonaFolders))
    os.mkdir('Results')
    os.chdir('Results')
    os.mkdir('Profiles')
    ScoreList = []
    Data = []
    for folderID in  range(0, len(SimonaFolders)):
        #Energy DATA
        DataInFolder = [] #only Energy values
        with open('../{}/trajectory.pdb'.format(SimonaFolders[folderID]), 'r') as infile:
            for line in infile:
                if '<E>' in line:
                    #print(line)
                    line = line.split()
                    DataInFolder.append(np.format_float_scientific(float(line[-1]), precision=4)) #last value is the total energy of the calculation
        #Angle DATA
        step_angles= []
        angles = [] #only angle values
        temp_traj = '../{}/trajectory.pdb'.format(SimonaFolders[folderID])
        cmd.load(temp_traj, "traj")
        cmd.set("retain_order",1)
        num_frames = cmd.count_frames ("traj")
        #print(SimonaFolders[folderID], num_frames)

        def special_round(n, level):
            round0 = round(n, level)
            if round0 == -0.0:
                return abs(round0)
            else:
               return round0 

        for idx in range(1, num_frames):
            angle = cmd.get_dihedral('name ' + DHList[folderID][0],'name '+ DHList[folderID][1] ,'name ' + DHList[folderID][2],'name ' + DHList[folderID][3], idx)
            #step_angles.append([idx, round(angle,0)])
            step_angles.append([idx, special_round(angle, -1)])
            angles.append(special_round(angle, -1))
        cmd.delete("traj")

        #write a file with angle vs energy in simulation directory
        with open('angle_energy_{}.csv'.format(SimonaFolders[folderID]), 'w') as f:
            writer = csv.writer(f)
            writer.writerows(zip(angles, DataInFolder))

        #Angle vs Energy: The lowest energy values per angle step must be choosen...
        FinalUsefulData = []# angle vs energy to plot

        AnglesLib = list(sorted(set(angles))) #types of angles, this sort the angles by group.

        step_angle_energy = []
        for ai in range(0, len(angles)):
            step_angle_energy.append([ai, angles[ai], DataInFolder[ai]]) 
        for AngleType in AnglesLib:
            tmp_list = []
            for candidate in step_angle_energy:
                if candidate[1]  == AngleType:
                    tmp_list.append(candidate)

            #find the lowest energy between AngleType
            def SORT_inAngleType(sublist):
                sublist.sort(key = lambda x: x[2])
                return sublist

            ChoosenCandidate = SORT_inAngleType(tmp_list)[-1]
            #print('In {} was chose: {}'. format(SimonaFolders[folderID], ChoosenCandidate))
            FinalUsefulData.append(ChoosenCandidate) 
        #print(FinalUsefulData)

        #CREATE PLOT
        x = []
        y = []
        SelectedFrames = []
        for idx in range(0, len(FinalUsefulData)):
            SelectedFrames.append(int(FinalUsefulData[idx][0]))
            x.append(FinalUsefulData[idx][1])
            y.append(float(FinalUsefulData[idx][2]))
        PlotDHProfiles(SimonaFolders[folderID], x, y, '{} in {}'.format(DHList[folderID], SimonaFolders[folderID]))
        #Extract PDB and xyz coordinates for next steps
        ExtractConf(SelectedFrames, SimonaFolders[folderID])
        os.system('mv angle_energy_{}.csv -t {}'.format(SimonaFolders[folderID], SimonaFolders[folderID]))

        #Final Score data
        def SORT_EnergyInSimulation(sublist):
            sublist.sort(key = lambda x: x[2])
            return sublist
        SortedUsefulData = SORT_EnergyInSimulation(FinalUsefulData)
        tmpMin = float(SortedUsefulData[-1][2])
        #print(type(SortedUsefulData[-1][2]))
        tmpMax = float(SortedUsefulData[1][2])
        tmpDeltaE = abs(tmpMax - tmpMin)
        ScoreList.append([SimonaFolders[folderID], tmpMin, tmpMax, tmpDeltaE])
        Data.append(DataInFolder)


    #After getting all the data, get the DataFrame file
    MyData = pd.DataFrame(Data) #pandas convert the Data nested list into a DataFrame
    MyData = MyData.transpose() #pandas do the transpose
    #print(MyData)
    MyData.to_csv(r'AllEnergies_outfile.txt', header=None, index=None, sep=' ',mode='a') #save the all the values (per step) as columns per simulation.
    
    def SortedScoreList(sublist):            
        sublist.sort(key = lambda x: x[3])            
        return sublist 
    print("\n Ranking of Dihedrals: \n")
    with open('FinalScore.dat', 'w') as outfile:
        for ranking, candidate in enumerate(SortedScoreList(ScoreList)[::-1]):

            print("{}. {} --> min: {}, max : {}, delta E : {} kcal/mol".format(ranking+1, candidate[0], round(candidate[1],3), round(candidate[2],3), round(candidate   [3], 3)))
            outfile.write("{}. {} --> min: {}, max : {}, delta E : {} kcal/mol \n".format(ranking+1, candidate[0], round(candidate[1],3), round(candidate[2],3), round  (candidate[3], 3)))
    os.system('mv *.png -t Profiles')
    os.chdir('../')

    #Get the right tar in the main DIR and create the DH_index file
    BestCandidate = SortedScoreList(ScoreList)[::-1][0][0]
    print('Extracting {} structures and dihedral info'.format(BestCandidate))
    os.system('cp Results/{}/structures.tar.gz .'.format(BestCandidate))
    print(DHindex.get(BestCandidate))

    tmpDictionary = DHindex.get(BestCandidate)
    id_list = []
    for ix in tmpDictionary:
        #print(ix)
        id_list.append(int(tmpDictionary.get(ix)))
    DHdic = {'id_list': id_list}
    with open('Result_DHInfo.yml', 'w') as outfile:
        yaml.dump(DHdic, outfile, default_flow_style=False)

    #with open('Result_DHInfo.txt', 'w') as outfile:
    #    outfile.write(str(DHindex.get(BestCandidate)))



def ExtractConf(SelectedFrames, SimonaFolder): #Selected Frames are the lowest energy structures from Analysis in PostProcessSIMONA2
    os.mkdir(SimonaFolder)
    os.chdir(SimonaFolder)
    os.mkdir('Coords_PDB')
    os.mkdir('Molecules')
    TrajectoryFile = "../../{}/trajectory.pdb".format(SimonaFolder)
    cmd.load(TrajectoryFile, "traj")
    cmd.set("retain_order", 1)
    #frames = cmd.count_frames ("traj")
    for idx, frame in enumerate(SelectedFrames):
        coord_filename = '{}step_frame{}.pdb'.format(idx, frame)
        cmd.save(coord_filename, "traj", frame)
        subprocess.run(["obabel", "-i", "pdb", coord_filename, "-o", "xyz", "-O",'{}step.xyz'.format(idx)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    #MAKE THE TAR
    #AllGoodFiles = glob.glob('*.xyz')
    #with tarfile.open('structures.tar.gz','w:gz') as tar:
    #    for file in AllGoodFiles:
    #        tar.add(file)

    os.system('mv *step*.pdb -t Coords_PDB')
    os.system('mv *step*.xyz -t Molecules')
    with tarfile.open('structures.tar.gz','w:gz') as tar:
        tar.add('Molecules')
    cmd.delete("traj")
    os.chdir('../')



def DFT_TurbomolePreprocessor(MOL, Step, DH_index):
    """
    Process 
    -------------------------
    Generates conformations after Absolute rotation of dihedral provided by user. Here there is no pre-optimization step (pymol).

    INPUT
    -------------------------
    1. MOL name
    2. number of rotation Step
    3. DH_index list provided by user

    OUTPUT
    -------------------------
    1. PDB structures after rotation
    2. XYZ coordinates inputs for DFT-Turbomole WANO
    3. XYZ tar file
    4. DFT-TURBOMOLE DH_index file (just in case)
    """
    #Pymol, angle must be in degrees, no radians
    CoordFile = '{}.pdb'.format(MOL)


    #Identify atom names for DH_index
    AtomsLib = {}
    DH_indexLIB = {}
    with open(CoordFile, 'r') as inputfile:
        for line in inputfile:
            if 'ATOM' in line:
                line = line.split()
                index = line[1]
                AtomName = line[2]
                AtomsLib[AtomName] = index
                for idx in DH_index:
                    if index == idx :
                        DH_indexLIB[AtomName] = index
            elif 'HETATM' in line:
                line = line.split()
                index = line[1]
                AtomName = line[2]
                AtomsLib[AtomName] = index
                for idx in DH_index:
                    if index == idx :
                        DH_indexLIB[AtomName] = index
            else:
                print('PDB format is not right, please check your PDB file input')

#    with open('Result_DHInfo.txt', 'w') as outfile:
#        outfile.write(str(DH_indexLIB))

    tmpDictionary = DH_indexLIB
    id_list = []
    for ix in tmpDictionary:
        #print(ix)
        id_list.append(int(tmpDictionary.get(ix)))
    DHdic = {'id_list': id_list}
    with open('Result_DHInfo.yml', 'w') as outfile:
        yaml.dump(DHdic, outfile, default_flow_style=False)


    cmd.load(CoordFile, "MOL")
    cmd.set("retain_order", 1)

    MyAngle = 360/Step 

    for step in range(1, Step+1):

        cmd.set_dihedral('index {}'.format(DH_index[0]), 'index {}'.format(DH_index[1]), 'index {}'.format(DH_index[2]), 'index {}'.format(DH_index[3]), step * MyAngle)
        cmd.save('{}_step.pdb'.format(step), "MOL")
        subprocess.run(["obabel", "-i", "pdb", '{}_step.pdb'.format(step), "-o", "xyz", "-O",'{}step.xyz'.format(idx)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    #MAKE THE TAR
    #AllGoodFiles = glob.glob('*.xyz')
    #with tarfile.open('structures.tar.gz','w:gz') as tar:
    #    for file in AllGoodFiles:
    #        tar.add(file)

    os.mkdir('Coords_PDB')
    os.mkdir('Coords_xyz')
    os.system('mv *step*.pdb -t Coords_PDB')
    os.system('mv *step*.xyz -t Molecules')

    with tarfile.open('structures.tar.gz','w:gz') as tar:
        tar.add('Molecules')


#-----------------__MAIN__-----------------#

if __name__ == '__main__': #TODO: make the workflow control work

    #CLUSTER SELECTION for SIMONAPATHS
    SIMONAPATHS = {}
    CLUSTER = os.uname()[1]

    print("USING CUSTOM SIMONA PATH")
    with open('settings.sh','r') as infile:
        for line in infile:
            if 'CUSTOMSIMONAPYTHONPATH' in line:
                line = line.split()
                SIMONAPATHS['PythonDir'] = '{}/'.format(line[-1][1:-1])
            if 'CUSTOMSIMONAPOEMPPPATH' in line:
                line = line.split()
                SIMONAPATHS['PoemppDir'] = '{}/'.format(line[-1][1:-1])

    print('your are working in: {} with this simona path: {}'.format(CLUSTER,SIMONAPATHS.get('PoemppDir')[0:-5]))
    
    #get the inputs from SimStack
    settings = {}
   
    with open('rendered_wano.yml', 'r') as file:
        wano_file = yaml.full_load(file)
        #1.Scan settings
        settings['MoleculeName'] = wano_file['Molecule Name']
        settings['MoleculeNetCharge'] = wano_file['Molecule net charge']
        settings['RotationSteps'] = wano_file['Rotation steps']

        #2.Dihedral Score (SIMONA)
        settings['DihedralScore'] = wano_file['Dihedral Score to indetify dihedrals (SIMONA)']  #TRUE or FALSE

        #2.1 Smile input
        settings['SmileON_OFF'] = wano_file['Input Dihedral-Score']['SMILE code'] #TRUE or FALSE
        settings['SmileCode'] = wano_file['Input Dihedral-Score']['Input SMILE Code'] #code

        #2.2 PDB input
        settings['PDBInputON_OFF'] = wano_file['Input Dihedral-Score']['Coordinates'] #TRUE or FALSE
        settings['SimonaCoord_File'] = wano_file['Input Dihedral-Score']['Structure'] #PDB file

        #3. Turbomol Scan calculation  (later)
        settings['TurbomoleON_OFF'] = wano_file['Generate structure inputs for DFT-Turbomole'] #TRUE or FALSE
        settings['TurboCoord_File'] = wano_file['Coordinates and dihedral IDs inputs']['Structure']  #PDB file
        
        id_list = []
        AllatomIDS = wano_file['Coordinates and dihedral IDs inputs']['Torsion atom ids']
        for ids in AllatomIDS:
            ai = ids['Atom id']
            id_list.append(ai)
        settings['DH_index'] = id_list
        
    #general values
    MoleculeName = settings.get('MoleculeName')
    Charge = settings.get('MoleculeNetCharge')
    Step = settings.get('RotationSteps')
    Angle = (360 / int(Step)) * 0.0174533 #radians
    print('Working with : {}, Net charge : {}, Rotation steps : {}, Delta angle: {} radians '.format(MoleculeName, Charge, Step, Angle))

    #Check the name and convert it to capitals if it is needed
    if MoleculeName.isupper() == False:
        MoleculeName = MoleculeName.upper()
        #print("Molecule id name: ", MoleculeName[:3])

    else:
        print("Molecule id name: ", MoleculeName[:3])
        
    MOL = MoleculeName[:3]

    if settings.get('DihedralScore') == True: #run simona scan
        
        if settings.get('SmileON_OFF') == True: #SMILE code input
            #Defensive statement when SmileCode is not provided for some reason:
            assert not settings.get('SmileCode') == 'NONE', 'SMILE is NONE, please try again'
            #functions , DihedralInfo(AtomTypes) = DH_atoms
            SMILE_Preprocessor(MOL, settings.get('SmileCode'))
            AcPyPeProcess(MOL, Charge)
            AtomTypes = SIMONA_files(MOL,SIMONAPATHS) 
            #XMLScanModification(DihedralInfo(AtomTypes), Angle, Step)
            DH_atoms, DH_index = DihedralInfo(AtomTypes)
            XMLScanModification2(DH_atoms, angle=0.174533, step=36, sigma_0=3.14159, MC_sigma=1.523599, repeat=10000, t_0=300, t_i=300)
            SIMONAPoempp(SIMONAPATHS)
            PostProcessSIMONA2(DH_atoms, DH_index)
        
        elif settings.get('PDBInputON_OFF') == True: #coordinate file

            #change PDB name to MOL.pdb
            os.system("cp {} {}.pdb".format(settings.get('SimonaCoord_File'), MOL))
            AcPyPeProcess(MOL, Charge)
            AtomTypes = SIMONA_files(MOL,SIMONAPATHS) 
            DH_atoms, DH_index = DihedralInfo(AtomTypes)
            XMLScanModification2(DH_atoms, angle=0.174533, step=36, sigma_0=3.14159, MC_sigma=1.523599, repeat=10000, t_0=300, t_i=300)
            SIMONAPoempp(SIMONAPATHS)
            PostProcessSIMONA2(DH_atoms, DH_index)
        else:

            with open('Error_report', 'w') as output:
                output.write("Input is wrong. Check error file")
                    
    elif settings.get('TurbomoleON_OFF') == True: #run Turbomole scan
        # HOW to conect to DFT-TURBOMOLE WANO?
        print('This part is on going')
        os.system("cp {} {}.pdb".format(settings.get('SimonaCoord_File'), MOL))
        DFT_TurbomolePreprocessor(MOL, Step, settings.get('DH_index'))

    else:
        with open('Error_report', 'w') as output:
            output.write("Input is wrong. Check error file")
