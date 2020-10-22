//    Class for BAT_builder, a program to convert a molecular dynamics trajectory from Cartesian to internal bond-angle-torsion coordinates
//    Copyright (C) 2015  Markus Fleck (member of the laboratory of Bojan Zagrovic, University of Vienna)
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License  version 3
//    as published by the Free Software Foundation.

//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.

//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.





//    A scientific publication about this program has been released in the Journal of Chemical Theory and Computation:
//		"PARENT: A Parallel Software Suite for the Calculation of Configurational Entropy in Biomolecular Systems"
//		DOI: 10.1021/acs.jctc.5b01217

//   We kindly ask you to include this citation in works that publish
//   results generated using this program or any modifications of it.







#ifndef TOPOLOGY_H
#define TOPOLOGY_H


#include <vector>
#include <string>


class Topology {

public:
    Topology(std::string topfileIN,std::vector <std::string> backboneAtomNamesIN);
    ~Topology();

    std::vector< std::vector< std::vector <int> > > bonds_table; //contains a list of all bonds in the system with proper offset, with bonds_table[#atom-1][#bond-1][0] the bonded atom and bonds_table[atom-1][#bond-1][1] the type of the bond (0==physical, 1==pseudo)
    std::vector <int> backbone; //backbone atoms in the system       
    std::vector< std::vector <int> > roots; //root atoms for all molecules in the system. root atoms need to be connected root[0] to root[1] and root[1] to root[2]. root[0] needs to be terminal, root[1] needs to be connected only to root[1] and terminal atoms. root[2] needs to be connected to non-terminal atoms.
    std::vector<float> masses; //diagonal of the mass matrix of the system    
    std::vector <std::string> residues; //contains the name of the residue for every atom in the system
    std::vector <int> residueNumbers; //contains the number of the residue for every atom in the system
    std::vector <std::string> atomNames; //contains the name of every atom in the system
    std::vector <std::string> belongsToMolecule; //the name of the molecule each atom belongs to

private:
    std::vector< std::vector<int> > moleculetypes;  //contains a list of the atoms in each moleculetype
    std::vector< std::vector<int> > molecules; //contains a list of atoms for each molecule in the system with correct offset, as defined in the [ molecules ] section
    std::vector< std::vector<int> > angles; //contains a list of angles for each molecule in the system with correct offset, as defined in the [ molecules ] section
    std::vector< std::vector< std::vector <int> > > moleculetypebonds; //contains a list  of all bonds for every moleculetype
    std::vector< std::vector< std::vector <int> > > bonds;  //contains a list of atoms for each molecule in the system with correct offset, as defined in the [ molecules ] section
    std::vector<std::string> moleculetypenames;  //names of the moleculetypes
    std::vector<std::string> moleculenames; //contains a full list of molecule names of the system as defined in the [ molecules ] section (water 1000 => water, water, water.......)
    std::vector< std::vector <std::string> > moleculetypeResidues; //contains the name of the residue for every moleculetype and every atom
    std::vector< std::vector <int> > moleculetypeResidueNumbers; //contains the name of the residue number for every moleculetype and every atom
    std::vector< std::vector <std::string> > moleculetypeAtomNames; //contains the name of the atom for every moleculetype and every atom
    std::vector< std::vector<float> >  moleculetypemasses; //masses for each moleculetype
    std::vector< std::vector<int> >  moleculetypebackbone; //backboneatoms for each moleculetype
    std::vector <std::string> backboneAtomNames; //names of the atomtypes of the backbone for building phase angles from dihedrals, which is done in BAT_trajectory.cpp
    long int current;
    long int total_atoms;
    unsigned short int newclustflag,newclustmassflag,newclustatomflag,newclustbackboneflag;
    std::string top_dir, top_file;

    unsigned short int add_bond(int atom1, int atom2);
    unsigned short int add_mass(float mass);
    unsigned short int add_atom(std::string resname, int resnumber, std::string atomname);
    unsigned short int add_backbone(int atom);
    void new_moleculetype(std::string mol);
    int sort();
    std::string get_dir_file(std::string file);
    unsigned short int check_comment(std::string line);
    int preprocess(std::string infile,std::stringstream* myVirtualOFile) ;
    int read_topology();
    void make_bonds_table();
    void get_roots();
    void get_roots_override(double perc);
};

#endif
