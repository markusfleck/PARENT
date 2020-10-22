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







#ifndef BAT_TOPOLOGY_H
#define BAT_TOPOLOGY_H


#include <vector>



class BAT_Topology {

public:
    BAT_Topology(std::vector< std::vector< std::vector <int> > > *bonds_table_in, std::vector <int> *backboneIN, std::vector< std::vector <int> > *roots_in);
    std::vector< std::vector <int> > dihedrals; //contains a list of all dihedrals found

private:
    std::vector< std::vector< std::vector <int> > > bonds_table; //contains a list of all bonds in the system with proper offset, with bonds_table[#atom-1][#bond-1][0] the bonded atom and bonds_table[atom-1][#bond-1][1] the type of the bond (0==physical, 1==pseudo)
    std::vector< std::vector <int> > roots; //contains the root atoms for every molecule in the system
    std::vector <int> backbone; //contains a list of the backbone atoms in the system
    bool *isBackbone;
    int *map_BAT_to_real;
    int *map_real_to_BAT;
    int current_BAT;//counter for the next BAT number to be assigned
    
    void do_mapping_depth_first();
    void add_node_recursive(int real_node);
    void do_mapping_breadth_first();
    void add_node_breadth(int real_node, std::vector <int> *to_be_processed_in);
    int find_closest_bond(int real_node, std::vector<int> exclude_list);
    void create_BAT();
    void add_phases();

};

#endif

