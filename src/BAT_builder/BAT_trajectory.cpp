//    BAT_builder, a program to convert a molecular dynamics trajectory from Cartesian to internal bond-angle-torsion coordinates
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







#include "../util/io/io.h"



using namespace std;

class BAT_Trajectory {

public:
    //constructor for the .xtc to .bat conversion
    BAT_Trajectory(string trjFileI,string trjFileO, vector< vector <int> > *dihedralsIn, vector <float> *massvec_in, vector <string> *residues_in,vector <int> *residueNumbers_in,vector <string> *atomNames_in,vector <string> *belongsToMolecule_in, int double_prec_in=0) {

        //definition, allcation etc
        pi=acos(-1); //define pi
        dihedrals_top=(*dihedralsIn); //a list of all dihedrals in the system (index0=atom1, index1=atom2, index2=atom3, index3=atom4, index4=type(physical=0,pseudo=1,improper=-1), index5="parent"dihedral(for phaseangle))
        massvec=(*massvec_in); //contains the masses all atoms of the system
        residues=(*residues_in); //contains the name of the residue for every atom in the system
        residueNumbers=(*residueNumbers_in); //contains the number of the residue for every atom in the system
        atomNames=(*atomNames_in); //contains the name of every atom in the system
        belongsToMolecule=(*belongsToMolecule_in); //the name of the molecule each atom belongs to
        x=(rvec*) malloc((dihedrals_top.size()+3)*sizeof(rvec)); //allocate memory for N Cartesian coordinate vectors
        if(x==NULL) {
            cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
            exit(EXIT_FAILURE);
        }
        trjFileIn=trjFileI; //name of the input file
        trjFileOut=trjFileO; //name of the output file
        double_prec=double_prec_in; //should the .bat trajectory be written in double precision?
        numframes=0;

        convert_xtc_to_BAT(); //start the procedure
    }
    //constructor for the .bat to .xtc conversion
    BAT_Trajectory(string trjFileI,string trjFileO) {
        pi=acos(-1); //define pi
        trjFileIn=trjFileI; //name of the input file
        trjFileOut=trjFileO; //name of the output file

        convert_BAT_to_xtc(); //start the procedure
    }

    int double_prec;
    vector< vector <int> > dihedrals_top;
    vector <float> massvec;
    vector <string> residues;
    vector <int> residueNumbers;
    vector <string> atomNames;
    vector <string> belongsToMolecule;
    string trjFileIn,trjFileOut;
    rvec *x;


    double root_origin_cartesian[3];
    double root_origin_phi;
    double root_origin_theta;
    double root_origin_dihedral;

    double rootbond[2];
    double rootangle;

    double *bonds;
    double *angles;
    double *dihedrals;
    double *bondsFull,*anglesFull;
    float *masses;

    int step;
    float time,prec;
    matrix box;

    double pi;
    int numframes;




    int convert_xtc_to_BAT() {
        XDRFILE *xdi;
        int natoms;
        read_xtc_natoms((char*)trjFileIn.c_str(),&natoms);//check the number of atoms in the .xtc file (function from XTC Library (GROMACS) ))
        if((unsigned int)(natoms)!=dihedrals_top.size()+3) {
            cerr << "ERROR: NUMBER OF ATOMS IN THE TRAJECTORY("<<natoms<<") DOES NOT FIT THE TOPOLOGY("<<dihedrals_top.size()+3<<")!\n";
            exit(EXIT_FAILURE);
        }
        bondsFull=(double*)malloc((dihedrals_top.size()+2)*sizeof(double)); //allocate memory for bonds
        anglesFull=(double*)malloc((dihedrals_top.size()+1)*sizeof(double)); //allocate memory for angles
        if((bondsFull==NULL)||(anglesFull==NULL)) {
            cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
            exit(EXIT_FAILURE);
        }
        bonds=&(bondsFull[2]);
        angles=&(anglesFull[1]);
        dihedrals=(double*)malloc(dihedrals_top.size()*sizeof(double)); //allocate memory for N-3 dihedrals
        if(dihedrals==NULL) {
            cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
            exit(EXIT_FAILURE);
        }
        xdi = xdrfile_open(trjFileIn.c_str(),"r");//open the .xtc input file
        if(xdi!=NULL) {
            ofstream outfile(trjFileOut.c_str(), ios::binary | ios::out);//open the .bat binary output file (function from XTC Library (GROMACS) )
            if(outfile.is_open()) {
                if(write_BAT_header(&outfile,double_prec,numframes,&dihedrals_top, &massvec, &residues,&residueNumbers,&atomNames,&belongsToMolecule)!=0) {//write the header information to the .bat file
                    cerr << "ERROR WRITING THE HEADER TO "<<trjFileOut.c_str()<<" !\n";
                    exit(EXIT_FAILURE);
                }
                //for all frames in the trajectory
                while(read_xtc(xdi,natoms,&step,&time,box,x,&prec)==exdrOK) { //read a frame from the trajectory (function from XTC Library (GROMACS) )
                    for(unsigned int j=0; j<dihedrals_top.size()+3; j++) {
                        x[j][0]*=10.0; //for comparison with previous work, use Angstroem
                        x[j][1]*=10.0;
                        x[j][2]*=10.0;
                    }
                    convert_to_BAT(); //convert the current frame to BAT coordinates
                    convert_dihedrals_to_phaseangles(); //calculate the phase angles
                    shift_dihedrals_0_2pi_range(); //shift the range of the dihedrals from [-pi,pi] to [0,2pi]
                    bondsFull[0]=rootbond[0];
                    bondsFull[1]=rootbond[1];
                    anglesFull[0]=rootangle;
                    if(write_BAT_frame(&outfile,double_prec, dihedrals_top.size(), time, prec, (float**)box, root_origin_cartesian,root_origin_theta,root_origin_phi,root_origin_dihedral,bondsFull,anglesFull,dihedrals)!=0) { //attach the current frame to the .bat trajectory
                        cerr << "ERROR WRITING FRAME "<<numframes+1<<" TO "<<trjFileOut.c_str()<<" !\n";
                        exit(EXIT_FAILURE);
                    }
                    numframes++; //count the number of frames
                }
                outfile.seekp(3*sizeof(int), ios::beg);//update "numframes" in the header
                outfile.write((char*)&numframes, sizeof(int));
                outfile.close();//close the output file

            }
            else {
                cerr << "ERROR: UNABLE TO OPEN FILE "<<trjFileOut.c_str()<<" !\n";
                exit(EXIT_FAILURE);
            }
            xdrfile_close(xdi);
        }
        else {
            cerr << "ERROR: UNABLE TO OPEN FILE "<<trjFileIn.c_str()<<" !\n";
            exit(EXIT_FAILURE);
        }
        cout<<"Processed "<<numframes<<" frames."<<endl;
        return 0;
    }


    int convert_BAT_to_xtc() {

        XDRFILE *xdo;
        ifstream infile(trjFileIn.c_str(), ios::binary | ios::in); //open the binary .bat input file
        if(infile.is_open()) {
            if (read_BAT_header(&infile,&double_prec,&numframes,&dihedrals_top,&massvec, &residues,&residueNumbers,&atomNames,&belongsToMolecule)) {
                exit(EXIT_FAILURE);   //read the header of the .bat file
            }
            xdo = xdrfile_open(trjFileOut.c_str(),"w"); //open the .xtc input file (function from XTC Library (GROMACS) )
            if(xdo==NULL) {
                cerr << "ERROR: UNABLE TO OPEN FILE "<<trjFileOut.c_str()<<" !\n";
                exit(EXIT_FAILURE);
            }
            int natoms=dihedrals_top.size()+3;
            bondsFull=(double*)malloc((dihedrals_top.size()+2)*sizeof(double)); //allocate memory bonds
            anglesFull=(double*)malloc((dihedrals_top.size()+1)*sizeof(double)); //allocate memory for angles
            dihedrals=(double*)malloc(dihedrals_top.size()*sizeof(double)); //allocate memory for dihedrals
            if((bondsFull==NULL)||(anglesFull==NULL)||(dihedrals==NULL)) {
                cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
                exit(EXIT_FAILURE);
            }
            bonds=&(bondsFull[2]);
            angles=&(anglesFull[1]);
            x=(rvec*) malloc((dihedrals_top.size()+3)*sizeof(rvec)); //allocate memory for N Cartesian coordinate vectors
            if(x==NULL) {
                cerr<<"ERROR ALLOCATNG MEMORY! ABORTING.\n";
                exit(EXIT_FAILURE);
            }
            step=0;
            while(read_BAT_frame(&infile,double_prec,dihedrals_top.size(), &time, &prec, (float **)box, root_origin_cartesian,&root_origin_theta,&root_origin_phi,&root_origin_dihedral,bondsFull,anglesFull,dihedrals)==0) { //for all frames in the .bat file
                rootbond[0]=bondsFull[0];
                rootbond[1]=bondsFull[1];
                rootangle=anglesFull[0];
                convert_phaseangles_to_dihedrals(); //calculate dihedrals from phaseangles
                convert_to_Cartesian(); //and convert to Cartesian coordinates
                for(unsigned int j=0; j<dihedrals_top.size()+3; j++) {
                    x[j][0]/=10.0; //for comparison with previous work, use Angstroem
                    x[j][1]/=10.0;
                    x[j][2]/=10.0;
                }
                write_xtc(xdo,natoms,step,time,box,x,prec); //write the current frame to the .xtc file
                step++;
            }
            infile.close();
            if(step!=numframes) {
                cerr << "ERROR: COULD ONLY PROCESS "<<step<<"FRAMES OUT OF "<<numframes<<" !\n";
                exit(EXIT_FAILURE);
            }
        } else {
            cerr << "ERROR: UNABLE TO OPEN FILE "<<trjFileIn.c_str()<<" !\n";
            exit(EXIT_FAILURE);
        }
        xdrfile_close(xdo);
        return 0;
    }

    void shift_dihedrals_0_2pi_range() {
        //shift the range of the dihedrals from [-pi,pi] to [0,2pi]
        for(unsigned int i=0; i<dihedrals_top.size(); i++) {
            dihedrals[i]=dihedrals[i]+2*pi-int((dihedrals[i]+2*pi)/(2*pi))*2*pi;
        }
    }




    void convert_to_BAT() {
        double b1vec[3];
        double b2vec[3];
        double b3vec[3];
        double ABxBC[3];
        double BCxCD[3];
        double prevBond;


        root_origin_cartesian[0]=x[dihedrals_top[0][0]-1][0];//define the first three external coordinates as the Cartesians of the first root atom
        root_origin_cartesian[1]=x[dihedrals_top[0][0]-1][1];
        root_origin_cartesian[2]=x[dihedrals_top[0][0]-1][2];

        b1vec[0]=root_origin_cartesian[0]; //define a vector from the origin to the first root atom
        b1vec[1]=root_origin_cartesian[1];
        b1vec[2]=root_origin_cartesian[2];

        b2vec[0]=x[dihedrals_top[0][1]-1][0]-x[dihedrals_top[0][0]-1][0];//define a vector from the first root atom to the second root atom
        b2vec[1]=x[dihedrals_top[0][1]-1][1]-x[dihedrals_top[0][0]-1][1];
        b2vec[2]=x[dihedrals_top[0][1]-1][2]-x[dihedrals_top[0][0]-1][2];

        b3vec[0]=x[dihedrals_top[0][2]-1][0]-x[dihedrals_top[0][1]-1][0];//define a vector from the second root atom to the third root atom
        b3vec[1]=x[dihedrals_top[0][2]-1][1]-x[dihedrals_top[0][1]-1][1];
        b3vec[2]=x[dihedrals_top[0][2]-1][2]-x[dihedrals_top[0][1]-1][2];


        root_origin_phi=atan2(b2vec[1],b2vec[0]); //define an additional external coordinate as the polar phi angle connecting root1 to root2

        rootbond[0]=sqrt(b2vec[0]*b2vec[0]+b2vec[1]*b2vec[1]+b2vec[2]*b2vec[2]); //calculate the bondlengths between the root atoms (internal coordinates)
        rootbond[1]=sqrt(b3vec[0]*b3vec[0]+b3vec[1]*b3vec[1]+b3vec[2]*b3vec[2]);

        root_origin_theta=acos(b2vec[2]/rootbond[0]);//define another external coordinate as the polar theta angle connecting root1 to root2


        ABxBC[0]=b1vec[1]*b2vec[2]-b1vec[2]*b2vec[1]; //define crossproducts  (b1vec x b2vec) and (b2vec x b3vec)
        ABxBC[1]=b1vec[2]*b2vec[0]-b1vec[0]*b2vec[2];
        ABxBC[2]=b1vec[0]*b2vec[1]-b1vec[1]*b2vec[0];
        BCxCD[0]=b2vec[1]*b3vec[2]-b2vec[2]*b3vec[1];
        BCxCD[1]=b2vec[2]*b3vec[0]-b2vec[0]*b3vec[2];
        BCxCD[2]=b2vec[0]*b3vec[1]-b2vec[1]*b3vec[0];


        //define a sixth external coordinate as the dihedral angle origin-root1-root2-root3
        root_origin_dihedral=atan2((ABxBC[2]*BCxCD[0] - ABxBC[0]*BCxCD[2])*b2vec[1]/rootbond[0] - (ABxBC[2]*BCxCD[1] - ABxBC[1]*BCxCD[2])*b2vec[0]/rootbond[0] - (ABxBC[1]*BCxCD[0] - ABxBC[0]*BCxCD[1])*b2vec[2]/rootbond[0], ABxBC[0]*BCxCD[0] + ABxBC[1]*BCxCD[1] + ABxBC[2]*BCxCD[2]);

        //calculate the angle between the root atoms (internal coordinate)
        rootangle=acos((b2vec[0]*b3vec[0]+b2vec[1]*b3vec[1]+b2vec[2]*b3vec[2])/rootbond[0]/rootbond[1]);



        //for all "later" dihedrals, perform somthing similar (only for internal coordinates)
        for(unsigned int i=0; i<dihedrals_top.size(); i++) {

            b1vec[0]=x[dihedrals_top[i][1]-1][0]-x[dihedrals_top[i][0]-1][0];//define a vector  from dihedral atom1 to atom2
            b1vec[1]=x[dihedrals_top[i][1]-1][1]-x[dihedrals_top[i][0]-1][1];
            b1vec[2]=x[dihedrals_top[i][1]-1][2]-x[dihedrals_top[i][0]-1][2];

            b2vec[0]=x[dihedrals_top[i][2]-1][0]-x[dihedrals_top[i][1]-1][0];//another from atom2 to atom3
            b2vec[1]=x[dihedrals_top[i][2]-1][1]-x[dihedrals_top[i][1]-1][1];
            b2vec[2]=x[dihedrals_top[i][2]-1][2]-x[dihedrals_top[i][1]-1][2];

            b3vec[0]=x[dihedrals_top[i][3]-1][0]-x[dihedrals_top[i][2]-1][0];//and a third from atom3 to atom4
            b3vec[1]=x[dihedrals_top[i][3]-1][1]-x[dihedrals_top[i][2]-1][1];
            b3vec[2]=x[dihedrals_top[i][3]-1][2]-x[dihedrals_top[i][2]-1][2];

            bonds[i]=sqrt(b3vec[0]*b3vec[0]+b3vec[1]*b3vec[1]+b3vec[2]*b3vec[2]);// the bondlength associated with dihedral[i] is the the distance between atom4 and atom3, so the length of b3vec
            prevBond=sqrt(b2vec[0]*b2vec[0]+b2vec[1]*b2vec[1]+b2vec[2]*b2vec[2]); //we also need the bondlength of b2vec for the calculation

            angles[i]=acos((b2vec[0]*b3vec[0]+b2vec[1]*b3vec[1]+b2vec[2]*b3vec[2])/bonds[i]/prevBond); //the angle associated with dihedral[i] is between atom4, atom3 and atom2, so between b2vec and b3vec.

            ABxBC[0]=b1vec[1]*b2vec[2]-b1vec[2]*b2vec[1]; //the formula for the dihedral angle contains the x-products between the vectors
            ABxBC[1]=b1vec[2]*b2vec[0]-b1vec[0]*b2vec[2];
            ABxBC[2]=b1vec[0]*b2vec[1]-b1vec[1]*b2vec[0];
            BCxCD[0]=b2vec[1]*b3vec[2]-b2vec[2]*b3vec[1];
            BCxCD[1]=b2vec[2]*b3vec[0]-b2vec[0]*b3vec[2];
            BCxCD[2]=b2vec[0]*b3vec[1]-b2vec[1]*b3vec[0];

            //the dihedral angle is then a function of the Cartesians of all four of its atoms
            dihedrals[i]=atan2((ABxBC[2]*BCxCD[0] - ABxBC[0]*BCxCD[2])*b2vec[1]/prevBond - (ABxBC[2]*BCxCD[1] - ABxBC[1]*BCxCD[2])*b2vec[0]/prevBond - (ABxBC[1]*BCxCD[0] - ABxBC[0]*BCxCD[1])*b2vec[2]/prevBond, ABxBC[0]*BCxCD[0] + ABxBC[1]*BCxCD[1] + ABxBC[2]*BCxCD[2]);
        }

    }

    //a phaseangle is formed by its difference to its "parent" dihedral
    void convert_dihedrals_to_phaseangles() {
        for(unsigned int i=0; i<dihedrals_top.size(); i++) {//for all dihedrals
            if(dihedrals_top[i][5]!=-1) { //dihedrals_top[i][5] indicates the "parent" dihedral of dihedral[i], so if the dihedral has a "parent" dihedral
                dihedrals[i]-=dihedrals[dihedrals_top[i][5]-1];//subtract from dihedral[i] the value of its "parent" dihedral
            }
        }
    }

    void convert_to_Cartesian() {

        double n[3];
        double bc[3];
        double D2[3];
        double AB[3];
        double dummy,sinus;


        x[dihedrals_top[0][0]-1][0]=root_origin_cartesian[0]; //use the stored external Cartesian atoms to specify the position of the first root atom
        x[dihedrals_top[0][0]-1][1]=root_origin_cartesian[1];
        x[dihedrals_top[0][0]-1][2]=root_origin_cartesian[2];

        x[dihedrals_top[0][1]-1][0]=root_origin_cartesian[0]+rootbond[0]*sin(root_origin_theta)*cos(root_origin_phi); //with the two external polar coordinates, as well as the first (internal) rootangle, the position of the second root atom can be calculated
        x[dihedrals_top[0][1]-1][1]=root_origin_cartesian[1]+rootbond[0]*sin(root_origin_theta)*sin(root_origin_phi);
        x[dihedrals_top[0][1]-1][2]=root_origin_cartesian[2]+rootbond[0]*cos(root_origin_theta);


        //the calculation of the third root atom is based on the NeRF algorithm of Rosetta (Rohl CA,Strauss CE,Misura KM,Baker D. Protein structure prediction using Rosetta. Methods Enzymol 2004; 383: 66-93.),
        //see e.g. Parsons J,Holmes JB,Rojas JM,Tsai J,Strauss CE. Practical conversion from torsion space to Cartesian space for in silico protein synthesis. J Comput Chem 2005; 26: 1063-1068
        //the external dihedral is formed by the three root atoms plus the origin of the coordinate system (0,0,0)
        bc[0]=x[dihedrals_top[0][1]-1][0]-x[dihedrals_top[0][0]-1][0];
        bc[1]=x[dihedrals_top[0][1]-1][1]-x[dihedrals_top[0][0]-1][1];
        bc[2]=x[dihedrals_top[0][1]-1][2]-x[dihedrals_top[0][0]-1][2];
        dummy=sqrt(bc[0]*bc[0]+bc[1]*bc[1]+bc[2]*bc[2]);
        bc[0]/=dummy;
        bc[1]/=dummy;
        bc[2]/=dummy;

        n[0]=root_origin_cartesian[1]*bc[2]-root_origin_cartesian[2]*bc[1];
        n[1]=root_origin_cartesian[2]*bc[0]-root_origin_cartesian[0]*bc[2];
        n[2]=root_origin_cartesian[0]*bc[1]-root_origin_cartesian[1]*bc[0];
        dummy=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
        n[0]/=dummy;
        n[1]/=dummy;
        n[2]/=dummy;

        sinus=sin(rootangle);
        D2[0]=rootbond[1]*cos(rootangle);
        D2[1]=rootbond[1]*cos(root_origin_dihedral)*sinus;
        D2[2]=rootbond[1]*sin(root_origin_dihedral)*sinus;

        x[dihedrals_top[0][2]-1][0]=bc[0]*D2[0] + (n[1]*bc[2] - n[2]*bc[1])*D2[1] + n[0]*D2[2] + x[dihedrals_top[0][1]-1][0];
        x[dihedrals_top[0][2]-1][1]=bc[1]*D2[0] + (n[2]*bc[0] - n[0]*bc[2])*D2[1] + n[1]*D2[2] + x[dihedrals_top[0][1]-1][1];
        x[dihedrals_top[0][2]-1][2]=bc[2]*D2[0] + (n[0]*bc[1] - n[1]*bc[0])*D2[1] + n[2]*D2[2] + x[dihedrals_top[0][1]-1][2];

        //now for all the non-root atoms also use the NeRF-based calculation
        for(unsigned int i=0; i<dihedrals_top.size(); i++) {

            bc[0]=x[dihedrals_top[i][2]-1][0]-x[dihedrals_top[i][1]-1][0];
            bc[1]=x[dihedrals_top[i][2]-1][1]-x[dihedrals_top[i][1]-1][1];
            bc[2]=x[dihedrals_top[i][2]-1][2]-x[dihedrals_top[i][1]-1][2];
            dummy=sqrt(bc[0]*bc[0]+bc[1]*bc[1]+bc[2]*bc[2]);
            bc[0]/=dummy;
            bc[1]/=dummy;
            bc[2]/=dummy;

            AB[0]=x[dihedrals_top[i][1]-1][0]-x[dihedrals_top[i][0]-1][0];
            AB[1]=x[dihedrals_top[i][1]-1][1]-x[dihedrals_top[i][0]-1][1];
            AB[2]=x[dihedrals_top[i][1]-1][2]-x[dihedrals_top[i][0]-1][2];

            n[0]=AB[1]*bc[2]-AB[2]*bc[1];
            n[1]=AB[2]*bc[0]-AB[0]*bc[2];
            n[2]=AB[0]*bc[1]-AB[1]*bc[0];

            dummy=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
            n[0]/=dummy;
            n[1]/=dummy;
            n[2]/=dummy;

            sinus=sin(angles[i]);
            D2[0]=bonds[i]*cos(angles[i]);
            D2[1]=bonds[i]*cos(dihedrals[i])*sinus;
            D2[2]=bonds[i]*sin(dihedrals[i])*sinus;

            x[dihedrals_top[i][3]-1][0]=bc[0]*D2[0] + (n[1]*bc[2] - n[2]*bc[1])*D2[1] + n[0]*D2[2] + x[dihedrals_top[i][2]-1][0];
            x[dihedrals_top[i][3]-1][1]=bc[1]*D2[0] + (n[2]*bc[0] - n[0]*bc[2])*D2[1] + n[1]*D2[2] + x[dihedrals_top[i][2]-1][1];
            x[dihedrals_top[i][3]-1][2]=bc[2]*D2[0] + (n[0]*bc[1] - n[1]*bc[0])*D2[1] + n[2]*D2[2] + x[dihedrals_top[i][2]-1][2];


        }
    }

    //to convert all phaseangles back to dihedrals
    void convert_phaseangles_to_dihedrals() {
        for(unsigned int i=0; i<dihedrals_top.size(); i++) {//for all dihedrals
            if(dihedrals_top[i][5]!=-1) {//check if this dihedral was implemented as a phase angle
                dihedrals[i]+=dihedrals[dihedrals_top[i][5]-1]; //if so, add the value of the "parent" dihedral back up
            }
        }
    }

};





