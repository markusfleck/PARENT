//    The class code for handling an EntropyMatrix from a .par file of the PARENT suite 
//    Copyright (C) 2016  Markus Fleck (member of the laboratory of Bojan Zagrovic, University of Vienna)
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


#pragma pack(1)

using namespace std;


EntropyMatrix::EntropyMatrix(char const * infileInput){
  infile.exceptions(std::ifstream::badbit);
  infile.open(infileInput, ios::binary | ios::out);//open the .par file;
  if(!infile.good())
    {
      MyError myError((string("ERROR OPENING FILE ")+string(infileInput)+string("! ABORTING.\n")).c_str());
      throw myError;
    }
  infile.exceptions(std::ifstream::badbit | std::ifstream::failbit | std::ifstream::eofbit);
  
  
  try{
    read_PAR_header();
  }
  catch(MyError myError)
  {
    throw myError;
  }
  catch(...){
      MyError myError((string("ERROR WHILE READING THE HEADER OF THE FILE ")+string(infileInput)+string("! ABORTING.")).c_str());
      throw myError;    
  }
  
  nBonds=nDihedrals+2;
  nAngles=nDihedrals+1;

  try{
    read_PAR_body();
  }
  catch(MyError myError)
  {
    throw myError;
  }
  catch(...){
      MyError myError((string("ERROR WHILE READING THE BODY OF THE FILE ")+string(infileInput)+string("! ABORTING.")).c_str());
      throw myError;    
  }
  
}



//to create a blank EntropyMatrix with nAtoms atoms
EntropyMatrix::EntropyMatrix(int nAtoms){
  
  nBonds=nAtoms-1;
  nAngles=nAtoms-2;
  nDihedrals=nAtoms-3;


  bondsEntropy1D = new double[nBonds]; // allocate storage for reading the .par file
  anglesEntropy1D=new double[nAngles];
  dihedralsEntropy1D=new double[nDihedrals];
  bbEntropy=new double[nBonds*(nBonds-1)/2];
  baEntropy=new double[nBonds*nAngles];
  bdEntropy=new double[nBonds*nDihedrals];
  aaEntropy=new double[nAngles*(nAngles-1)/2];
  adEntropy=new double[nAngles*nDihedrals];
  ddEntropy=new double[nDihedrals*(nDihedrals-1)/2];

  if(!((bondsEntropy1D!=NULL)&&(anglesEntropy1D!=NULL)&&(dihedralsEntropy1D!=NULL)&&(bbEntropy!=NULL)&&(baEntropy!=NULL)&&(bdEntropy!=NULL)&&(aaEntropy!=NULL)&&(adEntropy!=NULL)&&(ddEntropy!=NULL))) {
      MyError myError("ERROR: ALLOCATION FOR PAR FILE BODY FAILED. MORE MEMORY NEEDED?");
      throw myError;
  }
  
  double_prec=-1; //set the precision of the .bat file which underlies this .par file (EntropyMatrix) to "unknown"
  numFrames=0; //initialize the number of frames of the underlying .bat trajectory to 0, indicating unknown
  version=4; //probably unnecessary
  bDens1D=0; // and initialized all the used bins with 0
  aDens1D=0;
  dDens1D=0;
  bDens=0;
  aDens=0;
  dDens=0; 
  
  vector <int> dummyintvec;
  for(int i=0;i<5;i++){dummyintvec.push_back(0);}
  for(int i=0;i<nDihedrals;i++){dihedrals_top.push_back(dummyintvec);}
  for(int i=0;i<nAtoms;i++){
    masses.push_back(0);
    residues.push_back(string(""));
    residueNumbers.push_back(0);
    atomNames.push_back(string(""));
    belongsToMolecule.push_back(string(""));
    }
}


EntropyMatrix::~EntropyMatrix(){
  
  delete[] bondsEntropy1D; //deallocate storage for reading the .par file
  delete[] anglesEntropy1D;
  delete[] dihedralsEntropy1D;
  delete[] bbEntropy;
  delete[] baEntropy;
  delete[] bdEntropy;
  delete[] aaEntropy;
  delete[] adEntropy;
  delete[] ddEntropy;
  
  if(infile.is_open()){
		infile.close();
	}
	if(outfile.is_open()){
		outfile.close();
	}
}


void EntropyMatrix::write(char const * outfileInput){
  outfile.exceptions(std::ofstream::badbit);
  outfile.open(outfileInput, ios::binary | ios::out);//open the .par file
	if(!outfile.good())
    {
      MyError myError((string("ERROR OPENING FILE ")+string(outfileInput)+string("! ABORTING.\n")).c_str());
      throw myError;
    }
    outfile.exceptions(std::ofstream::badbit | std::ofstream::failbit | std::ofstream::eofbit);
  
  
   
  try{
    write_PAR_header();
  }
  catch(MyError myError)
  {
    throw myError;
  }
  catch(...){
      MyError myError((string("ERROR WHILE WRITING THE HEADER OF THE FILE ")+string(outfileInput)+string("! ABORTING.")).c_str());
      throw myError;    
  }
  
    
  try{
    write_PAR_body();
  }
  catch(MyError myError)
  {
    throw myError;
  }
  catch(...){
    MyError myError((string("ERROR WRITING THE BODY OF THE FILE ")+string(outfileInput)+string("! ABORTING.")).c_str());
    throw myError;    
  }
	
		outfile.close();
}



//to get 1D entropies 
double EntropyMatrix::getEntropy(int type, int index) {
    if(index>0){
			if((type==TYPE_B)&&(index<=nBonds)) {
					return bondsEntropy1D[index-1];
			}
			if((type==TYPE_A)&&(index<=nAngles)) {
					return anglesEntropy1D[index-1];
			}
			if((type==TYPE_D)&&(index<=nDihedrals)) {
					return dihedralsEntropy1D[index-1];
			}
		}
    return 0;
}


//to get 2D entropy values 
double EntropyMatrix::get2DEntropy(int type1,int type2,int index1, int index2) {
    int smaller,bigger,index;
    index1--;
    index2--;
    if((index1>=0)&&(index2>=0)){        
			if((type1==TYPE_B)&&(type2==TYPE_B)&&(index1!=index2)&&(index1<nBonds)&&(index2<nBonds)) { // for 2D entropy of two different bonds
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nBonds-smaller)*(nBonds-smaller-1)/2+smaller-bigger;//the 2D-entropies bonds-bonds (also angles-angles and dihedrals-dihedrals) were stored in reverse order, as documented in "Parent.cpp"
					return bbEntropy[index];
			}
			if((type1==TYPE_B)&&(type2==TYPE_A)&&(index1<nBonds)&&(index2<nAngles)) { // for 2D entropy of a bond and an angle
					return baEntropy[index1*nAngles+index2];
			}
			if((type1==TYPE_A)&&(type2==TYPE_B)&&(index1<nAngles)&&(index2<nBonds)) { // for 2D entropy of an angle and a bond
					return baEntropy[index2*nAngles+index1];
			}
			if((type1==TYPE_B)&&(type2==TYPE_D)&&(index1<nBonds)&&(index2<nDihedrals)) { // for 2D entropy of a bond and a dihedral
					return bdEntropy[index1*nDihedrals+index2];
			}
			if((type1==TYPE_D)&&(type2==TYPE_B)&&(index1<nDihedrals)&&(index2<nBonds)) { // for 2D entropy of a dihedral and a bond
					return bdEntropy[index2*nDihedrals+index1];
			}
			if((type1==TYPE_A)&&(type2==TYPE_A)&&(index1!=index2)&&(index1<nAngles)&&(index2<nAngles)) { // for 2D entropy of two different angles
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nAngles-smaller)*(nAngles-smaller-1)/2+smaller-bigger;//the 2D-entropies angles-angles were stored in reverse order, as documented in "Parent.cpp"
					return aaEntropy[index];
			}
			if((type1==TYPE_A)&&(type2==TYPE_D)&&(index1<nAngles)&&(index2<nDihedrals)) { // for 2D entropy of a bond and an angle
					return adEntropy[index1*nDihedrals+index2];
			}
			if((type1==TYPE_D)&&(type2==TYPE_A)&&(index1<nDihedrals)&&(index2<nAngles)) { // for 2D entropy of an angle and a bond
					return adEntropy[index2*nDihedrals+index1];
			}
			if((type1==TYPE_D)&&(type2==TYPE_D)&&(index1!=index2)&&(index1<nDihedrals)&&(index2<nDihedrals)) { // for 2D entropy of two different dihedrals
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nDihedrals-smaller)*(nDihedrals-smaller-1)/2+smaller-bigger;//the 2D-entropies dihedrals-dihedrals were stored in reverse order, as documented in "Parent.cpp"
					return ddEntropy[index];
			}
		}
		return 0;
}


//to get mutual information values
double EntropyMatrix::getMutual(int type1,int type2,int index1, int index2) {
    int smaller,bigger,index;
    index1--;
    index2--;
    if((index1>=0)&&(index2>=0)){        
			if((type1==TYPE_B)&&(type2==TYPE_B)&&(index1!=index2)&&(index1<nBonds)&&(index2<nBonds)) { // for mutual information between two different bonds
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nBonds-smaller)*(nBonds-smaller-1)/2+smaller-bigger;//the 2D-entropies bonds-bonds (also angles-angles and dihedrals-dihedrals) were stored in reverse order, as documented in "Parent.cpp"
					return bondsEntropy1D[smaller]+bondsEntropy1D[bigger]-bbEntropy[index];
			}
			if((type1==TYPE_B)&&(type2==TYPE_A)&&(index1<nBonds)&&(index2<nAngles)) { // for mutual information between a bond and an angle
					return bondsEntropy1D[index1]+anglesEntropy1D[index2]-baEntropy[index1*nAngles+index2];
			}
			if((type1==TYPE_A)&&(type2==TYPE_B)&&(index1<nAngles)&&(index2<nBonds)) { // for mutual information between an angle and a bond
					return bondsEntropy1D[index2]+anglesEntropy1D[index1]-baEntropy[index2*nAngles+index1];
			}
			if((type1==TYPE_B)&&(type2==TYPE_D)&&(index1<nBonds)&&(index2<nDihedrals)) { // for mutual information between a bond and a dihedral
					return bondsEntropy1D[index1]+dihedralsEntropy1D[index2]-bdEntropy[index1*nDihedrals+index2];
			}
			if((type1==TYPE_D)&&(type2==TYPE_B)&&(index1<nDihedrals)&&(index2<nBonds)) { // for mutual information between a dihedral and a bond
					return bondsEntropy1D[index2]+dihedralsEntropy1D[index1]-bdEntropy[index2*nDihedrals+index1];
			}
			if((type1==TYPE_A)&&(type2==TYPE_A)&&(index1!=index2)&&(index1<nAngles)&&(index2<nAngles)) { // for mutual information between two different angles
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nAngles-smaller)*(nAngles-smaller-1)/2+smaller-bigger;//the 2D-entropies angles-angles were stored in reverse order, as documented in "Parent.cpp"
					return anglesEntropy1D[smaller]+anglesEntropy1D[bigger]-aaEntropy[index];
			}
			if((type1==TYPE_A)&&(type2==TYPE_D)&&(index1<nAngles)&&(index2<nDihedrals)) { // for mutual information between a bond and an angle
					return anglesEntropy1D[index1]+dihedralsEntropy1D[index2]-adEntropy[index1*nDihedrals+index2];
			}
			if((type1==TYPE_D)&&(type2==TYPE_A)&&(index1<nDihedrals)&&(index2<nAngles)) { // for mutual information between an angle and a bond
					return anglesEntropy1D[index2]+dihedralsEntropy1D[index1]-adEntropy[index2*nDihedrals+index1];
			}
			if((type1==TYPE_D)&&(type2==TYPE_D)&&(index1!=index2)&&(index1<nDihedrals)&&(index2<nDihedrals)) { // for mutual information between two different dihedrals
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nDihedrals-smaller)*(nDihedrals-smaller-1)/2+smaller-bigger;//the 2D-entropies dihedrals-dihedrals were stored in reverse order, as documented in "Parent.cpp"
					return dihedralsEntropy1D[smaller]+dihedralsEntropy1D[bigger]-ddEntropy[index];
			}
		}
	return 0;
}

//to set 1D entropies 
void EntropyMatrix::setEntropy(int type, int index, double value) {
    if(index>0){
			if((type==TYPE_B)&&(index<=nBonds)) {
					bondsEntropy1D[index-1]=value;
			}
			if((type==TYPE_A)&&(index<=nAngles)) {
					anglesEntropy1D[index-1]=value;
			}
			if((type==TYPE_D)&&(index<=nDihedrals)) {
					dihedralsEntropy1D[index-1]=value;
			}
	}
}

//to set 2D entropy values 
void EntropyMatrix::set2DEntropy(int type1,int type2,int index1, int index2, double value) {
    int smaller,bigger,index;
    index1--;
    index2--;
    if((index1>=0)&&(index2>=0)&&(value>=0)){    
			if((type1==TYPE_B)&&(type2==TYPE_B)&&(index1!=index2)&&(index1<nBonds)&&(index2<nBonds)) { // for 2D entropy of two different bonds
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nBonds-smaller)*(nBonds-smaller-1)/2+smaller-bigger;//the 2D-entropies bonds-bonds (also angles-angles and dihedrals-dihedrals) are stored in reverse order, as documented in "Parent.cpp"
					bbEntropy[index]=value;
			}
			if((type1==TYPE_B)&&(type2==TYPE_A)&&(index1<nBonds)&&(index2<nAngles)) { // for 2D entropy of a bond and an angle
					baEntropy[index1*nAngles+index2]=value;
			}
			if((type1==TYPE_A)&&(type2==TYPE_B)&&(index1<nAngles)&&(index2<nBonds)) { // for 2D entropy of an angle and a bond
					baEntropy[index2*nAngles+index1]=value;
			}
			if((type1==TYPE_B)&&(type2==TYPE_D)&&(index1<nBonds)&&(index2<nDihedrals)) { // for 2D entropy of a bond and a dihedral
					bdEntropy[index1*nDihedrals+index2]=value;
			}
			if((type1==TYPE_D)&&(type2==TYPE_B)&&(index1<nDihedrals)&&(index2<nBonds)) { // for 2D entropy of a dihedral and a bond
					bdEntropy[index2*nDihedrals+index1]=value;
			}
			if((type1==TYPE_A)&&(type2==TYPE_A)&&(index1!=index2)&&(index1<nAngles)&&(index2<nAngles)) { // for 2D entropy of two different angles
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nAngles-smaller)*(nAngles-smaller-1)/2+smaller-bigger;//the 2D-entropies angles-angles are stored in reverse order, as documented in "Parent.cpp"
					aaEntropy[index]=value;
			}
			if((type1==TYPE_A)&&(type2==TYPE_D)&&(index1<nAngles)&&(index2<nDihedrals)) { // for 2D entropy of a bond and an angle
					adEntropy[index1*nDihedrals+index2]=value;
			}
			if((type1==TYPE_D)&&(type2==TYPE_A)&&(index1<nDihedrals)&&(index2<nAngles)) { // for 2D entropy of an angle and a bond
					adEntropy[index2*nDihedrals+index1]=value;
			}
			if((type1==TYPE_D)&&(type2==TYPE_D)&&(index1!=index2)&&(index1<nDihedrals)&&(index2<nDihedrals)) { // for 2D entropy of two different dihedrals
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nDihedrals-smaller)*(nDihedrals-smaller-1)/2+smaller-bigger;//the 2D-entropies dihedrals-dihedrals are stored in reverse order, as documented in "Parent.cpp"
					ddEntropy[index]=value;
			}
		}
}


//to set mutual information values by changing 2D entropy values, without modyfing 2D entropy values
void EntropyMatrix::setMutual(int type1,int type2,int index1, int index2, double value) {
    int smaller,bigger,index;
    index1--;
    index2--;
    if((index1>=0)&&(index2>=0)&&(value>=0)){
			if((type1==TYPE_B)&&(type2==TYPE_B)&&(index1!=index2)&&(index1<nBonds)&&(index2<nBonds)) { // for mutual information between two different bonds
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nBonds-smaller)*(nBonds-smaller-1)/2+smaller-bigger;//the 2D-entropies bonds-bonds (also angles-angles and dihedrals-dihedrals) were stored in reverse order, as documented in "Parent.cpp"
					bbEntropy[index]=bondsEntropy1D[smaller]+bondsEntropy1D[bigger]-value;
			}
			if((type1==TYPE_B)&&(type2==TYPE_A)&&(index1<nBonds)&&(index2<nAngles)) { // for mutual information between a bond and an angle
					baEntropy[index1*nAngles+index2]=bondsEntropy1D[index1]+anglesEntropy1D[index2]-value;
			}
			if((type1==TYPE_A)&&(type2==TYPE_B)&&(index1<nAngles)&&(index2<nBonds)) { // for mutual information between an angle and a bond
					baEntropy[index2*nAngles+index1]=bondsEntropy1D[index2]+anglesEntropy1D[index1]-value;
			}
			if((type1==TYPE_B)&&(type2==TYPE_D)&&(index1<nBonds)&&(index2<nDihedrals)) { // for mutual information between a bond and a dihedral
					bdEntropy[index1*nDihedrals+index2]=bondsEntropy1D[index1]+dihedralsEntropy1D[index2]-value;
			}
			if((type1==TYPE_D)&&(type2==TYPE_B)&&(index1<nDihedrals)&&(index2<nBonds)) { // for mutual information between a dihedral and a bond
					bdEntropy[index2*nDihedrals+index1]=bondsEntropy1D[index2]+dihedralsEntropy1D[index1]-value;
			}
			if((type1==TYPE_A)&&(type2==TYPE_A)&&(index1!=index2)&&(index1<nAngles)&&(index2<nAngles)) { // for mutual information between two different angles
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nAngles-smaller)*(nAngles-smaller-1)/2+smaller-bigger;//the 2D-entropies angles-angles were stored in reverse order, as documented in "Parent.cpp"
					aaEntropy[index]=anglesEntropy1D[smaller]+anglesEntropy1D[bigger]-value;
			}
			if((type1==TYPE_A)&&(type2==TYPE_D)&&(index1<nAngles)&&(index2<nDihedrals)) { // for mutual information between a bond and an angle
					adEntropy[index1*nDihedrals+index2]=anglesEntropy1D[index1]+dihedralsEntropy1D[index2]-value;
			}
			if((type1==TYPE_D)&&(type2==TYPE_A)&&(index1<nDihedrals)&&(index2<nAngles)) { // for mutual information between an angle and a bond
					adEntropy[index2*nDihedrals+index1]=anglesEntropy1D[index2]+dihedralsEntropy1D[index1]-value;
			}
			if((type1==TYPE_D)&&(type2==TYPE_D)&&(index1!=index2)&&(index1<nDihedrals)&&(index2<nDihedrals)) { // for mutual information between two different dihedrals
					smaller=index1<index2?index1:index2;
					bigger=index1<index2?index2:index1;
					index=(nDihedrals-smaller)*(nDihedrals-smaller-1)/2+smaller-bigger;//the 2D-entropies dihedrals-dihedrals were stored in reverse order, as documented in "Parent.cpp"
					ddEntropy[index]=dihedralsEntropy1D[smaller]+dihedralsEntropy1D[bigger]-value;
			}
		}
}

void EntropyMatrix::write_PAR_header() {
    int dummy=dihedrals_top.size();
    int myversion=4;
    char dummystring[31];


    outfile.write((char*)&myversion, sizeof(int)); //first write the .par version number as an integer
    outfile.write((char*)&double_prec, sizeof(int)); //then write an integer declaring if the trajectory was stored in double precision (-1 means unknown)
    outfile.write((char*)&dummy, sizeof(int)); //write an integer containing the number of dihedrals
    outfile.write((char*)&numFrames, sizeof(int)); //and an integer containing the number of frames of the trajectory used for calculation
    outfile.write((char*)&bDens1D, sizeof(int)); //write the used number of bins for 1D histograms
    outfile.write((char*)&aDens1D, sizeof(int));
    outfile.write((char*)&dDens1D, sizeof(int));
    outfile.write((char*)&bDens, sizeof(int));//and the used number of bins for 2D histograms
    outfile.write((char*)&aDens, sizeof(int));
    outfile.write((char*)&dDens, sizeof(int));

    for(int i=0; i<dummy+3; i++) { //for ever atom in the system
        for(int j=0; j<8; j++) {
            dummystring[j]='\0';
        }
        sprintf(dummystring, "%.7s", residues[i].c_str());
        outfile.write(dummystring, 8*sizeof(char));//write the name of the residue it belongs to
        outfile.write((char*)&(residueNumbers[i]), sizeof(int));//write the number of the residue it belongs to
        for(int j=0; j<8; j++) {
            dummystring[j]='\0';
        }
        sprintf(dummystring, "%.7s", atomNames[i].c_str());
        outfile.write(dummystring, 8*sizeof(char));//write the name of the atom
        for(int j=0; j<31; j++) {
            dummystring[j]='\0';
        }
        sprintf(dummystring, "%.30s", belongsToMolecule[i].c_str());
        outfile.write(dummystring, 31*sizeof(char));//write the name of the molecule it belongs to
    }


    for(int i=0; i<dummy; i++) {//then for all dihedrals
        outfile.write((char*)&(dihedrals_top[i][0]), sizeof(int)); //write the atomnumber of the first atom
        outfile.write((char*)&(dihedrals_top[i][1]), sizeof(int)); //second atom
        outfile.write((char*)&(dihedrals_top[i][2]), sizeof(int)); //third atom
        outfile.write((char*)&(dihedrals_top[i][3]), sizeof(int)); //fourth atom
        outfile.write((char*)&(dihedrals_top[i][4]), sizeof(int)); // an integer containing the type of the dihedral(physical=0,pseudo=1,improper=-1)
        outfile.write((char*)&(dihedrals_top[i][5]), sizeof(int)); // and an integer containing the "parent"dihedral(for phaseangles, -1 if no "parent")
    }
    for(int i=0; i<dummy+3; i++) {
        outfile.write((char*)&(masses[i]), sizeof(float)); //then write the whole massestor of the atoms in the system in single precision (float)
    }
}



void EntropyMatrix::read_PAR_header() {
    char dummystring[31];

    infile.read((char*)&version, sizeof(int)); //first read the .par version number as an integer
    infile.read((char*)&double_prec, sizeof(int));//then read an integer declaring if the trajectory was stored in double precision (-1 means unknown)
    infile.read((char*)&nDihedrals,sizeof(int));//read an integer containing the number of dihedrals
    infile.read((char*)&numFrames,sizeof(int));//and an integer containing the number of frames of the trajectory used for calculation
    infile.read((char*)&bDens1D,sizeof(int)); //read the used number of bins for 1D histograms
    infile.read((char*)&aDens1D,sizeof(int));
    infile.read((char*)&dDens1D,sizeof(int));
    infile.read((char*)&bDens,sizeof(int)); //and the used number of bins for 2D histograms
    infile.read((char*)&aDens,sizeof(int));
    infile.read((char*)&dDens,sizeof(int));


    if(version<0) {
        MyError myError((string("ERROR: FILE HEADER CORRUPTED. VERSION NUMBER (")+to_string(version)+string(") < 0! ABORTING.")).c_str());
        throw myError;
    }
    if((double_prec!=0)&&(double_prec!=1)&&(double_prec!=-1)) {
        MyError myError((string("ERROR: FILE HEADER CORRUPTED. DOUBLE PRECISION VALUE (")+to_string(double_prec)+string(") NEITHER 0, 1 NOR -1! ABORTING.")).c_str());
        throw myError;
    }
    if(nDihedrals>19997) {
        cerr<<"WARNING: "<<nDihedrals+3<<" ATOMS DECLARED IN THE FILE HEADER (CORRUPTED?). THIS WILL LEAD TO LARGE OUTPUT."<<endl;
    }
    if(nDihedrals<0) {
        MyError myError((string("ERROR: FILE HEADER CORRUPTED. NUMBER OF DIHEDRALS (")+to_string(nDihedrals)+string(") < 0! ABORTING.")).c_str());
        throw myError;
    }
    if(numFrames<0) {
        MyError myError((string("ERROR: FILE HEADER CORRUPTED. NUMBER OF FRAMES (")+to_string(numFrames)+string(") < 0! ABORTING.")).c_str());
        throw myError;
    }
    if(bDens1D<1) {
        MyError myError((string("ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR BONDS IN 1D (")+to_string(bDens1D)+string(") < 1! ABORTING.")).c_str());
        throw myError;
    }
    if(aDens1D<1) {
        MyError myError((string("ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR ANGLES IN 1D (")+to_string(aDens1D)+string(") < 1! ABORTING.")).c_str());
        throw myError;
    }
    if(dDens1D<1) {
        MyError myError((string("ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR DIHEDRALS IN 1D (")+to_string(dDens1D)+string(") < 1! ABORTING.")).c_str());
        throw myError;
    }
    if(bDens<1) {
        MyError myError((string("ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR BONDS IN 2D (")+to_string(bDens)+string(") < 1! ABORTING.")).c_str());
        throw myError;
    }
    if(aDens<1) {
        MyError myError((string("ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR ANGLES IN 2D (")+to_string(aDens)+string(") < 1! ABORTING.")).c_str());
        throw myError;
    }
    if(dDens<1) {
        MyError myError((string("ERROR: FILE HEADER CORRUPTED. NUMBER BINS FOR DIHEDRALS IN 2D (")+to_string(dDens)+string(") < 1! ABORTING.")).c_str());
        throw myError;
    }
		
		
		
		if (version>=3) {
        for(int i=0; i<nDihedrals+3; i++) { //for every atom in the system
            infile.read(dummystring, 8*sizeof(char));//read the name of the residue it belongs to
            residues.push_back(dummystring);
            residueNumbers.push_back(0);//the number of the residue it belongs to
            infile.read((char*)&(residueNumbers[i]), sizeof(float));
            infile.read(dummystring, 8*sizeof(char));//the name of the atom
            atomNames.push_back(dummystring);
            infile.read(dummystring, 31*sizeof(char));
            belongsToMolecule.push_back(dummystring); //and the name of the molecule it belongs to
        }
    }
		


    vector<int>dummyvec;
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    dummyvec.push_back(0);
    for(int i=0; i<nDihedrals; i++) { //then for all dihedrals
        dihedrals_top.push_back(dummyvec);
        infile.read((char*)&(dihedrals_top[i][0]),sizeof(int));//write the atomnumber of the first atom
        if(dihedrals_top[i][0]<1) {
          MyError myError("ERROR: FILE HEADER CORRUPTED. DIHEDRAL ID < 1 !");
          throw myError;          
        }
        infile.read((char*)&(dihedrals_top[i][1]),sizeof(int));//second atom
        if(dihedrals_top[i][1]<1) {
          MyError myError("ERROR: FILE HEADER CORRUPTED. DIHEDRAL ID < 1 !");
          throw myError;   
        }
        infile.read((char*)&(dihedrals_top[i][2]),sizeof(int));//third atom
        if(dihedrals_top[i][2]<1) {
          MyError myError("ERROR: FILE HEADER CORRUPTED. DIHEDRAL ID < 1 !");
          throw myError;   
        }
        infile.read((char*)&(dihedrals_top[i][3]),sizeof(int));//fourth atom
        if(dihedrals_top[i][3]<1) {
          MyError myError("ERROR: FILE HEADER CORRUPTED. DIHEDRAL ID < 1 !");
          throw myError;   
        }
        infile.read((char*)&(dihedrals_top[i][4]),sizeof(int));// an integer containing the type of the dihedral(physical=0,pseudo=1,improper=-1)
        if(!((dihedrals_top[i][4]==1)||(dihedrals_top[i][4]==0)||(dihedrals_top[i][4]==-1))) {
          MyError myError("ERROR: FILE HEADER CORRUPTED. DIHEDRAL TYPE NEITHER 0, 1 NOR -1 !");
          throw myError;   
        }
        infile.read((char*)&(dihedrals_top[i][5]),sizeof(int));// and an integer containing the "parent" dihedral(for phaseangles, -1 if no "parent")
        if(dihedrals_top[i][5]<-1) {
          MyError myError("ERROR: FILE HEADER CORRUPTED. DIHEDRAL IS PHASEANGLE OF A DIHEDRAL WITH NEGATIVE ID !");
          throw myError;   
        }

        if((version==1)&&(dihedrals_top[i][5]!=-1)) {
            dihedrals_top[i][5]++; //for backwards compatibility
        }
    }
    for(int i=0; i<nDihedrals+3; i++) {
        masses.push_back(0);
        infile.read((char*)&(masses)[i], sizeof(float));//then read the whole massestor of the atoms in the system in single precision (float)
    }

}

void EntropyMatrix::write_PAR_body() {

    int nBonds=nDihedrals+2;
    int nAngles=nDihedrals+1;

    outfile.write((char*)bondsEntropy1D, nBonds*sizeof(double)); //write the 1D entropy bonds array
    outfile.write((char*)anglesEntropy1D, nAngles*sizeof(double)); //write the 1D entropy angles array
    outfile.write((char*)dihedralsEntropy1D, nDihedrals*sizeof(double)); //write the 1D entropy dihedrals array
    outfile.write((char*)bbEntropy, nBonds*(nBonds-1)/2*sizeof(double)); //write the 2D bonds-bonds half-matrix as an array
    outfile.write((char*)baEntropy, nBonds*nAngles*sizeof(double)); //write the 2D bonds-angles matrix
    outfile.write((char*)bdEntropy, nBonds*nDihedrals*sizeof(double)); //write the 2D bonds-dihedrals matrix
    outfile.write((char*)aaEntropy, nAngles*(nAngles-1)/2*sizeof(double)); //write the 2D angles-angles half-matrix as an array
    outfile.write((char*)adEntropy, nAngles*nDihedrals*sizeof(double)); //write the 2D angles-dihedrals matrix
    outfile.write((char*)ddEntropy, nDihedrals*(nDihedrals-1)/2*sizeof(double)); //write the 2D dihedrlas-dihedrals half-matrix as an array
}

void EntropyMatrix::read_PAR_body() {

    int nBonds=nDihedrals+2;
    int nAngles=nDihedrals+1;

    bondsEntropy1D = new double[nBonds]; // allocate storage for reading the .par file
    anglesEntropy1D=new double[nAngles];
    dihedralsEntropy1D=new double[nDihedrals];
    bbEntropy=new double[nBonds*(nBonds-1)/2];
    baEntropy=new double[nBonds*nAngles];
    bdEntropy=new double[nBonds*nDihedrals];
    aaEntropy=new double[nAngles*(nAngles-1)/2];
    adEntropy=new double[nAngles*nDihedrals];
    ddEntropy=new double[nDihedrals*(nDihedrals-1)/2];

    if(!((bondsEntropy1D!=NULL)&&(anglesEntropy1D!=NULL)&&(dihedralsEntropy1D!=NULL)&&(bbEntropy!=NULL)&&(baEntropy!=NULL)&&(bdEntropy!=NULL)&&(aaEntropy!=NULL)&&(adEntropy!=NULL)&&(ddEntropy!=NULL))) {
        MyError myError("ERROR: ALLOCATION FOR PAR FILE BODY FAILED. MORE MEMORY NEEDED?");
        throw myError;
    }

    infile.read((char*)bondsEntropy1D, nBonds*sizeof(double)); // and the actually reads it in (see write_PAR_body for details)
    infile.read((char*)anglesEntropy1D, nAngles*sizeof(double));
    infile.read((char*)dihedralsEntropy1D, nDihedrals*sizeof(double));
    infile.read((char*)bbEntropy, nBonds*(nBonds-1)/2*sizeof(double));
    infile.read((char*)baEntropy, nBonds*nAngles*sizeof(double));
    infile.read((char*)bdEntropy, nBonds*nDihedrals*sizeof(double));
    infile.read((char*)aaEntropy, nAngles*(nAngles-1)/2*sizeof(double));
    infile.read((char*)adEntropy, nAngles*nDihedrals*sizeof(double));
    infile.read((char*)ddEntropy, nDihedrals*(nDihedrals-1)/2*sizeof(double));
}

int EntropyMatrix::getNBonds(){
  return nBonds;
}


int EntropyMatrix::getNAngles(){
  return nAngles;
}

int EntropyMatrix::getNDihedrals(){
  return nDihedrals;
}


#pragma pack(0)

