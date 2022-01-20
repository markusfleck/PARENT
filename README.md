# PARENT: a parallel software suite for the calculation of configurational entropy in biomolecular systems

PARENT suite, a collection of programs to compute the configurational entropy from a molecular dynamics trajectory
Copyright (C) 2015  Markus Fleck (member of the laboratory of Bojan Zagrovic, University of Vienna) 

These programs are free software: you can redistribute them and/or modify
them under the terms of the GNU General Public License version 3 
as published by the Free Software Foundation.

These programs are distributed in the hope that they will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

A scientific publication about this program has been released in the Journal of Chemical Theory and Computation:  
"PARENT: A Parallel Software Suite for the Calculation of Configurational Entropy in Biomolecular Systems"  
DOI: 10.1021/acs.jctc.5b01217  

We kindly ask you to include this citation in works that publish
results generated using this program or modifications of it.
<br />  
<br />  
<br />  
0) QUICK AND DIRTY

  Unzip, and navigate to the top folder. In the file "parameters" change the input
  files (TRJ,TOP) according to the trajectory you want to compute the configurational 
  entropy for. Make sure you have a recent version of Open MPI and gcc installed on
  your Linux system. 
  For example, if you want to use a single CPU with four cores for the calculation, issue 
  
  	export OMP_NUM_THREADS=4; ./run.sh parameters
  
  in a bash shell. The output of the configurational entropy terms using the MIST approximation 
  is contained in "output/PARENT_MIST.txt". If you are looking for the pairwise mutual
  information terms, they are in "output/PARENT_MIE.txt", with their indexing in
  "output/PARENT_MIE_topology.txt".
  
  If you run the code for testing purposes without any modifications,
  the calculation should take something like 2-5 minutes. The resulting values at the end of the
  files should match those in "test_system/sample_output".
	 
  If you get permission errors that means that your program has resided (or still resides) on a
  filesystem which does not support Linux permissions. Move to a supported filesystem and 
  unzip again or "chmod +x" the files which throw errors.
	
  If you get warnings about a node not being able to bind a process and performance is not as expected, 
  installing numactl (maybe devel) packages and recompiling/installing MPI should help (for
  using Docker, as outlined below, this behavior is expected and should be fine on a single CPU, which is
  what the Dockerfiles are intended for).
	
  You are strongly encouraged to read the rest of this document, but at least the next section.





## 1) INSTALLATION AND TESTING

### 2.1) System installation

  This code uses MPI as well as openMP parallelization, so you need to make
  sure that your system supports this. The code was developed and tested 
  in a GNU/Linux environment using Open MPI version 1.7.4 as a MPI implementation 
  and gcc version 4.8.2 as a compiler, but different versions/compilers should work 
  just as well given my experience over the years.
  
  Sample trajectory and topology files are shipped with this package.
  To check out if your system correctly compiles and runs all of the provided
  programs, proceed as follows:

  If, for example, you want to test this code on a machine with a single 
  processor offering 8 CPU threads, from a bash shell run:
  
    export OMP_NUM_THREADS=8; ./run.sh parameters
    
  Otherwise, if you want to use the code with a batch system (job scheduler)
  in e.g. a cluster environment, just use your custom submission script
  and in it issue the above command, with OMP_NUM_THREADS set to the number of CPU threads per
  node. The run.sh script invokes mpirun via "--map-by ppr:1:node:pe=$OMP_NUM_THREADS".
  For example, if you want to calculate on 4 nodes with 8 physical cores and 16 virtual cores (i.e.  16 CPU threads),
  you would issue "export OMP_NUM_THREADS=16; ./run.sh parameters". The according "--map-by ppr:1:node:pe=16" tells 
  MPI to spawn one MPI process on each node, with 16 openMP threads each. If the 4 physical cores reside on two 
  CPUs in each node, you could change to mapping per socket, i.e. use "export OMP_NUM_THREADS=8" and, in run.sh,
  "--map-by ppr:1:socket:pe=$OMP_NUM_THREADS". Alternatively, you could spawn two MPI processes per CPU
  socket: "--map-by ppr:2:socket:pe=4". Certainly, many combinations are possibly, and you could (and should) 
  check for performance gains. From my personal experience with my machines, one MPI process per node/socket is a safe 
  bet. If you happen to use the TORQUE Resource Manager, a sample 
  submission script called "run.pbs" is supplied in the top directory.
  The script is set up for 4 nodes and 16 cores. You might want to modify these numbers
  to fit your environment.
  
  Furthermore, depending on your MPI version and cluster architecture,
  in the file "run.sh" you might want to change this part of the two MPI invocation lines:
  
    mpirun --use-hwthread-cpus --bind-to hwthread --mca btl vader,self,tcp 

  "--use-hwthread-cpus --bind-to hwthread" tells MPI to consider virtual rather than physical cores.
  During my benchmarking, virtual cores showed almost the same computational power than physical cores,
  leading to a speedup of a factor of two (each physical core features two virutal cores), 
  so using virtual cores should be the go-to approach in almost all but
  the most exotic CPU architectures. Next, if you are using InfiniBand, you might want to change the tcp part  
  to openib (e. g. to --mca btl vader,self,openib). If you don't have "vader", you could try with "sm" instead.
  
  After executing the "run.sh" script, issue the following commands to check if 
  everything went correctly:
  
    tail -n 16 output/*_MIE.txt > compare_MIE.txt
    tail -n 10 output/*_MIST.txt > compare_MIST.txt
    
    diff compare_MIE.txt test_system/sample_output/sample_output_MIE.txt
    diff compare_MIST.txt test_system/sample_output/sample_output_MIST.txt

  You should get no output at all here. If you do get output, first check if the files 
  "compare_MIE.txt" and "compare_MIST.txt" do exist in your top directory and 
  are non-empty. If this is not the case, check the file output/log.txt for further hints
  on what went wrong. If however you are using more exotic hardware and the files only show minor numerical differences 
  to "test_system/sample_output/sample_output_MIE.txt" and 
  "test_system/sample_output/sample_output_MIST.txt" respectively,
  this might well be due to machine precision and therefore be perfectly okay.
  Generally perform multiple of these tests and check that you are constantly getting 
  the same results. If this is not the case, that means that your system is not executing 
  the code in a proper manner (run condition). 
  
  With that said, I have been using the code now for ten years and never observed any deviations from
  previously calculated results. Neither on a XEON cluster, a Threadripper CPU nor even on
  a Raspberry Pi using ARM architecture.
  
  After executing the "run.sh" script, the compiled executables are located in the folder
  "exec". If this is not the case, you should check the file "output/log.txt" to see what went wrong.
  For compilation, the "Makefile" in the top directory is processed by "run.sh", so this is where you might 
  want to start troubleshooting (or maybe fine tuning).
    
### 2.2) Docker

  For building Docker images, two Dockerfiles have been added in the dockerfiles directory, one setting up 
  Ubuntu 18.04, the other setting up CentOS7. A symbolic link named "Dockerfile", pointing to the Ubuntu Dockerfile,
  has been provided in the main directory for convenience. You can switch to CentOS via 

    rm Dockerfile; ln -s dockerfiles/Dockerfile_centos7 ./Dockerfile

  if you prefer. Note that during testing, I had issues with MPI memory binding. On a single CPU, this should 
  be perfectly fine. On a cluster however, this can lead to problems, which is why such usage is discouraged. 
  While the provided Dockerfiles might work in a distributed computation environment, I did not perform any testing
  in this perspective.

  Building and running the Docker image works in the following manner. Note that you might need root privileges
  on the host machine for the following commands to be able to talk to the Docker daemon. Build the image via

    docker build -t parent

  PARENT is compiled during the build. A user named "docker" with sudo privileges is created in the image (if you want 
  to become root, issue "sudo su"). The password for the user is set to "docker", you might want to change this. 
  You can start a container from the image issuing

    docker run -it -v $(pwd):/PARENT --name parent parent

  This logs you in the container as the docker user. The "-v $(pwd):/PARENT" maps your current host directory to the
  /PARENT directory in the container. Changes in the host directory are reflected in the container /PARENT directory
  and vice versa. Note that this means that the PARENT binaries compiled during the image building step are overlayed
  by the host working directory. In particular, if your current mapped host working directory contains binaries
  compiled on the host (or for the other Dockerfile, Ubuntu vs. CentOS), problems will arise. In this case, from within 
  the container, issue 

    make clean; make

  If you do not want such a mapping to the host directory, remove the "-v $(pwd):/PARENT" switch. In this case, you
  will acccess the binaries compiled during the image building stage. Note, however, that you need to transfer your
  calculated files manually from the container to the host (refer to the according docker documentation in this case).
  If you delete the container without having the files transferred and without directory mapping, your calulations are 
  lost.

  After stopping the container using the "exit" command inside the container shell, you can restart it via

    docker container restart parent; docker container attach parent

  The Dockerfiles contain a short summary of the commands given here. Also, you might want to install
  additional tools during bulding the docker image. The according location to do so is marked with a comment.

  Note that binding more than one MPI process per CPU using Docker turned out problematic during my testing, so it is
  discouraged.  
  
  
## 2) RUNNING YOUR OWN TRAJECTORIES
  
  The easiest way to do this is by just modifying the file "parameters" in the top directory.
  
  The input files (GROMACS .top and .xtc) as well as the output folder should be specified either 
  relative to the top directory of the package or as an absolute path.
  
  The "NAME" parameter specifies the prefix of all output files and is set by default to the name 
  of the top folder containing the package. If you want to do multiple calculations using the same 
  output folder, you should change this parameter every time. Otherwise your previous results will 
  be overwritten without warning.
  
  The BINS parameters control how many bins are used for building the histograms which sample the 
  probability densities of the degrees of freedom(B=bonds, A=angles, D=dihedrals(torsions)).
  For the 2D values, a multiple of these values will be taken for building the histograms.
  Meaning if BBINS2D=40 and ABINS2D=60, 2400 bins will be taken for every 2D histogram
  relating a bond to an angle.
  These parameters are used (exclusively) during the Mutual Information Expansion (MIE) calculations, performed 
  by the program PARENT.x, which can be considered the core of this suite. Using 50 of them seems 
  to be a good starting point for ~10^6 frames of the MD trajectory. If you have considerably more
  (less) frames to process, you might want to tweak these parameters slightly up (down).
  In doubt, just leave these parameters as they are.
  
  The "BACKBONE_ATOMS" parameter is designed to find rigid torsion angles in your topology, 
  e. g. at a protein backbone. The names of the atoms here should match the names from 
  the GROMACS .top file. Other dihedrals are declared relative to the backbone dihedrals if they 
  share 3 atoms with them (and are then termed phaseangles). This improves the accuracy of the
  entropy calculations.
  
  In the case you want to run a considerable amount of trajectories so that compilation time 
  (only a few seconds) is an issue for you, you might want to uncomment the lines
  
    ./scripts/init.sh
    
  and
  
    make >> ${OUTDIR}/log.txt 2>&1
    
  in run.sh (after you compiled successfully for the first time), but keep in mind to recompile with
  	
    make clean; make
  
  when you change computer architecture or opering system (e. g. heterogeneous cluster, docker) and 
  also make sure that you are not overwriting previous results (change NAME in the parameterfile).
      
## 3) EXPLANATION OF THE PROGRAMS
  
  
### 3.1) BAT_builder.x
    
This program converts your GROMACS topology (.top) and trajectory (.xtc) files to binary .bat
trajectory files, which are then processed by PARENT.x, the core program of
this suite. The main purpose of this program is to convert every frame in the .xtc file,
which is stored in Cartesian coordinates, to internal bond-angle-torsion (BAT, Z-matrix) coordinates.
Furthermore additional information is attached to the header of the resulting .bat file, namely 

	-a version number
	
	-the precision of the file (single or double precision)
	
	-the number of non-redundant torsion angles(dihedrals) in the system 
	(which relates to the number of atoms by #atoms = #torsions + 3)
	
	-the number of frames in the trajectory
	
	-a list of all non-redundant torsion angles
	(specifying their constituent atoms using the atom numbers from the .top file)
	
	-atom weights of all atoms in the system (not used yet)
	
	-atom names, residue names, residue numbers and molecule names for every atom in the system



The program is used in the following manner:

	./BAT_builder.x -t input.top -x input.xtc -o output.bat -bb "BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ..." [--single_precision]

input.top, input.xtc and output.bat are self-explanatory.

"BackboneAtomName1 BackboneAtomName2 BackboneAtomName3 ..." lists the names of atoms belonging to a rigid backbone
as stated in the .top file, e. g.  "CA C N H1 O1" for a protein. Phaseangles are defined relative to a rigid dihedral. Also see section 2 for further information.

\[single_precision\] (the square brackets indicate optional), if set, writes the .bat trajectory in single precision instead of double precision,
which is discouraged, since all calculation is done in double precision anyway. Only use single precision if you are short of harddisk storage.




Additionally, the program can perform a back-conversion from .bat to .xtc, which is done by issuing the following command:

	./BAT_builder.x -b input.bat -o output.xtc


When the trajectory of a complex consisting of more than one molecule is converted, non-physical bonds (termed pseudo-bonds)
are added in order for the complex to have a connected topology. This is done automatically. Also the algorithm guarantees
that the chosen topology for every molecule in the complex is consistent with the topology which would be chosen if the molecules would be treated separately, in isolation. 


### 3.2) PARENT.x

This program can be considered the core of this suite. It calculates the configurational entropy according to the pairwise 
Mutual Information Expansion (MIE). For more detailed information, please read our article in the Journal of Chemical Theory and Computation  

"PARENT: A Parallel Software Suite for the Calculation of Configurational Entropy in Biomolecular Systems"  
DOI: 10.1021/acs.jctc.5b01217 

and additionally

B. J. Killian, J. Y. Kravitz, and M. K. Gilson, J. Chem. Phys. 127: 024107 (2007).  
J. Numata, and E.-W. Knapp, J. Chem. Theory Comput. 8: 1235 (2012).  

Calculations are done using MPI as well as openMP parallelization, making feasible the calculation of the 
configurational entropy using MIE for large, biologically realistic molecules, considerable sampling or large sets of structures.
Calculation time is scaling quadratically with the number of atoms. 

PARENT.x takes a .bat file as an input. Make sure that you provide (in total, from all compute nodes) reasonably more RAM 
than the size of the .bat file, otherwise Linux uses the hard-disk for temporary storage, which in the best case slows down the calculation considerably.
In the worst case scenario, the program might be killed due to memory shortage or crash without a warning.


The program is used in the following manner:

	./PARENT.x input.bat entropy.par \#bondsbins1D \#anglesbins1D \#dihedralsbins1D \#bondsbins2D \#anglesbins2D \#dihedralsbins2D

input.bat is the result from the conversion to internal BAT coordinates done with BAT_builder.x.

entropy.par is the binary output file, containing all the 1D and 2D entropy terms (from which the mutual information terms can be easiliy calculated).
The header of the file includes the same information as for the .bat file and additionally

\#bondsbins1D \#anglesbins1D \#dihedralsbins1D \#bondsbins2D \#anglesbins2D \#dihedralsbins2D
,the numbers of bins which were used for the entropy calculation. See section 2 for further information.


To have more control over the MPI process bindings, use e. g.

	mpirun --use-hwthread-cpus --bind-to hwthread --mca btl vader,self,tcp --report-bindings --map-by ppr:1:node:pe=16 ./PARENT.x input.bat entropy.par #bondsbins1D #anglesbins1D #dihedralsbins1D #bondsbins2D #anglesbins2D #dihedralsbins2D


### 3.3) get_PAR_MIE.x

The output of PARENT.x is a binary file (.par), so get_PAR_MIE.x is provided to output its contents to stdout in text form.
First it extracts a list of all 1D entropy terms of all degrees of freedom (bonds, angles, torsions). Then follows
a list of all 2D entropy terms between all pairs of degrees of freedom, with an additional last column, which
represents the pairwise mutual information. In the end summed values per category (bonds, angles, dihedrals)
are written, as well as the total entropy computed.


The program is used in the following manner:

	./get_PAR_MIE.x input.par

It can be started via MPI using an analogue command as for PARENT.x


### 3.4) get_PAR_info.x

The previous program (get_PAR_MIE.x) outputs all entropy/mutual information terms, but gives no information about
which atoms constitute a specific degree of freedom (e.g. bond 729). For that purpose get_PAR_info.x is provided. 

The output starts with the bonds indexing. The columns are indicating:
	
	#bond #atom1 #atom2 

In the same manner follow the angles:
	
	#angle #atom1 #atom2 #atom3

The dihedrals contain (additionally to #atom4) two special columns:

	#dihedral #atom1 #atom2 #atom3 #atom4    dihedral_type   phaseangle_of

The column "dihedral_type" is set to 0 for a common dihedral, to -1 for an improper dihedral, and to 1 if the dihedral contains a pseudo bond (see subsection 3.1. For every pseudo bond, there will be 3 pseudo dihedrals. The 2 atoms which all of them share
constitute the pseudo bond.).

The last column is set to -1 if the dihedral is not a phase angle (see subsection 3.1). Otherwise it will be set to the index of the dihedral to which it is relative to.

The program is used in the following manner:

	./get_PAR_info.x input.par


### 3.5) get_PAR_MIST.x

Although based on a different mathematical framework than MIE, the Maximum Information Spanning Tree (MIST) approximation relies on the same 
terms to be computed as for MIE. Empirically it seems to demonstrate superior convergence properties, so from a computational perspective
one is tempted to consider MIST a refinement of MIE. We highly recommend applying get_PAR_MIST.x to your output .par file of PARENT.x (at least if you are 
interested in total configurational entropy values). This program also makes use of MPI/openMP hybrid parallelization, which significantly reduces
computation time. In addition to our article in the Journal of Chemical Theory and Computation   

"PARENT: A Parallel Software Suite for the Calculation of Configurational Entropy in Biomolecular Systems"  
DOI: 10.1021/acs.jctc.5b01217 

please read

B. M. King, N. W. Silver, and B. Tidor. J. Phys. Chem. B 116: 2891 (2012).

The program is used in the following manner:

	./get_PAR_MIST.x input.par output.txt

Remember that you need to invoke it using MPI (see subsection 3.2 for further details).

## 4) IO library

The PARENT suite comprises a library for easily accessing its binary file formats. To use it, link with 

	obj/io.o
	
and use

	#include src/util/io/io.h
	
in your c++ source files. The Makefile in the top directory additionally hints the proper usage.



## 5) CONTACT INFORMATION

If you have any questions, feel free to contact: 

markus.fleck@univie.ac.at


    
    
    
    
    
    
    
    
    
    
    
    
    


      
      
      
    
    
    
  
  
  
  
  
  

