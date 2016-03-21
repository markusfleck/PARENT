all : exec/BAT_builder.x exec/PARENT.x exec/get_PAR_MIST.x par_mie

clean :
	rm -r exec obj




COPS = -Wall


obj :
	mkdir obj

exec :
	mkdir exec
 
	
obj/io_binary.o : src/util/io/io_binary.cpp src/util/io/io_binary.h | obj
	g++ -c src/util/io/io_binary.cpp -o obj/io_binary.o $(COPS)

obj/io_text.o : src/util/io/io_text.cpp src/util/io/io_text.h | obj
	g++ -c src/util/io/io_text.cpp -o obj/io_text.o $(COPS)

obj/io.o : obj/io_binary.o obj/io_text.o | obj
	ld -r obj/io_binary.o obj/io_text.o -o obj/io.o

# the utility library of the PARENT suite:
obj/util.o : src/util/util.cpp src/util/util.h | obj
	g++ -c src/util/util.cpp -o obj/util.o $(COPS)
	

# this program converts the GROMACS .xtc and .top files to a binary .bat file in bond-angle-torsion coordinates
exec/BAT_builder.x : src/util/io/xdrfile/xdrfile.c src/util/io/xdrfile/xdrfile_xtc.c src/util/io/xdrfile/xdrfile.h src/util/io/xdrfile/xdrfile_xtc.h src/BAT_builder/bat.cpp src/BAT_builder/bat.h src/BAT_builder/topology.cpp src/BAT_builder/BAT_topology.cpp src/BAT_builder/BAT_trajectory.cpp obj/io.o src/util/io/io.h src/util/io/io_binary.h src/util/io/io_text.h | obj exec
	gcc -c src/util/io/xdrfile/xdrfile.c -o obj/xdrfile.o
	gcc -c src/util/io/xdrfile/xdrfile_xtc.c -o obj/xdrfile_xtc.o
	g++ src/BAT_builder/bat.cpp obj/xdrfile.o  obj/xdrfile_xtc.o obj/io.o -o exec/BAT_builder.x $(COPS)


# the "main" program: calculates the entropy terms from the .bat file and puts out a binary .par file
exec/PARENT.x : src/PARENT/PARENT.cpp obj/io.o src/util/io/io.h src/util/io/io_binary.h src/util/io/io_text.h | exec
	mpiCC -fopenmp src/PARENT/PARENT.cpp obj/io.o -o exec/PARENT.x $(COPS)

# this program calculates the MIST approximation from the binary .par file and outputs it in text form
exec/get_PAR_MIST.x : src/process_output/get_PAR_MIST.cpp  obj/io.o src/util/io/io.h src/util/io/io_binary.h src/util/io/io_text.h obj/util.o src/util/util.h | exec
	mpiCC -fopenmp src/process_output/get_PAR_MIST.cpp  obj/io.o obj/util.o -o exec/get_PAR_MIST.x $(COPS)


par_mie : exec/get_PAR_info.x exec/get_PAR_MIE.x


# this program reads the binary .par file and extracts the bond-angle-torsion topology from it (in text format to stdout)
exec/get_PAR_info.x : src/process_output/get_PAR_info.cpp obj/io.o src/util/io/io.h src/util/io/io_binary.h src/util/io/io_text.h | exec
	g++ src/process_output/get_PAR_info.cpp obj/io.o -o exec/get_PAR_info.x $(COPS)

# this program reads the binary .par file and extracts the entropy and mutual information terms from it (in text format to stdout)
exec/get_PAR_MIE.x : src/process_output/get_PAR_MIE.cpp obj/io.o src/util/io/io.h src/util/io/io_binary.h src/util/io/io_text.h obj/util.o src/util/util.h | exec
	g++ src/process_output/get_PAR_MIE.cpp obj/io.o obj/util.o -o exec/get_PAR_MIE.x $(COPS)




