# Makefile created by myfmkmf.pl at 09:46:00 on Jan 5, 2012
#
# Modification of fmkmf.pl by Rick D. Saylor, NOAA/ARL/ATDD, Nov 2011
#
FC=gfortran 
CFLAGS=-O3 
LFLAGS=  

# create executable
access: GlobalData.o EnvironData.o PhysChemData.o Output.o Utils.o Initialize.o Emissions.o CanopyPhysics.o DryDep.o GasChem.o Chem.o DVODE_f90_m.o Integrate.o VertTransport.o Main.o 
	 $(FC) -o ../access $(LFLAGS) GlobalData.o EnvironData.o PhysChemData.o Output.o Utils.o Initialize.o Emissions.o CanopyPhysics.o DryDep.o GasChem.o Chem.o DVODE_f90_m.o Integrate.o VertTransport.o Main.o

# compilation dependencies

GlobalData.o:./GlobalData.f90  
	 $(FC) $(CFLAGS) -c ./GlobalData.f90 

PhysChemData.o:./PhysChemData.f90 GlobalData.o
	 $(FC) $(CFLAGS) -c ./PhysChemData.f90 

Utils.o:./Utils.f90 GlobalData.o PhysChemData.o
	 $(FC) $(CFLAGS) -c ./Utils.f90 

Output.o:./Output.f90 GlobalData.o Utils.o
	 $(FC) $(CFLAGS) -c ./Output.f90 

VertTransport.o:./VertTransport.f90 GlobalData.o Utils.o
	 $(FC) $(CFLAGS) -c ./VertTransport.f90 

CanopyPhysics.o:./CanopyPhysics.f90 GlobalData.o PhysChemData.o VertTransport.o Utils.o
	 $(FC) $(CFLAGS) -c ./CanopyPhysics.f90 

DryDep.o:./DryDep.f90 GlobalData.o CanopyPhysics.o PhysChemData.o Utils.o
	 $(FC) $(CFLAGS) -c ./DryDep.f90 

EnvironData.o:./EnvironData.f90 GlobalData.o CanopyPhysics.o PhysChemData.o Utils.o
	 $(FC) $(CFLAGS) -c ./EnvironData.f90 

Initialize.o:./Initialize.f90 GlobalData.o EnvironData.o PhysChemData.o DryDep.o Utils.o Output.o 
	 $(FC) $(CFLAGS) -c ./Initialize.f90 

Emissions.o:./Emissions.f90 GlobalData.o Utils.o
	 $(FC) $(CFLAGS) -c ./Emissions.f90 

GasChem.o:./GasChem.f90 GlobalData.o Utils.o
	 $(FC) $(CFLAGS) -c ./GasChem.f90 

DVODE_f90_m.o:./DVODE_f90_m.f90
	 $(FC) $(CFLAGS) -c ./DVODE_f90_m.f90

Chem.o:./Chem.f90 GlobalData.o GasChem.o DVODE_f90_m.o
	 $(FC) $(CFLAGS) -c ./Chem.f90 

Integrate.o:./Integrate.f90 GlobalData.o EnvironData.o Emissions.o DryDep.o Chem.o GasChem.o VertTransport.o Utils.o
	 $(FC) $(CFLAGS) -c ./Integrate.f90 

Main.o:Main.f90 GlobalData.o Initialize.o Integrate.o Output.o Utils.o 
	 $(FC) $(CFLAGS) -c Main.f90 

clean:
	rm -f -r *.o *.mod 

veryclean:
	rm -f -r *.o *.mod ../access 
