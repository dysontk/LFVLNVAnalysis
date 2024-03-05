CXX = g++
FC = g77

# ROOTINC= -I/home/gangli/Documents/Simulation/MG5_aMC_v2_6_6/ExRootAnalysis/ExRootAnalysis
# ROOTLIB= -L/home/gangli/Documents/Simulation/MG5_aMC_v2_6_6/ExRootAnalysis -lExRootAnalysis

ROOTINC= -I/home/dysontravis/MG5_aMC_v3_4_0/ExRootAnalysis/ExRootAnalysis
ROOTLIB= -L/home/dysontravis/MG5_aMC_v3_4_0/ExRootAnalysis -lExRootAnalysis

# ROOTFLAG=${shell root-config --cflags --libs}
ROOTFLAG=${shell root-config --cflags --glibs} -lEG

FJFLAG=${shell /home/dysontravis/MG5_aMC_v3_4_0/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins}

# ROOTINC2= -I/home/gangli/Documents/Simulation/MG5_aMC_v2_6_6/Delphes/external \
#          -I/home/gangli/Documents/Simulation/MG5_aMC_v2_6_6/Delphes/display \
#          -I/home/gangli/Documents/Simulation/MG5_aMC_v2_6_6/Delphes
# ROOTLIB2= -L/home/gangli/Documents/Simulation/MG5_aMC_v2_6_6/Delphes -lDelphes

ROOTINC2= -I/home/dysontravis/MG5_aMC_v3_4_0/Delphes/external \
         -I/home/dysontravis/MG5_aMC_v3_4_0/Delphes/display \
         -I/home/dysontravis/MG5_aMC_v3_4_0/Delphes
ROOTLIB2= -L/home/dysontravis/MG5_aMC_v3_4_0/Delphes -lDelphes



CPPFLAGS=${ROOTFLAG} ${ROOTINC} ${FJFLAG} ${ROOTINC2}  #${ROOTINC2}  #${FJFLAG} ${ROOTINC}
LDFLAGS= ${ROOTFLAG} ${ROOTLIB} ${FJFLAG} ${ROOTLIB2} #${ROOTLIB} #	${ROOTLIB2} #${FJFLAG} ${ROOTLIB}
#FOBJS := $(patsubst %.f,%.o,$(wildcard *.f)) 
CCOBJS := $(patsubst %.cpp,%.o,$(wildcard *.cpp))

main:  ${CCOBJS}
	$(CXX) *.o  ${LDFLAGS} -o  $@

%.o : %.cpp
	$(CXX) -c $(CPPFLAGS) $& $<

clean:
	rm *.o main
