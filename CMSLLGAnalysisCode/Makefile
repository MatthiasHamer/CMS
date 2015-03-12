ROOTGLIBS =`root-config --glibs` -lMinuit
LIBS = `root-config --glibs` -L/usr/X11R6/lib -lm -ldl -lstdc++ -Wl,--no-as-needed $(ROOTGLIBS)

CXXFLAGS = -O3 -Wall --exceptions `root-config --ldflags --cflags` -I./

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

RunAnalysis: LLGAnalysis.o RunAnalysis.o
	$(CXX) -o RunAnalysis RunAnalysis.o LLGAnalysis.o $(LDFLAGS) $(LIBS)
clean:
	rm -f RunAnalysis
	rm -f *.o
	rm -f *~
	rm -f *.so
	rm -f Loader_C.d
