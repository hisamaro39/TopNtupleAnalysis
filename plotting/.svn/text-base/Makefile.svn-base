
ROOTINC = $(shell root-config --cflags)
ROOTLIB = $(shell root-config --libs)

#CXXFLAGS = -c -g -std=c++11 $(ROOTINC)
#LDFLAGS = -g -std=c++11 -rdynamic $(ROOTLIB)

CXXFLAGS = -c -O4 -std=c++11 $(ROOTINC)
LDFLAGS = -O4 -std=c++11 -rdynamic $(ROOTLIB)

DEPS = ParseUtils.o utils.o Hist.o SystematicImplementation.o SampleSet.o SystematicCalculation.o

all: plot compare plotCompareNominal plotChannelRatio

%.o:%.cxx
	g++ $(CXXFLAGS) -o $@ $<


plot: $(DEPS) plot.o
	g++ -o $@ $(DEPS) plot.o $(LDFLAGS)

plotChannelRatio: $(DEPS) plotChannelRatio.o
	g++ -o $@ $(DEPS) plotChannelRatio.o $(LDFLAGS)

plotCompareNominal: $(DEPS) plotCompareNominal.o
	g++ -o $@ $(DEPS) plotCompareNominal.o $(LDFLAGS)

compare: $(DEPS) compare.o
	g++ -o $@ $(DEPS) compare.o $(LDFLAGS)

clean:
	rm -f *.o plot compare plotCompareNominal

