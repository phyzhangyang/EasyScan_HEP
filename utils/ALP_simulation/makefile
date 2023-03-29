all: mg5 alp_model ma5 sig 

# madgraph/madevent 5 with pythia 8 and delphes 3
MG5 = ./MG5
mg5:
	[ -d $(MG5) ] || (wget  && tar -xf MG5_aMC_v3.4.2.tar.gz && mv MG5_aMC_v3_4_2 $(SUSYHIT) && rm MG5_aMC_v3.4.2.tar.gz)	
	( cd $(MG5) && echo "install pythia8" | ./bin/mg5_aMC )
	( cd $(MG5) && echo "install Delphes" | ./bin/mg5_aMC )

# alp 
ALP = ./ALP
alp:
	[ -d $(ALP) ] || (wget http://feynrules.irmp.ucl.ac.be/raw-attachment/wiki/ALPsEFT/ALP_linear_UFO.tar.gz && tar -xf ALP_linear_UFO.tar.gz && mv ALP_linear_UFO $(ALP) && rm ALP_linear_UFO.tar.gz)	
	(cp -r $(ALP) $(MG5)/models)

# ma5 
MA5 = ./MA5
ma5:
	[ -d $(MA5) ] || (wget MadAnalysis5_v1.9.60.tgz && tar -xf MadAnalysis.tar.gz && mv madanalysis5 $(MA5) && rm MadAnalysis_v1.9.60.tgz)	
	[ -f $(HiggsBounds)/build/example_programs/HBSLHAinputblocksfromFH ] || (cd $(HiggsBounds) && mkdir build && cd build && cmake -DFeynHiggs_ROOT=../../FeynHiggs/ .. && make)


# GM2Calc
GM2Calc = ./GM2Calc
gm2calc:
	[ -d $(GM2Calc) ] || (wget https://github.com/GM2Calc/GM2Calc/archive/v2.1.0.tar.gz && tar -xf v2.1.0.tar.gz && mv GM2Calc-2.1.0 $(GM2Calc) && rm v2.1.0.tar.gz)	
	[ -f $(GM2Calc)/bin/gm2calc.x ] || (cd $(GM2Calc) && cmake . && make)




.PHONY: clean
clean: 
	rm -rf $(SUSYHIT) $(FeynHiggs) $(HiggsBounds) $(GM2Calc)
