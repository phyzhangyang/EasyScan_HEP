all: susyhit gm2calc feynhiggs higgsbounds

# SUSYHIT
SUSYHIT = ./susyhit
susyhit:
	[ -d $(SUSYHIT) ] || (wget https://www.itp.kit.edu/~maggie/SUSY-HIT/susyhit.tar.gz && mkdir $(SUSYHIT) && tar -xf susyhit.tar.gz -C $(SUSYHIT) && rm susyhit.tar.gz)	
	[ -f $(SUSYHIT)/run ] || (cd $(SUSYHIT) && make)
	(cp -r for_MSSM_mW.susyhit $(SUSYHIT)/suspect2_lha.in)


# FeynHiggs
FeynHiggs = ./FeynHiggs
feynhiggs:
	[ -d $(FeynHiggs) ] || (wget --no-check-certificate https://wwwth.mpp.mpg.de/members/heinemey/feynhiggs/newversion/FeynHiggs-2.19.0.tar.gz && tar -xf FeynHiggs-2.19.0.tar.gz && mv FeynHiggs-2.19.0 $(FeynHiggs) && rm FeynHiggs-2.19.0.tar.gz)	
	[ -f $(FeynHiggs)/"$(shell uname -m)-$(shell uname -s)"/lib64/libFH.a ] || (cd $(FeynHiggs) && ./configure && make && make install)

# HiggsBounds
HiggsBounds = ./HiggsBounds
higgsbounds:
	[ -d $(HiggsBounds) ] || (wget https://gitlab.com/higgsbounds/higgsbounds/-/archive/5.10.2/higgsbounds-5.10.2.tar.gz && tar -xf higgsbounds-5.10.2.tar.gz && mv higgsbounds $(HiggsBounds) && rm higgsbounds-5.10.2.tar.gz)	
	[ -f $(HiggsBounds)/build/example_programs/HBSLHAinputblocksfromFH ] || (cd $(HiggsBounds) && mkdir build && cd build && cmake -DFeynHiggs_ROOT=../../FeynHiggs/ .. && make)


# GM2Calc
GM2Calc = ./GM2Calc
gm2calc:
	[ -d $(GM2Calc) ] || (wget https://github.com/GM2Calc/GM2Calc/archive/v2.1.0.tar.gz && tar -xf v2.1.0.tar.gz && mv GM2Calc-2.1.0 $(GM2Calc) && rm v2.1.0.tar.gz)	
	[ -f $(GM2Calc)/bin/gm2calc.x ] || (cd $(GM2Calc) && cmake . && make)




.PHONY: clean
clean: 
	rm -rf $(SUSYHIT) $(FeynHiggs) $(HiggsBounds) $(GM2Calc)
