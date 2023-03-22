all: susyhit gm2calc feynhiggs

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
	[ -f $(FeynHiggs)/bin/gm2calc.x ] || (cd $(FeynHiggs) && cmake . && make)

# HiggsBounds
HiggsBounds = ./HiggsBounds
higgsbounds:
	[ -d $(HiggsBounds) ] || (wget https://higgsbounds.hepforge.org/downloads?f=HiggsBounds-5.3.2beta.tar.gz && tar -xf downloads?f=HiggsBounds-5.3.2beta.tar.gz && mv HiggsBounds-5.3.2beta $(HiggsBounds) && rm downloads?f=HiggsBounds-5.3.2beta.tar.gz)	
	[ -f $(HiggsBounds)/bin/gm2calc.x ] || (cd $(HiggsBounds) && cmake . && make)


# GM2Calc
GM2Calc = ./GM2Calc
gm2calc:
	[ -d $(GM2Calc) ] || (wget https://github.com/GM2Calc/GM2Calc/archive/v2.1.0.tar.gz && tar -xf v2.1.0.tar.gz && mv GM2Calc-2.1.0 $(GM2Calc) && rm v2.1.0.tar.gz)	
	[ -f $(GM2Calc)/bin/gm2calc.x ] || (cd $(GM2Calc) && cmake . && make)




.PHONY: clean
clean: 
	rm -rf $(SUSYHIT) $(GM2Calc)
