OBJS1 = bsg.o twoloophiggs.o suspect2.o
OBJS2 = sdecay.o lightst4bod.o Xvegas.o lightst3bod.o
OBJS3 = susylha.o hgaga.o hdecay.o

FC=gfortran

.f.o: 
	$(FC) -c -finit-local-zero -fbacktrace -fno-align-commons $*.f

susyhit:$(OBJS1) $(OBJS2) $(OBJS3) 
	$(FC) $(OBJS1) $(OBJS2) $(OBJS3) -o run

clean: 
	rm -f *.o