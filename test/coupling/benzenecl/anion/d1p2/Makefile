include ../../../../Makefile.inc
all: benzene.out cl.out coupling.out
clean: 
	rm -rf *.out out 

benzene.out: benzene.in
	$(PW_CMD) < $^ > $@

cl.out: cl.in
	$(PW_CMD) < $^ > $@

coupling.out: coupling.in
	$(CUP_CMD) < $^ > $@
