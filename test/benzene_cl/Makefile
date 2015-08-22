include ../Makefile.inc
all: left.out  right.out coupling.out
clean: 
	rm -rf left.out right.out out coupling.out

left.out: left.in
	$(PW_CMD) < $^ > $@

right.out: right.in
	$(PW_CMD) < $^ > $@

coupling.out: coupling.in
        sh ../../scripts/setup_coupling_input.sh left.out right.out > $^
	$(CUP_CMD) < $^ > $@
