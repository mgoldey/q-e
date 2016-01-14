include ../Makefile.inc
all: left.out right.out coupling.out
clean:
	rm -rf out *.out tmp.pp *.cube

# the meat
left.out: left.in
	$(PW_CMD) < $^ > $@

right.out: right.in
	$(PW_CMD) < $^ > $@

coupling.out: coupling.in
	../../scripts/setup_coupling_input.sh left.out right.out > coupling.in
	$(CUP_CMD) < $^ > $@
