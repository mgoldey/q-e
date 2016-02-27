include ../Makefile.inc
all: left.out right.out coupling.in coupling.out left_pp.out right_pp.out lefts.out rights.out couplings.in couplings.out left_pps.out right_pps.out

clean: 
	rm -rf *.out out coupling.in

left.out: left.in
	$(PW_CMD) < $^ > $@

right.out: right.in
	$(PW_CMD) < $^ > $@

coupling.in: left.out right.out
	sh ../../scripts/setup_coupling_input.sh left.out right.out > coupling.in

coupling.out: coupling.in
	$(CUP_CMD) < $^ > $@

left_pp.out: left_pp.in
	$(PP_CMD)  < $^ > $@

right_pp.out: right_pp.in
	$(PP_CMD)  < $^ > $@

lefts.out: lefts.in
	$(PW_CMD) < $^ > $@

rights.out: rights.in
	$(PW_CMD) < $^ > $@

couplings.in: lefts.out rights.out
	sh ../../scripts/setup_coupling_input.sh lefts.out rights.out > couplings.in

couplings.out: couplings.in
	$(CUP_CMD) < $^ > $@

left_pps.out: left_pps.in
	$(PP_CMD)  < $^ > $@

right_pps.out: right_pps.in
	$(PP_CMD)  < $^ > $@
