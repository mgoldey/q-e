include ../Makefile.inc
all: left.out 
clean: 
	rm -rf left.out # right.out out

left.out: left.in
	$(PW_CMD) < $^ > $@

#right.out: right.in
#	$(PW_CMD) < $^ > $@
#
#coupling.out: coupling.in
#	$(CUP_CMD) < $^ > $@
