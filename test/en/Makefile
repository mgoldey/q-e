include ../Makefile.inc
all: zero.out rerun.out
clean: 
	rm -rf out *.out tmp.pp *.cube

# the meat
zero.out: zero.in
	$(PW_CMD) < $^ > $@
rerun.out: rerun.in
	$(PW_CMD) < $^ > $@
