testdirs:= zn silane zero
all:
	for d in $(testdirs); do (cd $$d; $(MAKE)); done
clean:
	for d in $(testdirs); do (cd $$d; $(MAKE) clean); done

