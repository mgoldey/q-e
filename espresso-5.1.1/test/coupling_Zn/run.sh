#!/bin/bash

echo "running pw.in"
../../bin/pw.x < pw.in > pw.out
echo "done"
echo ""
echo "running pw2.in ( it is a copy of pw.in but with outdir = out2 )"
../../bin/pw.x < pw2.in > pw2.out
echo "done"
echo ""
echo "running coupling.in between pw.in and pw2.in"
../../bin/epcdft_coupling.x < coupling.in > coupling.out
echo "done"
