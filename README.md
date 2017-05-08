# External Potential Constrained Density Functional Theory (EPCDFT)

This is CDFT implemented in QE.

## Installation

Clone the repository and compile pw and pp following the normal compilation guidelines for QE. 

## Documentation

A tutorial can be found [here](http://www.nicholasbrawand.com/constrained-density-functional-theory). Documentation for CDFT in PW can be found in PW/Doc/INPUT_PW.def and the coupling code documentation is in PP/Doc/INPUT_CDFT.def

NOTA BENE:  Until the new XML format is stable, add -D__OLDXML to MANUAL_DFLAGS in your make.inc

## Tests

Nightly tests of the master branch are run on UofC RCC Midway. 

## Contacts

Nick Brawand nicholasbrawand@gmail.com or Matthew Goldey matthew.goldey@gmail.com

## License

Please do not distribute.
