#!/bin/bash

#source $1 

TTT=`grep "(Failed)" /scratch/albany/nightly_log_clangdbg.txt -c`
TTTT=`grep "(Not Run)" /scratch/albany/nightly_log_clangdbg.txt -c`
TTTTT=`grep "(Timeout)" /scratch/albany/nightly_log_clangdbg.txt -c`
TT=`grep "...   Passed" /scratch/albany/nightly_log_clangdbg.txt -c`

#/bin/mail -s "Albany ($ALBANY_BRANCH): $TTT" "albany-regression@software.sandia.gov" < $ALBOUTDIR/albany_runtests.out
/bin/mail -s "Albany64BitClangDbg, cee-compute011: $TT tests passed, $TTT tests failed, $TTTT tests not run, $TTTTT timeouts" "ikalash@sandia.gov, daibane@sandia.gov, tjfulle@sandia.gov, mperego@sandia.gov, lbertag@sandia.gov" < /scratch/albany/results_clangdbg
#/bin/mail -s "Albany64BitClang, cee-compute016: $TT tests passed, $TTT tests failed, $TTTT tests not run, $TTTTT timeouts" "ikalash@sandia.gov" < /scratch/albany/results_clangdbg
