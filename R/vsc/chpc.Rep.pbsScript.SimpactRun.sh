A = 1
B = 50
while [ $B -le 1000 ]
do
  echo $A
  cat simpact.wrapper.easyABC.run.R | sed -e "s/comp <- .*/comp <- 'chpc'/" -e  "s/min.chunk <- .*/min.chunk <- ${A}/" -e "s/max.chunk <- .*/max.chunk <- ${B}/" -e "s/ncluster.use <- .*/ncluster.use <- 16/" > simpact.wrapper.easyABC.run.${A}.R

  cat chpc.testABCreject.SimpactRun.pbs | sed -e "s/simpact.wrapper.easyABC.run.100/simpact.wrapper.easyABC.run.${A}/" -e "s/walltime=00:30:00/walltime=04:00:00/" -e "s/simpact.wrapper.easyABC.run.R/simpact.wrapper.easyABC.run.${A}.R/" > testABCreject.SimpactRun.${A}.pbs

  qsub testABCreject.SimpactRun.${A}.pbs

  sleep 4

  A = `expr $B + 1`
  B = `expr $B + 50`
done
