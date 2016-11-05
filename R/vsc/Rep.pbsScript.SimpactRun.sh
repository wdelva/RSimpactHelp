A=1
B=10
while [ $B -le 100 ]
do
  echo $A
  cat simpact.wrapper.easyABC.run.R | sed  -e  "s/main.filename <- .*/main.filename <- 'INPUT.df-100Points10Par2016-11-03.csv'/" -e  "s/min.chunk <- .*/min.chunk <- ${A}/" -e "s/max.chunk <- .*/max.chunk <- ${B}/" -e "s/sim_repeat <- .*/sim_repeat <- 2/" -e "s/ncluster <- .*/ncluster <- 16/" > simpact.wrapper.easyABC.run.${A}.R

  cat testABCreject.SimpactRun.pbs | sed -e "s/simpact.wrapper.easyABC.run.100/simpact.wrapper.easyABC.run.${A}/" -e "s/walltime=00:30:00/walltime=04:00:00/" -e "s/simpact.wrapper.easyABC.run.R/simpact.wrapper.easyABC.run.${A}.R/" > testABCreject.SimpactRun.${A}.pbs
  qsub testABCreject.SimpactRun.${A}.pbs
  sleep 2
  A=`expr $B + 1`
  B=`expr $B + 10`
done
