A = 1
B = 2000
while [ $B -le 100000 ]
do
  echo $A
  cat simpact.wrapper.easyABC.infectiousness.run.R | sed -e  "s/min.chunk <- 1/min.chunk <- ${A}/" -e "s/max.chunk <- 2000/max.chunk <- ${B}/" -e "s/design.points = 2000/design.points = 100000/" > simpact.wrapper.easyABC.infectiousness.run.${A}.R

  cat gent.testABCreject.SimpactRun.pbs | sed -e "s/simpact.wrapper.easyABC.run.1.to.2000/simpact.wrapper.easyABC.2run.${A}.to.${B}/" -e "s/simpact.wrapper.easyABC.infectiousness.run.1.R/simpact.wrapper.easyABC.infectiousness.2run.${A}.R/" > gent.testABCreject.SimpactRun.${A}.pbs
  qsub gent.testABCreject.SimpactRun.${A}.pbs

  sleep 5

  A=`expr $B + 1`
  B=`expr $B + 2000`
done
