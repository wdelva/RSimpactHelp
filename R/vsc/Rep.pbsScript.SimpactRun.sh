A = 1
B = 400
while [ $B -le 8000 ]
do
  echo $A
  cat simpact.wrapper.easyABC.run.R | sed  -e  "s/main.filename <- .*/main.filename <- 'INPUT.df-10Points10Par2016-11-02.csv'/" -e  "s/min.chunk <- .*/min.chunk <- ${A}/" -e "s/max.chunk <- .*/max.chunk <- ${B},/" -e "s/sim_repeat <- .*/sim_repeat <- 16/" -e "s/ncluster <- .*/ncluster <- 5/" > simpact.wrapper.easyABC.run.${A}.R
  
  cat testABCrejectSimpactRun.pbs | sed -e "s/simpact.wrapper.easyABC.run.8000/simpact.wrapper.easyABC.run.${A}/" -e "s/walltime=00:30:00/walltime=04:00:00/" -e "s/simpact.wrapper.easyABC.run.R/simpact.wrapper.easyABC.run.${A}.R/" > testABCrejectSimpactRun.${A}.pbs
  qsub testABCreject.AIMS.${A}.pbs
  sleep 2
  A =`expr $B + 1`
  B =`expr $B + 400`
done

cd /user/data/gent/vsc400/vsc40070/simpact-test/code/