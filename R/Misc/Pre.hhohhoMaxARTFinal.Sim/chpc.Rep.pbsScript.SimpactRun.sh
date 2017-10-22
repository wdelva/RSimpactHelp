A=1
B=8000
while [ $B -le 10000 ]
do
  echo $A
  cat pre.hho.simpact.wrapper.easyABC.run.R | sed -e "s/min.chunk <- 4001/min.chunk <- ${A}/" -e "s/max.chunk <- 6000/max.chunk <- ${B}/" > simpact.wrapper.easyABC.run.${A}.R

  cat pre.hho.chpcSimpactRun.pbs | sed -e "s/simpact.wrapper.easyABC.run.4001/simpact.wrapper.easyABC.run.${A}/" -e "s/simpact.wrapper.easyABC.run.4001.R/simpact.wrapper.easyABC.run.${A}.R/" > chpcSimpactRunTest.${A}.pbs

  qsub chpcSimpactRunTest.${A}.pbs

  sleep 4

  A = `expr $B + 1`
  B = `expr $B + 500`
done
