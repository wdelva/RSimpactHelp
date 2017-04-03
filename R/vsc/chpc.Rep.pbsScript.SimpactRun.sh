A = 6001
B = 8000
while [ $B -le 22000 ]
do
  echo $A
  cat simpact.wrapper.easyABC.run.4001.R | sed -e "s/min.chunk <- 4001/min.chunk <- ${A}/" -e "s/max.chunk <- 6000/max.chunk <- ${B}/" > simpact.wrapper.easyABC.run.${A}.R

  cat chpcSimpactRunTest.pbs | sed -e "s/simpact.wrapper.easyABC.run.4001/simpact.wrapper.easyABC.run.${A}/" -e "s/simpact.wrapper.easyABC.run.4001.R/simpact.wrapper.easyABC.run.${A}.R/" > chpcSimpactRunTest.${A}.pbs

  qsub chpcSimpactRunTest.${A}.pbs

  sleep 4

  A = `expr $B + 1`
  B = `expr $B + 2000`
done
