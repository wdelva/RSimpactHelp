A=1
B=2
while [ $B -le 2 ]
do
  echo $A
  qsub testABCreject.SimpactRunEMU.pbs
  sleep 5
  cat testABCreject.SimpactRunEMU.pbs | sed -e "s/simpact-emulator-knit-pca-CLUSTER/simpact-emulator-knit-CLUSTER/" -e "s/simpact-emulator-knit-pca-CLUSTER.R/simpact-emulator-knit-CLUSTER.R/" > testABCreject.SimpactRun.NONEPC.pbs
  qsub testABCreject.SimpactRun.NONEPC.pbs
  sleep 2
  A=`expr $B + 1`
  B=`expr $B + 10`
done
