#!/bin/sh


rm /tmp/invariants.txt

echo "compiling and installing ELINA" | colorize red
cd /home/cav2018/ELINA/
make clean
make 
make install

cd /home/cav2018/CAV2018/Benchmarks/ 
for file in *.c
do 
  basefile=`basename ${file}`
  echo "Running ${basefile} with ELINA" | colorize red
  timeout 3600 /home/cav2018/seahorn/build/run/bin/sea pf --crab --crab-dom=opt-pk-elina $file > /home/cav2018/CAV2018/results/ELINA/${basefile}_runtime.txt
  mv /tmp/invariants.txt /home/cav2018/CAV2018/results/ELINA/${basefile}_invariants.txt
done


echo "compiling and installing Poly-RL" | colorize blue
cd /home/cav2018/ELINA/poly_rl
make clean
make 
make install
cd /home/cav2018/CAV2018/Benchmarks/ 
for file in *.c
do
  basefile=`basename ${file}`
  echo "Running ${basefile} with Poly-RL" | colorize blue
  timeout 3600 /home/cav2018/seahorn/build/run/bin/sea pf --crab --crab-dom=opt-pk-elina $file > /home/cav2018/CAV2018/results/Poly-RL/${basefile}_runtime.txt
  mv /tmp/invariants.txt /home/cav2018/CAV2018/results/Poly-RL/${basefile}_invariants.txt
done

# echo "compiling and installing Poly-Fixed" | colorize green
# cd /home/cav2018/ELINA/poly_fixed
# make clean
# make 
# make install
# cd /home/cav2018/CAV2018/Benchmarks/ 
# for file in *.c
# do
#   basefile=`basename ${file}`
#   echo "Running ${basefile} with Poly-Fixed" | colorize green
#   timeout 3600 /home/cav2018/seahorn/build/run/bin/sea pf --crab --crab-dom=opt-pk-elina $file > /home/cav2018/CAV2018/results/Poly-Fixed/${basefile}_runtime.txt
#   mv /tmp/invariants.txt /home/cav2018/CAV2018/results/Poly-Fixed/${basefile}_invariants.txt
# done

# echo "compiling and installing Poly-Init" | colorize yellow
# cd /home/cav2018/ELINA/poly_init
# make clean
# make 
# make install
# cd /home/cav2018/CAV2018/Benchmarks/ 
# for file in *.c
# do
#   basefile=`basename ${file}`
#   echo "Running ${basefile} with Poly-Init" | colorize yellow
#   timeout 3600 /home/cav2018/seahorn/build/run/bin/sea pf --crab --crab-dom=opt-pk-elina $file > /home/cav2018/CAV2018/results/Poly-Init/${basefile}_runtime.txt
#   mv /tmp/invariants.txt /home/cav2018/CAV2018/results/Poly-Init/${basefile}_invariants.txt
# done

echo "compiling and installing Poly-DQN" | colorize blue
cd ~/Desktop/elina_dqn/test_extended_nn
make clean
make 
make install
cd /home/cav2018/CAV2018/Benchmarks/ 
for file in *.c
do
  basefile=`basename ${file}`
  echo "Running ${basefile} with Poly-DQN" | colorize blue
  timeout 3600 /home/cav2018/seahorn/build/run/bin/sea pf --crab --crab-dom=opt-pk-elina $file > /home/cav2018/CAV2018/results/Poly-RL/${basefile}_runtime.txt
  mv /tmp/invariants.txt /home/cav2018/CAV2018/results/Poly-RL/${basefile}_invariants.txt
done

cd /home/cav2018/ELINA/elina_poly
make 
make install


cd /home/cav2018/CAV2018/
gcc -o check_precision check_precision.c -loptpoly -lelinaux -lelinalinearize -lgmp -lmpfr -lm


echo "Comparing  runtime" 
cd /home/cav2018/CAV2018/Benchmarks/ 
for file in *.c
do
  basefile=`basename ${file}`  
  echo "Runtime of ELINA for ${basefile}$" | colorize red
  grep -i --color=always "optpoly" /home/cav2018/CAV2018/results/ELINA/${basefile}_runtime.txt
  echo "Runtime of Poly-RL for ${basefile}$" | colorize blue
  grep -i --color=always "optpoly" /home/cav2018/CAV2018/results/Poly-RL/${basefile}_runtime.txt
  echo "Runtime of Poly-Fixed for ${basefile}$" | colorize green
  grep -i --color=always "optpoly" /home/cav2018/CAV2018/results/Poly-Fixed/${basefile}_runtime.txt
  echo "Runtime of Poly-Init for ${basefile}$" | colorize yellow
  grep -i --color=always "optpoly" /home/cav2018/CAV2018/results/Poly-Init/${basefile}_runtime.txt
done


echo "Comparing  precision"
for file in *.c
do
  basefile=`basename ${file}`  
  echo "precision of Poly-RL for ${basefile}$" | colorize blue
  ../check_precision /home/cav2018/CAV2018/results/Poly-RL/${basefile}_invariants.txt /home/cav2018/CAV2018/results/ELINA/${basefile}_invariants.txt
  echo "precision of Poly-Fixed for ${basefile}$" | colorize green
  ../check_precision /home/cav2018/CAV2018/results/Poly-Fixed/${basefile}_invariants.txt /home/cav2018/CAV2018/results/ELINA/${basefile}_invariants.txt
  echo "precision of Poly-Init for ${basefile}$" | colorize yellow
  ../check_precision /home/cav2018/CAV2018/results/Poly-Init/${basefile}_invariants.txt /home/cav2018/CAV2018/results/ELINA/${basefile}_invariants.txt
done

