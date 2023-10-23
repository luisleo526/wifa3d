NUM_OF_PROCESS=$1
echo $NUM_OF_PROCESS > /app/run/wt1.pro
cd /app
bash /app/scripts/preprocess.sh
bash /app/scripts/compile.sh
cd /app/run; mpirun -np $NUM_OF_PROCESS ./wifa3d; 
./wifa3dpost; rm RESULTS/*
