ROOT=$(basename $1)
NUM_OF_PROCESS=$2

find /app/$ROOT -type f -name "*.pro" -exec sh -c 'echo "${0}" > "${1}"' $NUM_OF_PROCESS {} \;
cp -r /app/src /app/$ROOT

cd /app/$ROOT
mkdir -p RESULT
mpif90 -o wifa3dhead src/wifa3dhead.f -O3 -mcmodel=large -fno-align-commons
./wifa3dhead
cp head.inc ./src/

mpif90 -o wifa3d ./src/wifa3d.f -O3 -mcmodel=large -fno-align-commons
mpif90 -o wifa3dpost ./src/wifa3dpost.f -O3 -mcmodel=large -fno-align-commons
mpirun -np $NUM_OF_PROCESS ./wifa3d; 
./wifa3dpost; rm RESULT/*
