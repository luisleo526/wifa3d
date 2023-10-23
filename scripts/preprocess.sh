mpif90 -o /app/run/wifa3dhead src/wifa3dhead.f -O3 -mcmodel=large -fno-align-commons
cd /app/run; ./wifa3dhead
cp /app/run/head.inc /app/src/
