NUM_OF_PROCESS=$1

if [ ! $( docker images -a | grep wifa3d | wc -l ) -gt 0 ]; then
    docker build -t wifa3d .
fi

docker run -itd --rm -v .:/app wifa3d bash scripts/run.sh $NUM_OF_PROCESS

