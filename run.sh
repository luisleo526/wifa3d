ROOT=$1
NUM_OF_PROCESS=$2

if [ ! $( docker images -a | grep wifa3d | wc -l ) -gt 0 ]; then
    docker build -t wifa3d .
fi

docker run -it --rm --user $UID  -v .:/app wifa3d bash /app/entrypoint.sh $ROOT $NUM_OF_PROCESS

