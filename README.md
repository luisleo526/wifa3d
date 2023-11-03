# Requirements
- [Docker](https://docs.docker.com/engine/install/)

and run the command `docker build -t wifa3d .` if you have not done before.

# Command
Make sure the input files (`.chk, .ctl, .grd, .pro, .wtd`) placed in YOUR `$run` directory, with same name stored in `$run/CASE`
and execute the following command:

```bash
bash run.sh $run $NUM_OF_PROCESS
```

`bash run.sh sample-input 32` for example.
