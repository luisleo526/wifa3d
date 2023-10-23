FROM ubuntu:20.04
WORKDIR /app
RUN apt update -y && apt upgrade -y 
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt install build-essential gfortran cmake mpich -y
