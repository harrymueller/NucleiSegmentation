# from ubuntu image
FROM ubuntu

# update apt
RUN apt-get update

# copy scripts req. for docker
COPY . ~/tools 

# ensure apt is not interactive
ARG DEBIAN_FRONTEND=noninteractive

# run install script
RUN ./install.sh --noninteractive