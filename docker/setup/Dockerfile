# from ubuntu image
#FROM ubuntu
FROM rocker/rstudio

# update apt
RUN apt-get update

# copy scripts req. for docker
#COPY install.sh . 
#COPY r_pkgs.R .

# add aliases
RUN cat aliases >> ~/.bashrc
RUN cat aliases >> /home/rstudio/.bashrc

# ensure apt is not interactive
ARG DEBIAN_FRONTEND=noninteractive

# run install script
RUN ./install.sh --noninteractive
RUN Rscript r_pkgs.R

# add rstudio user group to sudoers
RUN usermod -aG sudo rstudio

