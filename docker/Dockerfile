# from ubuntu image
#FROM ubuntu
FROM rocker/rstudio

# update apt
RUN apt-get update

# copy scripts req. for docker
COPY . . 

# ensure apt is not interactive
ARG DEBIAN_FRONTEND=noninteractive
# disable auth for rstudio server
#ENV DISABLE_AUTH=true
# root permissions
#ENV ROOT=TRUE

# run install script
RUN ./install.sh --noninteractive
RUN Rscript r_pkgs.R

#ENV USER="rstudio"
#CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize", "0", "--auth-none", "1"]