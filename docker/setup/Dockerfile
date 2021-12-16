# from rstudio-server image
FROM rocker/rstudio

# update apt
RUN apt-get update

# copy scripts to /install dir
RUN mkdir /install
COPY aliases /install
COPY install.sh /install
COPY r_pkgs.R /install

# add aliases
RUN cat /install/aliases >> ~/.bashrc
RUN cat /install/aliases >> /home/rstudio/.bashrc

# ensure apt is not interactive
ARG DEBIAN_FRONTEND=noninteractive

# run install scripts
RUN bash /install/install.sh
RUN Rscript /install/r_pkgs.R

# add rstudio user group to sudoers
RUN usermod -aG sudo rstudio