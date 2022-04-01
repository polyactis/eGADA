FROM ubuntu:20.04
LABEL AUTHOR="Yu Huang <polyactis@gmail.com>"
LABEL DESCRIPTION="https://github.com/polyactis/eGADA"
ADD bashrc /root/.bashrc
ENV DEBIAN_FRONTEND="noninteractive"
# apt-get will pop some dialogs. This disable them. 
# https://stackoverflow.com/questions/40890011/ubuntu-dockerfile-mailutils-install
RUN apt-get update
RUN apt-get install -y tzdata
# set your timezone
ENV TZ=Asia/Shanghai
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN dpkg-reconfigure --frontend noninteractive tzdata

RUN apt-get install -y apt-utils iputils-ping less lynx \
    libboost-program-options1.71.0 libboost-iostreams1.71.0 libboost-python1.71.0\
    man python3-minimal telnet tmux vim
# install pyflow, use absolute path. Default path is /.
# do not split it in two commands (cd and then run this. as every command starts from /)

ENV HOME /root
CMD source /root/.bashrc
