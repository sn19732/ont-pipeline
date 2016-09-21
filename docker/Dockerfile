FROM ubuntu:16.04
MAINTAINER Botond.Sipos@nanoporetech.com

# Upgrade and install necessary packages:
RUN apt-get update && apt-get upgrade &&\
    apt-get install -y sudo git make build-essential mummer python-numpy python-matplotlib time bwa samtools software-properties-common

# Install latest Java:
RUN add-apt-repository ppa:webupd8team/java &&\
    apt-get update &&\
    (echo "oracle-java8-installer shared/accepted-oracle-license-v1-1 select true" | sudo debconf-set-selections) &&\
     apt-get install -y oracle-java8-installer

# Start working in /opt
WORKDIR /opt

# Compile and install canu:
RUN git clone clone https://github.com/marbl/canu.git &&\
    cd canu/src && make -j && cp /opt/canu/Linux-amd64/bin/canu/* /usr/local/bin/ && cd - &&\
    rm -r canu
