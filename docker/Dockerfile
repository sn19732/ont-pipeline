FROM ubuntu:16.04
MAINTAINER Botond.Sipos@nanoporetech.com

RUN apt-get update && apt-get install -y sudo git make build-essential mummer python-numpy python-matplotlib time bwa-mem samtools software-properties-common
RUN add-apt-repository ppa:webupd8team/java &&\
    apt-get update &&\
    (echo "oracle-java8-installer shared/accepted-oracle-license-v1-1 select true" | sudo debconf-set-selections) &&\
     apt-get install -y oracle-java8-installer
