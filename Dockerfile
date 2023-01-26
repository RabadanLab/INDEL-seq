FROM ubuntu:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/fastq-multx:/opt/pear-0.9.11-linux-x86_64/bin:$PATH
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt update
RUN apt install -y build-essential \
                    git \
                    python3-pandas \
                    cutadapt \
                    emboss

#Install fastq-multx
RUN cd /opt/ && git clone https://github.com/brwnj/fastq-multx && cd fastq-multx && make

#Install PEAR
COPY pear-0.9.11-linux-x86_64.tar.gz /opt/
RUN cd /opt/ && tar -xvf /opt/pear-0.9.11-linux-x86_64.tar.gz
