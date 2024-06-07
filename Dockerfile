# Use a base image with the desired dependencies
FROM ubuntu:22.04 as builder
# FROM ubuntu:latest as builder

# Installations based on saragiuliani/mum-phinder
RUN apt-get update -qq && \
    apt-get install -y zlib1g-dev \
                    git \
                    cmake \
                    build-essential \
                    python3 \
                    gcc-9 \
                    g++-9 \
                    seqtk \
                    samtools \
                    && \
    update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-9 9 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-9 9

# Clone and build MONI
RUN git clone https://github.com/StephenHwang/moni.git && \
    cd moni && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    make install
# Add the build directory to the PATH
ENV PATH="/moni/build:${PATH}"


# Clone MEMO
RUN git clone https://github.com/StephenHwang/MEMO.git 
# Add the build directory to the PATH
ENV PATH="/MEMO/src:${PATH}"
