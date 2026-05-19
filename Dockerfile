#FROM mambaorg/micromamba:2.0.5-cuda12.4.1-ubuntu22.04
FROM mambaorg/micromamba@sha256:1f64e993d58570b4c4fa8889cd7518720586e4ef86d8624da61769f96dac29b5

USER root

RUN --mount=type=cache,target=/var/cache/apt \
    apt-get update -y &&\
    apt-get -qq -y install curl && \
    apt-get -y install python3-dev && \
    apt-get -y install gcc && \
    apt-get install -y zlib1g-dev && \
    apt-get install -y augustus && \
    apt-get install -y augustus-data augustus-doc && \
    apt-get install -y git

ARG MAMBA_DOCKERFILE_ACTIVATE=1 

COPY ./jump_orf.yaml /tmp/env.yaml

RUN --mount=type=cache,target=/opt/conda/pkgs \
    micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

RUN mkdir /module
RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install URI::Escape
WORKDIR /module
RUN git clone --branch TransDecoder-v5.7.1 --depth 1 https://github.com/TransDecoder/TransDecoder

WORKDIR /
COPY ./module /module
COPY ./hmm /hmm
COPY ./JumpORF.py /usr/local/bin/JumpORF
RUN chmod +x /usr/local/bin/JumpORF
RUN chmod +x /module/run_LTR_RTE_detection.py
RUN chmod +x /module/run_DNA_TE_detection.py
RUN chmod +x /module/run_Helitron_detection.py
RUN chmod +x /module/run_LINE_detection.py
RUN chmod +x /module/run_DIRS_detection.py
RUN chmod +x /module/run_PLE_detection.py
RUN chmod +x /module/run_aORF_detection.py

WORKDIR /home

# Set the compiler and library paths to use the system defaults
#ENV CC=/usr/bin/gcc
#ENV CXX=/usr/bin/g++
#ENV CPLUS_INCLUDE_PATH=/usr/include:/usr/include/x86_64-linux-gnu
#ENV LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:/usr/lib

SHELL ["/bin/bash", "-c"]

EXPOSE 8888
CMD /bin/bash

