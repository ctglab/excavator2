FROM ghcr.io/mamba-org/micromamba:latest
USER root
RUN apt-get update && apt-get install -y curl
COPY excavator2.yml /tmp/excavator2.yml
RUN micromamba env create -f /tmp/excavator2.yml && \
    micromamba clean --all --yes
RUN echo 'eval "$(micromamba shell hook --shell bash)"' >> /root/.bashrc
RUN echo 'micromamba activate excavator2' >> /root/.bashrc
RUN echo 'export PATH="/opt/excavator2:${PATH}"' >> /root/.bashrc
COPY excavator2 /opt/excavator2
WORKDIR /opt/excavator2
ENV PATH="/opt/excavator2:${PATH}"
CMD ["/opt/excavator2/excavator2.sh"]