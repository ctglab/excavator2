FROM ghcr.io/mamba-org/micromamba:latest
USER root
RUN apt-get update && apt-get install -y curl

COPY excavator2.yml /tmp/excavator2.yml
RUN micromamba env create -f /tmp/excavator2.yml && \
    micromamba clean --all --yes

# COPY your custom entrypoint script and make it executable
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

COPY excavator2 /opt/excavator2
WORKDIR /opt/excavator2

ENV PATH="/opt/excavator2:${PATH}"

# Set the ENTRYPOINT to your new script
ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]

# CMD provides the default command to the ENTRYPOINT
CMD ["/opt/excavator2/excavator2.sh"]
