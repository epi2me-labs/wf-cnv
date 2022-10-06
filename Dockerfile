ARG BASEIMAGE=ontresearch/base-workflow-image:v0.2.0
FROM $BASEIMAGE
ARG ENVFILE=environment.yaml

COPY $ENVFILE $HOME/environment.yaml
RUN \
    . $CONDA_DIR/etc/profile.d/micromamba.sh \
    && micromamba activate \
    && micromamba install -n base --file $HOME/environment.yaml \
    && fix-permissions $CONDA_DIR \
    && fix-permissions $HOME \
    && Rscript -e 'requireNamespace("BiocManager"); BiocManager::install("QDNAseq.hg19")' \
    && R -e "remotes::install_github('asntech/QDNAseq.hg38@main')"



USER $WF_UID
WORKDIR $HOME
