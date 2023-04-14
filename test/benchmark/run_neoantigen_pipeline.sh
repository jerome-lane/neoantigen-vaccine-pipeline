#!/bin/bash

sudo systemctl start docker && sudo systemctl status docker

WORDIR='/data/benchmarking'

OLD_IMAGE_TAG='neoantigen-vaccine-pipeline:latest'
NEW_IMAGE_TAG='neoantigen-vaccine-pipeline-new:latest'
CONFIGFILE=BENCHMARK_CONFIG.yaml

CONTAINER_NAME=OLD_VARIANT_CALLERS
# Create openvax container
docker run \
--entrypoint /bin/bash \
--name $CONTAINER_NAME \
-it \
-v "$WORDIR/inputs":/inputs \
-v "$WORDIR/old/$OUTPUTS":/outputs \
-v "$WORDIR/reference-data/grch38":/reference-genome \
-d "$OLD_IMAGE_TAG"

docker exec -t $CONTAINER_NAME sh -c "python run_snakemake.py --configfile=/inputs/$CONFIGFILE --run-qc --somatic-variant-only >> /outputs/global_snakemake.log"

CONTAINER_NAME=NEW_VARIANT_CALLERS

# Create openvax container
docker run \
--entrypoint /bin/bash \
--name $CONTAINER_NAME \
-it \
-v "$WORDIR/inputs":/inputs \
-v "$WORDIR/new/$OUTPUTS":/outputs \
-v "$WORDIR/reference-data/grch38":/reference-genome \
-d "$NEW_IMAGE_TAG"

docker exec -t $CONTAINER_NAME sh -c "python run_snakemake.py --configfile=/inputs/$CONFIGFILE --run-qc --somatic-variant-only >> /outputs/global_snakemake.log"
