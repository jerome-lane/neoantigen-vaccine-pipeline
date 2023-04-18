#!/bin/bash

sudo systemctl start docker && sudo systemctl status docker

WORKDIR='/data/benchmarking'

OLD_IMAGE_TAG='jeromelane/neoantigen-discovery-pipeline:latest'
NEW_IMAGE_TAG='jeromelane/neoantigen-discovery-pipeline-feature_update_mutect_strelka_variant_callers:latest'
CONFIGFILE=BENCHMARK_CONFIG.yaml

CONTAINER_NAME=OLD_VARIANT_CALLERS
# Create openvax container
docker run \
--entrypoint /bin/bash \
--name $CONTAINER_NAME \
-it \
-v "$WORKDIR/inputs":/inputs \
-v "$WORKDIR/old/$OUTPUTS":/outputs \
-v "/data/reference-data":/reference-genome \
-d "$OLD_IMAGE_TAG"

docker exec -t $CONTAINER_NAME sh -c "python run_snakemake.py --configfile=/inputs/$CONFIGFILE --run-qc --somatic-variant-calling-only >> /outputs/global_snakemake.log"

CONTAINER_NAME=NEW_VARIANT_CALLERS

# Create openvax container
docker run \
--entrypoint /bin/bash \
--name $CONTAINER_NAME \
-it \
-v "$WORKDIR/inputs":/inputs \
-v "$WORKDIR/new/$OUTPUTS":/outputs \
-v "$WORKDIR/reference-data/grch38":/reference-genome \
-d "$NEW_IMAGE_TAG"

docker exec -t $CONTAINER_NAME sh -c "python run_snakemake.py --configfile=/inputs/$CONFIGFILE --run-qc --somatic-variant-calling-only >> /outputs/global_snakemake.log"
