#!/bin/bash

sudo systemctl start docker && sudo systemctl status docker

WORKDIR='/data/benchmarking'
VARIANT_VCF_COPY_FOLDER="$WORKDIR/variant_calls"
PATIENT_ID='HCC_1395'

OLD_IMAGE_TAG='jeromelane/neoantigen-discovery-pipeline:latest'
#NEW_IMAGE_TAG='jeromelane/neoantigen-discovery-pipeline-feature_update_mutect_strelka_variant_callers:latest'
NEW_IMAGE_TAG='neoantigen-vaccine-pipeline-new'

CONFIGFILE=BENCHMARK_CONFIG.yaml

CONTAINER_NAME=OLD_VARIANT_CALLERS
# Create openvax container
docker run \
--entrypoint /bin/bash \
--name $CONTAINER_NAME \
-it \
-v "$WORKDIR/inputs":/inputs \
-v "$WORKDIR/old":/outputs \
-v "/data/reference-data":/reference-genome \
-d "$OLD_IMAGE_TAG"

docker exec -t $CONTAINER_NAME sh -c "python run_snakemake.py --configfile=/inputs/$CONFIGFILE --run-qc --somatic-variant-calling-only >> /outputs/global_snakemake.log"

cp -r $WORKDIR/old/$PATIENT_ID/mutect.vcf $WORKDIR/old/$PATIENT_ID/strelka.vcf $WORKDIR/old/$PATIENT_ID/strelka_output/ $VARIANT_VCF_COPY_FOLDER/old/
rm -r $WORKDIR/old/$PATIENT_ID/mutect.vcf $WORKDIR/old/$PATIENT_ID/strelka.vcf $WORKDIR/old/$PATIENT_ID/strelka_output/

CONTAINER_NAME=NEW_VARIANT_CALLERS

# Create openvax container
docker run \
--entrypoint /bin/bash \
--name $CONTAINER_NAME \
-it \
-v "$WORKDIR/inputs":/inputs \
-v "$WORKDIR/old":/outputs \
-v "/data/reference-data":/reference-genome \
-d "$NEW_IMAGE_TAG"

docker exec -t $CONTAINER_NAME sh -c "python run_snakemake.py --configfile=/inputs/$CONFIGFILE --run-qc --somatic-variant-calling-only >> /outputs/global_snakemake.log"

cp -r $WORKDIR/old/$PATIENT_ID/mutect.vcf $WORKDIR/old/$PATIENT_ID/strelka.vcf $WORKDIR/old/$PATIENT_ID/strelka_output/ $VARIANT_VCF_COPY_FOLDER/new/
