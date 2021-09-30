#!/bin/sh
# Affymetrix
tabix -h ${HOME}/data/godb/affy/chr${CHR}.vcf.gz ${TCHR}:${STARTP}-${ENDP} > ${HOME}/devel/godb/data/affy/chr${CHR}.vcf
bgzip -f ${HOME}/devel/godb/data/affy/chr${CHR}.vcf
tabix ${HOME}/devel/godb/data/affy/chr${CHR}.vcf.gz
gzcat ${HOME}/devel/godb/data/affy/chr${CHR}.vcf.gz | grep -v '#' | cut -f1-9 > ${HOME}/devel/godb/data/affy/chr${CHR}_prfx.txt
cat ${HOME}/devel/godb/data/affy/chr${CHR}_prfx.txt | cut -f3 > ${HOME}/devel/godb/data/affy/chr${CHR}_snps.txt
# Broad
tabix -h ${HOME}/data/godb/broad/chr${CHR}.vcf.gz ${TCHR}:${STARTP}-${ENDP} > ${HOME}/devel/godb/data/broad/chr${CHR}.vcf
bgzip -f ${HOME}/devel/godb/data/broad/chr${CHR}.vcf
tabix ${HOME}/devel/godb/data/broad/chr${CHR}.vcf.gz
gzcat ${HOME}/devel/godb/data/broad/chr${CHR}.vcf.gz | grep -v '#' | cut -f1-9 > ${HOME}/devel/godb/data/broad/chr${CHR}_prfx.txt
cat ${HOME}/devel/godb/data/broad/chr${CHR}_prfx.txt | cut -f3 > ${HOME}/devel/godb/data/broad/chr${CHR}_snps.txt
# Exome
tabix -h ${HOME}/data/godb/exome/chr${CHR}.vcf.gz ${TCHR}:${STARTP}-${ENDP} > ${HOME}/devel/godb/data/exome/chr${CHR}.vcf
bgzip -f ${HOME}/devel/godb/data/exome/chr${CHR}.vcf
tabix ${HOME}/devel/godb/data/exome/chr${CHR}.vcf.gz
gzcat ${HOME}/devel/godb/data/exome/chr${CHR}.vcf.gz | grep -v '#' | cut -f1-9 > ${HOME}/devel/godb/data/exome/chr${CHR}_prfx.txt
cat ${HOME}/devel/godb/data/exome/chr${CHR}_prfx.txt | cut -f3 > ${HOME}/devel/godb/data/exome/chr${CHR}_snps.txt
# Illumina
tabix -h ${HOME}/data/godb/illumina/chr${CHR}.vcf.gz ${TCHR}:${STARTP}-${ENDP} > ${HOME}/devel/godb/data/illumina/chr${CHR}.vcf
bgzip -f ${HOME}/devel/godb/data/illumina/chr${CHR}.vcf
tabix ${HOME}/devel/godb/data/illumina/chr${CHR}.vcf.gz
gzcat ${HOME}/devel/godb/data/illumina/chr${CHR}.vcf.gz | grep -v '#' | cut -f1-9 > ${HOME}/devel/godb/data/illumina/chr${CHR}_prfx.txt
cat ${HOME}/devel/godb/data/illumina/chr${CHR}_prfx.txt | cut -f3 > ${HOME}/devel/godb/data/illumina/chr${CHR}_snps.txt
# Metabochip
tabix -h ${HOME}/data/godb/metabo/chr${CHR}.vcf.gz ${TCHR}:${STARTP}-${ENDP} > ${HOME}/devel/godb/data/metabo/chr${CHR}.vcf
bgzip -f ${HOME}/devel/godb/data/metabo/chr${CHR}.vcf
tabix ${HOME}/devel/godb/data/metabo/chr${CHR}.vcf.gz
gzcat ${HOME}/devel/godb/data/metabo/chr${CHR}.vcf.gz | grep -v '#' | cut -f1-9 > ${HOME}/devel/godb/data/metabo/chr${CHR}_prfx.txt
cat ${HOME}/devel/godb/data/metabo/chr${CHR}_prfx.txt | cut -f3 > ${HOME}/devel/godb/data/metabo/chr${CHR}_snps.txt
