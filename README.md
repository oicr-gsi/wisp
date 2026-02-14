# wisp

WISP ctDNA tumor fraction estimation workflow. The pipeline proceeds through seven major steps: (1) Primary Tumor Characterization - AMBER calculates B-allele frequencies at known heterozygous germline SNP positions, COBALT measures read depth ratios across the genome for copy number analysis, while a pre-computed SAGE VCF provides somatic variants. (2) Variant Annotation - PAVE annotates the SAGE VCF with gene and transcript coding effects, protein impacts, and other functional annotations. (3) Purity and Ploidy Estimation - PURPLE integrates AMBER, COBALT, and PAVE outputs to estimate tumor purity, overall ploidy, and genome-wide copy number segments, establishing the reference somatic profile. (4) Plasma Sample Processing - AMBER and COBALT run on the plasma sample against the matched normal, generating allele frequency and read depth data from cfDNA. (5) Data Integration - AMBER and COBALT outputs from both primary and plasma analyses are merged into unified directories as required by WISP. (6) Variant Force-Calling - SAGE runs in append mode on the plasma sample, force-calling read evidence at variant positions from the primary tumor. (7) Tumor Fraction Estimation - WISP integrates all outputs to estimate ctDNA fraction by examining variant allele frequencies in plasma at known somatic sites. Requirements: whole-genome sequencing data with sufficient coverage; chromosome-subset or targeted panel data will not produce valid estimates.

## Overview

## Dependencies

* [PURPLE](https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md)
* [AMBER](https://github.com/hartwigmedical/hmftools/blob/master/amber/README.md)
* [COBALT](https://github.com/hartwigmedical/hmftools/blob/master/cobalt/README.md)
* [WISP](https://github.com/hartwigmedical/hmftools/tree/master/wisp)
* [SAGE](https://github.com/hartwigmedical/hmftools/tree/master/sage)
* [PAVE](https://github.com/hartwigmedical/hmftools/tree/master/pave)


## Usage

### Cromwell
```
java -jar cromwell.jar run wisp.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`tumour_bam`|File|Input tumor file (bam) of primary sample
`tumour_bai`|File|Input tumor file index (bai) of primary sample
`normal_bam`|File|Input normal file (bam) of primary sample
`normal_bai`|File|Input normal file index (bai) of primary sample
`plasma_bam`|File|Input of bam file of plasma from same donor as primary sample
`plasma_bai`|File|Input of bam index file of plasma from same donor
`donor`|String|Patient identifier
`sage_primary_vcf`|File|Pre-computed SAGE VCF from primary tumor/normal run
`sage_primary_vcf_index`|File|Index for pre-computed SAGE VCF


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`genomeVersion`|String|"hg38"|Genome Version, only 38 supported
`chromosomes`|Array[String]|["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22,chrX,chrY"]|List of chromosomes to process in parallel


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extractTumorName.memory`|Int|4|Memory allocated for this job (GB)
`extractTumorName.heapRam`|Int|1|Heap RAM allocation for GATK (GB)
`extractTumorName.timeout`|Int|4|Hours before task timeout
`extractNormalName.memory`|Int|4|Memory allocated for this job (GB)
`extractNormalName.heapRam`|Int|1|Heap RAM allocation for GATK (GB)
`extractNormalName.timeout`|Int|4|Hours before task timeout
`splitPonByChromosome.memory`|Int|4|Memory allocated for this job (GB)
`splitPonByChromosome.timeout`|Int|10|Hours before task timeout
`amberPrimaryChr.min_mapping_quality`|Int|30|Minimum mapping quality for an alignment to be used
`amberPrimaryChr.min_base_quality`|Int|25|Minimum quality for a base to be considered
`amberPrimaryChr.heapRam`|Int|32|Heap RAM allocation for AMBER (GB)
`amberPrimaryChr.additionalParameters`|String|""|Additional parameters to pass to AMBER
`amberPrimary.memory`|Int|8|Memory allocated for this job (GB)
`amberPrimary.timeout`|Int|10|Hours before task timeout
`cobaltPrimary.gamma`|String|"300"|gamma (penalty) value for segmenting
`cobaltPrimary.min_mapping_quality`|Int|30|Minimum mapping quality for an alignment to be used
`cobaltPrimary.threads`|Int|8|Requested CPU threads
`cobaltPrimary.memory`|Int|32|Memory allocated for this job (GB)
`cobaltPrimary.heapRam`|Int|8|Heap RAM allocation for COBALT (GB)
`cobaltPrimary.timeout`|Int|100|Hours before task timeout
`cobaltPrimary.additionalParameters`|String|""|Additional parameters to pass to COBALT
`pave.driverGenePanel`|String?|None|Driver gene panel file (optional)
`pave.threads`|Int|8|Number of threads
`pave.memory`|Int|32|Memory allocated for this job (GB)
`pave.heapRam`|Int|16|Heap RAM allocation for Pave (GB)
`pave.timeout`|Int|24|Hours before task timeout
`pave.additionalParameters`|String|""|Additional parameters to pass to Pave
`purple.solution_name`|String|"Primary"|Name of solution
`purple.outfilePrefix`|String|tumour_name + ".sol" + solution_name|Prefix of output file
`purple.min_diploid_tumor_ratio_count`|Int|60|smooth over contiguous segments
`purple.min_ploidy`|String?|None|minimum ploidy to consider
`purple.max_ploidy`|String?|None|maximum ploidy to consider
`purple.min_purity`|String?|None|minimum purity to consider
`purple.max_purity`|String?|None|maximum purity to consider
`purple.ploidy_penalty_factor`|String?|None|ploidy penalty factor
`purple.ploidy_penalty_standard_deviation`|String?|None|ploidy penalty standard deviation
`purple.threads`|Int|8|Requested CPU threads
`purple.memory`|Int|32|Memory allocated for this job (GB)
`purple.heapRam`|Int|8|Heap RAM allocation for PURPLE (GB)
`purple.timeout`|Int|100|Hours before task timeout
`purple.additionalParameters`|String|""|Additional parameters to pass to PURPLE
`amberPlasmaChr.min_mapping_quality`|Int|30|Minimum mapping quality for an alignment to be used
`amberPlasmaChr.min_base_quality`|Int|25|Minimum quality for a base to be considered
`amberPlasmaChr.heapRam`|Int|32|Heap RAM allocation for AMBER (GB)
`amberPlasmaChr.additionalParameters`|String|""|Additional parameters to pass to AMBER
`amberPlasma.memory`|Int|8|Memory allocated for this job (GB)
`amberPlasma.timeout`|Int|10|Hours before task timeout
`cobaltPlasma.gamma`|String|"300"|gamma (penalty) value for segmenting
`cobaltPlasma.min_mapping_quality`|Int|30|Minimum mapping quality for an alignment to be used
`cobaltPlasma.threads`|Int|8|Requested CPU threads
`cobaltPlasma.memory`|Int|32|Memory allocated for this job (GB)
`cobaltPlasma.heapRam`|Int|8|Heap RAM allocation for COBALT (GB)
`cobaltPlasma.timeout`|Int|100|Hours before task timeout
`cobaltPlasma.additionalParameters`|String|""|Additional parameters to pass to COBALT
`sage_append.driverGenePanel`|String?|None|Driver gene panel TSV (optional)
`sage_append.min_map_quality`|Int|10|Minimum mapping quality
`sage_append.hard_min_tumor_qual`|Int|50|Hard minimum tumor quality
`sage_append.hard_min_tumor_raw_alt_support`|Int|2|Minimum raw alt support
`sage_append.hard_min_tumor_vaf`|Float|0.002|Minimum tumor VAF
`sage_append.threads`|Int|8|Requested CPU threads
`sage_append.memory`|Int|40|Memory allocated for this job (GB)
`sage_append.heapRam`|Int|32|Heap RAM allocation for SAGE (GB)
`sage_append.timeout`|Int|24|Hours before task timeout
`sage_append.additionalParameters`|String|""|Additional parameters to pass to SAGE
`annotatePrimaryVcfWithPurple.modules`|String|"bcftools/1.9"|Required environment modules
`generateProbeVariants.outputFileName`|String|"probe_variants.csv"|Output CSV filename for probe variants
`generateProbeVariants.jobMemory`|Int|4|Memory allocated for this job (GB)
`generateProbeVariants.timeout`|Int|1|Hours before task timeout
`runWisp.additionalParameters`|String|""|Additional parameters to pass to WISP
`runWisp.threads`|Int|4|Requested CPU threads
`runWisp.heapRam`|Int|16|Heap RAM allocation for WISP (GB)
`runWisp.memory`|Int|16|Memory allocated for this job (GB)
`runWisp.timeout`|Int|50|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`wisp_directory`|File|Zipped WISP output directory|
`wisp_summary`|File|WISP tumor fraction summary|
`wisp_somatic_variants`|File|Per-variant SNV results|
`sage_plasma_vcf`|File|SAGE VCF with plasma sample appended|
`sage_plasma_vcf_index`|File|Index for SAGE plasma VCF|

## Commands
This section lists command(s) run by WORKFLOW workflow

* Running WORKFLOW


```
    set -euo pipefail

    if [ -f "~{inputBam}" ]; then
      gatk --java-options "-Xmx~{heapRam}g" GetSampleName -R ~{refFasta} -I ~{inputBam} -O input_name.txt -encode
    fi

    cat input_name.txt
```
```
    set -euo pipefail
    
    # Extract just this chromosome from PON
    bcftools view -r ~{chromosome} ~{PON} -O z -o pon.~{chromosome}.vcf.gz
    tabix -p vcf pon.~{chromosome}.vcf.gz
```
```
    set -euo pipefail

    mkdir ~{file_prefix}.amber
    genome_version=$(if [ "~{genomeVersion}" == "hg38" ]; then echo "38"; else echo "19"; fi)

    java -Xmx~{heapRam}G -cp $HMFTOOLS_ROOT/amber.jar com.hartwig.hmftools.amber.AmberApplication \
      -reference ~{normal_name} -reference_bam ~{normal_bam} \
      -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
      -output_dir ~{file_prefix}.amber/ \
      -loci ~{PON} \
      -ref_genome_version ${genome_version} \
      -min_mapping_quality ~{min_mapping_quality} \
      -min_base_quality ~{min_base_quality} \
      ~{additionalParameters}

    zip -r ~{file_prefix}.amber.zip ~{file_prefix}.amber/

```
```
    set -euo pipefail

    mkdir -p ~{tumour_name}.amber
    mkdir -p temp

    # Unzip all chromosome results
    for zip_file in ~{sep=' ' chr_zips}; do
      unzip -o ${zip_file} -d temp/
    done

    first_baf=$(ls temp/*/*.amber.baf.tsv.gz | sort -V | sed -n '1p')
    zcat "$first_baf" | sed -n '1p' > temp_header.txt
    zcat temp/*/*.amber.baf.tsv.gz | grep -v "^chromosome" | sort -k1,1V -k2,2n > temp_data.txt
    cat temp_header.txt temp_data.txt | gzip > ~{tumour_name}.amber/~{tumour_name}.amber.baf.tsv.gz
    rm temp_header.txt temp_data.txt

    first_pcf=$(ls temp/*/*.amber.baf.pcf | sort -V | sed -n '1p')
    sed -n '1p' "$first_pcf" > temp_header_pcf.txt
    tail -q -n +2 temp/*/*.amber.baf.pcf | sort -k2,2V -k4,4n > temp_data_pcf.txt
    cat temp_header_pcf.txt temp_data_pcf.txt > ~{tumour_name}.amber/~{tumour_name}.amber.baf.pcf
    rm temp_header_pcf.txt temp_data_pcf.txt

    cp $(find temp -name "*.amber.qc" -type f | sed -n '1p') ~{tumour_name}.amber/~{tumour_name}.amber.qc

    if ls temp/*/*.amber.contamination.tsv 1> /dev/null 2>&1; then
      {
        head -n1 $(ls temp/*/*.amber.contamination.tsv | sort -V | head -n1)
        tail -q -n +2 temp/*/*.amber.contamination.tsv | sort -k1,1V -k2,2n
      } > ~{tumour_name}.amber/~{tumour_name}.amber.contamination.tsv
    fi

    if ls temp/*/*.amber.contamination.vcf.gz 1> /dev/null 2>&1; then
      bcftools concat $(ls temp/*/*.amber.contamination.vcf.gz | sort -V) | \
        bcftools sort -O z -o ~{tumour_name}.amber/~{tumour_name}.amber.contamination.vcf.gz
      tabix -p vcf ~{tumour_name}.amber/~{tumour_name}.amber.contamination.vcf.gz
    fi

    if ls temp/*/*.amber.homozygousregion.tsv 1> /dev/null 2>&1; then
      {
        head -n1 $(ls temp/*/*.amber.homozygousregion.tsv | sort -V | head -n1)
        tail -q -n +2 temp/*/*.amber.homozygousregion.tsv | sort -k1,1V -k2,2n
      } > ~{tumour_name}.amber/~{tumour_name}.amber.homozygousregion.tsv
    fi

    if find temp -name "amber.version" -type f | grep -q .; then
      cp $(find temp -name "amber.version" -type f | head -n1) ~{tumour_name}.amber/amber.version
    fi

    zip -r ~{tumour_name}.amber.zip ~{tumour_name}.amber/
```
```
    set -euo pipefail

    mkdir ~{tumour_name}.cobalt

    java -Xmx~{heapRam}G -cp $HMFTOOLS_ROOT/cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication \
      -reference ~{normal_name} -reference_bam ~{normal_bam} \
      -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
      -output_dir ~{tumour_name}.cobalt/ \
      -gc_profile ~{gcProfile} \
      -pcf_gamma ~{gamma} \
      -min_quality ~{min_mapping_quality} \
      ~{additionalParameters}

    zip -r ~{tumour_name}.cobalt.zip ~{tumour_name}.cobalt/

```
```
    mkdir -p ~{output_name}
    unzip -o ~{primary_zip} -d ~{output_name}/
    unzip -o ~{plasma_zip} -d ~{output_name}/
    # Flatten if nested
    find ~{output_name} -mindepth 2 -type f -exec mv {} ~{output_name}/ \;
    find ~{output_name} -mindepth 1 -type d -empty -delete
    zip -r ~{output_name}.zip ~{output_name}/
```
```
    set -euo pipefail

    echo $SAGE_ROOT
    ls -la ${SAGE_ROOT}/

    java -Xmx~{heapRam}G -cp $SAGE_ROOT/sage.jar com.hartwig.hmftools.sage.append.SageAppendApplication \
      -reference ~{reference_name} \
      -reference_bam ~{reference_bam} \
      -ref_genome_version 38 \
      -ref_genome ~{refFasta} \
      ~{"-driver_gene_panel " + driverGenePanel} \
      -input_vcf ~{input_vcf} \
      -output_vcf ~{reference_name}.sage.vcf.gz \
      -threads ~{threads} \
      -min_map_quality ~{min_map_quality} \
      -hard_min_tumor_qual ~{hard_min_tumor_qual} \
      -hard_min_tumor_raw_alt_support ~{hard_min_tumor_raw_alt_support} \
      -hard_min_tumor_vaf ~{hard_min_tumor_vaf} \
      -skip_msi_jitter \
      ~{additionalParameters}

    mkdir -p ~{reference_name}.sage.bqr/
    mv *.sage.bqr.tsv  ~{reference_name}.sage.bqr/
    zip -r ~{reference_name}.sage.bqr.zip ~{reference_name}.sage.bqr/

```
```
    set -euo pipefail

    # Unzip input directories
    unzip ~{purple_dir}
    unzip ~{amber_dir}
    unzip ~{cobalt_dir}
    unzip ~{bqr_dir}

    # Create output directory
    mkdir -p ~{plasma_name}.wisp

    # Run WISP with skip_bqr
    java -Xmx~{heapRam}G -jar $WISP_ROOT/wisp.jar \
      -patient_id ~{donor} \
      -tumor_id ~{tumour_name} \
      -samples ~{plasma_name} \
      -purple_dir ~{tumour_name}.solPrimary.purple/ \
      -amber_dir ~{tumour_name}.amber/ \
      -cobalt_dir ~{tumour_name}.cobalt/ \
      -bqr_dir ~{plasma_name}.sage.bqr \
      -somatic_vcf ~{somatic_vcf} \
      -ref_genome ~{refFasta} \
      -output_dir ~{plasma_name}.wisp/ \
      -probe_variants_file ~{probe_variants_file} \
      -threads ~{threads} \
      -write_types ALL \
      ~{additionalParameters}

    # Zip output
    zip -r ~{plasma_name}.wisp.zip ~{plasma_name}.wisp/
    cp ~{plasma_name}.wisp/~{donor}_~{plasma_name}.wisp.summary.tsv ~{plasma_name}.wisp.summary.tsv 
    cp ~{plasma_name}.wisp/~{donor}_~{plasma_name}.wisp.somatic_variants.tsv ~{plasma_name}.wisp.somatic_variants.tsv

```
```
    set -euo pipefail

    unzip ~{amber_directory}
    unzip ~{cobalt_directory}
    mkdir ~{outfilePrefix}.purple
    genome_version=$(if [ "~{genomeVersion}" == "hg38" ]; then echo "38"; else echo "19"; fi)

    java -Xmx~{heapRam}G -jar $HMFTOOLS_ROOT/purple.jar \
      -ref_genome_version ${genome_version} \
      -ref_genome ~{refFasta}  \
      -gc_profile ~{gcProfile} \
      -ensembl_data_dir ~{ensemblDir}  \
      -reference ~{normal_name} -tumor ~{tumour_name}  \
      -amber ~{tumour_name}.amber -cobalt ~{tumour_name}.cobalt \
      ~{"-ploidy_penalty_factor " + ploidy_penalty_factor} \
      ~{"-ploidy_penalty_standard_deviation " + ploidy_penalty_standard_deviation} \
      ~{"-somatic_vcf " + somatic_vcf} \
      ~{"-min_ploidy " + min_ploidy} \
      ~{"-max_ploidy " + max_ploidy} \
      ~{"-min_purity " + min_purity} \
      ~{"-max_purity " + max_purity} \
      -no_charts \
      -min_diploid_tumor_ratio_count ~{min_diploid_tumor_ratio_count} \
      -output_dir ~{outfilePrefix}.purple \
      ~{additionalParameters}

    zip -r ~{outfilePrefix}.purple.zip ~{outfilePrefix}.purple/

```
```
    # Extract Purple INFO annotations and merge into plasma VCF
    bcftools annotate \
      -a ~{purple_vcf} \
      -c INFO/SUBCL,INFO/PURPLE_VCN,INFO/PURPLE_AF,INFO/PURPLE_CN \
      ~{sage_vcf} \
      -Oz -o annotated.vcf.gz
    
    tabix -p vcf annotated.vcf.gz
```
```
    set -euo pipefail

    # Generate probe variants CSV from PASS + HIGH_CONFIDENCE variants
    zcat ~{annotated_vcf} | \
      grep -v "^#" | \
      awk 'BEGIN{print "TumorId,Category,Variant"}
           $7=="PASS" && $8~"TIER=HIGH_CONFIDENCE" {
             printf "~{tumor_id},REPORTABLE_MUTATION,%s:%s %s>%s\n", $1, $2, $4, $5
           }' > ~{outputFileName}

    # Verify file has content
    line_count=$(wc -l < ~{outputFileName})
    if [ "$line_count" -le 1 ]; then
      echo "WARNING: No PASS + HIGH_CONFIDENCE variants found"
      echo "This may be expected for downsampled/test data"
    fi

    echo "Generated probe variants file with $((line_count - 1)) variants"
```
```
    set -euo pipefail

    genome_version=$(if [ "~{genomeVersion}" == "hg38" ]; then echo "38"; else echo "37"; fi)

    java -Xmx~{heapRam}G -jar $PAVE_ROOT/pave.jar \
      -sample ~{sample} \
      -input_vcf ~{vcf_file} \
      -ref_genome ~{refFasta} \
      -ref_genome_version ${genome_version} \
      -ensembl_data_dir ~{ensemblDir} \
      ~{"-driver_gene_panel " + driverGenePanel} \
      -output_dir ./ \
      -threads ~{threads} \
      ~{additionalParameters}
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
