version 1.0

struct GenomeResources {
    String PON
    String wispModules
    String hmfModules
    String gatkModules
    String gcProfile
    String ensemblDir
    String refFasta
    String refFai
    String pon_sgl_file
    String pon_sv_file
    String known_hotspot_file
    String repeat_mask_file
    String knownfusion
    String hotspots
    String driverGenePanel
    String highConfBed
}

workflow wisp {
  input {
    File tumour_bam
    File tumour_bai
    File normal_bam
    File normal_bai
    File plasma_bam
    File plasma_bai
    String donor
    String genomeVersion = "38"
    File sage_primary_vcf
    File sage_primary_vcf_index
    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22,chrX,chrY"]
  }

  parameter_meta {
    tumour_bam: "Input tumor file (bam) of primary sample"
    tumour_bai: "Input tumor file index (bai) of primary sample"
    normal_bam: "Input normal file (bam) of primary sample"
    normal_bai: "Input normal file index (bai) of primary sample"
    plasma_bam: "Input of bam file of plasma from same donor as primary sample"
    plasma_bai: "Input of bam index file of plasma from same donor"
    donor: "Patient identifier"
    genomeVersion: "Genome Version, only 38 supported"
    sage_primary_vcf: "Pre-computed SAGE VCF from primary tumor/normal run"
    sage_primary_vcf_index: "Index for pre-computed SAGE VCF"
    chromosomes: "List of chromosomes to process in parallel"
  }

  Map[String,GenomeResources] resources = {
    "38": {
      "wispModules": "wisp/v1.2-beta.3 hg38/p12",
      "hmfModules": "hmftools/1.1 hg38/p12 hmftools-data/53138",
      "gatkModules": "hg38-gridss-index/1.0 gatk/4.1.6.0",
      "refFasta": "$HG38_ROOT/hg38_random.fa",
      "refFai": "$HG38_GRIDSS_INDEX_ROOT/hg38_random.fa.fai",
      "PON" : "$SAGE_DATA_ROOT/GermlineHetPon.38.vcf.gz",
      "ensemblDir": "$HMFTOOLS_DATA_ROOT/ensembl_data",
      "gcProfile": "$HMFTOOLS_DATA_ROOT/copy_number/GC_profile.1000bp.38.cnp",
      "pon_sgl_file": "$HMFTOOLS_DATA_ROOT/sv/sgl_pon.38.bed.gz",
      "pon_sv_file": "$HMFTOOLS_DATA_ROOT/sv/sv_pon.38.bedpe.gz",
      "known_hotspot_file": "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe",
      "repeat_mask_file": "$HMFTOOLS_DATA_ROOT/sv/repeat_mask_data.38.fa.gz",
      "knownfusion": "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe",
      "hotspots": "$SAGE_DATA_ROOT/KnownHotspots.hg38.fixed.vcf.gz",
      "driverGenePanel": "$SAGE_DATA_ROOT/DriverGenePanel.hg38.tsv",
      "highConfBed": "$SAGE_DATA_ROOT/highConfidenceBed"
    }
  }

  call extractName as extractTumorName {
    input:
    refFasta = resources [ genomeVersion ].refFasta,
    refFai = resources [ genomeVersion ].refFai,
    modules = resources [ genomeVersion ].gatkModules,
    inputBam = tumour_bam,
    inputBai = tumour_bai
  }

  call extractName as extractNormalName {
    input:
    refFasta = resources [ genomeVersion ].refFasta,
    refFai = resources [ genomeVersion ].refFai,
    modules = resources [ genomeVersion ].gatkModules,
    inputBam = normal_bam,
    inputBai = normal_bai
  }

  call extractName as extractPlasmaName {
    input:
      refFasta = resources[genomeVersion].refFasta,
      refFai = resources[genomeVersion].refFai,
      modules = resources[genomeVersion].gatkModules,
      inputBam = plasma_bam,
      inputBai = plasma_bai
  }

  scatter (chr in chromosomes) {
    call splitPonByChromosome {
      input:
        PON = resources[genomeVersion].PON,
        chromosome = chr,
        modules = "bcftools/1.9 sage-data/1.0"
    }
  }

  
  scatter (idx in range(length(chromosomes))) {
    String chr_label = sub(chromosomes[idx], ",", "_")
    call amber as amberPrimaryChr {
      input:
        tumour_bam = tumour_bam,
        tumour_bai = tumour_bai,
        normal_bam = normal_bam,
        normal_bai = normal_bai,
        normal_name = extractNormalName.input_name,
        tumour_name = extractTumorName.input_name,  
        output_prefix = extractTumorName.input_name + "." + chr_label,  
        PON = splitPonByChromosome.chr_pon[idx],
        genomeVersion = genomeVersion,
        modules = resources[genomeVersion].hmfModules,
        threads = 2,
        memory = 8,
        timeout = 30
    }
}

  call mergeAmberChromosomes as amberPrimary {
    input:
      chr_zips = amberPrimaryChr.output_directory,
      tumour_name = extractTumorName.input_name,
      modules = "bcftools/1.9"
  }

  call cobalt as cobaltPrimary {
    input:
      tumour_bam = tumour_bam,
      tumour_bai = tumour_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = extractNormalName.input_name,
      tumour_name = extractTumorName.input_name,
      modules = resources[genomeVersion].hmfModules,
      gcProfile = resources[genomeVersion].gcProfile
  }

  call purple {
    input:
      normal_name = extractNormalName.input_name,
      tumour_name = extractTumorName.input_name,
      amber_directory = amberPrimary.output_directory,
      cobalt_directory = cobaltPrimary.output_directory,
      somatic_vcf = sage_primary_vcf,
      genomeVersion = genomeVersion,
      modules = resources[genomeVersion].hmfModules,
      gcProfile = resources[genomeVersion].gcProfile,
      ensemblDir = resources[genomeVersion].ensemblDir,
      refFasta = resources[genomeVersion].refFasta
  }

  scatter (idx in range(length(chromosomes))) {
    String chr_label_plasma = sub(chromosomes[idx], ",", "_")
    call amber as amberPlasmaChr {
      input:
        tumour_bam = plasma_bam,
        tumour_bai = plasma_bai,
        normal_bam = normal_bam,
        normal_bai = normal_bai,
        normal_name = extractNormalName.input_name,
        tumour_name = extractTumorName.input_name,
        output_prefix = extractTumorName.input_name + "." + chr_label_plasma,
        PON = splitPonByChromosome.chr_pon[idx],
        genomeVersion = genomeVersion,
        modules = resources[genomeVersion].hmfModules,
        threads = 2,  
        memory = 8,   
        timeout = 30  
    }
  }

  call mergeAmberChromosomes as amberPlasma {
    input:
      chr_zips = amberPlasmaChr.output_directory,  
      tumour_name = extractPlasmaName.input_name,  
      modules = "bcftools/1.9"
  }

  call cobalt as cobaltPlasma {
    input:
      tumour_bam = plasma_bam,
      tumour_bai = plasma_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = extractNormalName.input_name,
      tumour_name = extractPlasmaName.input_name,
      modules = resources[genomeVersion].hmfModules,
      gcProfile = resources[genomeVersion].gcProfile
  }

  call mergeDirs as mergeAmber {
    input:
      primary_zip = amberPrimary.output_directory,
      plasma_zip = amberPlasma.output_directory,
      output_name = "merged_amber"
  }

  call mergeDirs as mergeCobalt {
    input:
      primary_zip = cobaltPrimary.output_directory,
      plasma_zip = cobaltPlasma.output_directory,
      output_name = "merged_cobalt"
  }

  scatter (chr in chromosomes) {
    call sage {
      input:
        chromosome = chr,
        reference_name = extractPlasmaName.input_name,
        reference_bam = plasma_bam,
        reference_bai = plasma_bai,
        refFasta = resources[genomeVersion].refFasta,
        input_vcf = sage_primary_vcf,
        input_vcf_index = sage_primary_vcf_index,
        modules = "sage/3.4.4 hg38/p12",
        ensemblDir = resources[genomeVersion].ensemblDir,
        hotspots = resources[genomeVersion].hotspots,
        driverGenePanel = resources[genomeVersion].driverGenePanel,
        highConfBed = resources[genomeVersion].highConfBed
    }
  }

  # Merge plasma chromosome VCFs
  call mergeVcfs as mergePlasmaVcfs {
    input:
      vcfs = sage.output_vcf,
      vcf_indices = sage.output_vcf_index,
      sample_name = extractPlasmaName.input_name,
      modules = "bcftools/1.9"
  }

  # Merge plasma BQR directories
  call mergeBqrDirs as mergePlasmaBqr {
    input:
      bqr_zips = sage.output_bqr_directory,
      sample_name = extractPlasmaName.input_name
  }

  call annotatePlasmaVcfWithPurple {
    input:
      purple_vcf = purple.purple_somatic_vcf,
      purple_vcf_index = purple.purple_somatic_vcf_index,
      plasma_sage_vcf = mergePlasmaVcfs.merged_vcf,
      plasma_sage_vcf_index = mergePlasmaVcfs.merged_vcf_index
  }


  call runWisp {
    input:
      donor = donor,  
      tumour_name = extractTumorName.input_name,
      plasma_name = extractPlasmaName.input_name,
      purple_dir = purple.purple_directory,
      amber_dir = mergeAmber.output_zip,
      cobalt_dir = mergeCobalt.output_zip,
      somatic_vcf = annotatePlasmaVcfWithPurple.annotated_vcf,
      somatic_vcf_index = annotatePlasmaVcfWithPurple.annotated_vcf_index,
      bqr_dir = mergePlasmaBqr.merged_bqr_zip,
      refFasta = resources[genomeVersion].refFasta,
      genomeVersion = genomeVersion,
      modules = resources[genomeVersion].wispModules
  }

  meta {
    author: "Gavin Peng"
    email: "gpeng@oicr.on.ca"
    description: "WISP tumor fraction estimation with chromosome-split SAGE plasma append"
    dependencies: [
        {
            name: "PURPLE",
            url: "https://github.com/hartwigmedical/hmftools/blob/master/purple/README.md"
        },
        {
            name: "AMBER",
            url: "https://github.com/hartwigmedical/hmftools/blob/master/amber/README.md"
        },
        {
            name: "COBALT",
            url: "https://github.com/hartwigmedical/hmftools/blob/master/cobalt/README.md"
        },
        {
            name: "WISP",
            url: "https://github.com/hartwigmedical/hmftools/tree/master/wisp"
        },
        {
            name: "SAGE",
            url: "https://github.com/hartwigmedical/hmftools/tree/master/sage"
        }
    ]
    output_meta: {
      wisp_directory: "Zipped WISP output directory",
      wisp_summary: "WISP tumor fraction summary",
      wisp_snv_results: "Per-variant SNV results",
      sage_plasma_vcf: "SAGE VCF with plasma sample appended",
      sage_plasma_vcf_index: "Index for SAGE plasma VCF",
      sage_plasma_bqr: "Zipped SAGE BQR results directory for plasma"
    }
  }

  output {
    File wisp_directory = runWisp.wisp_directory
    File wisp_summary = runWisp.wisp_summary
    File? wisp_snv_results = runWisp.wisp_snv_results
    File sage_plasma_vcf = mergePlasmaVcfs.merged_vcf
    File sage_plasma_vcf_index = mergePlasmaVcfs.merged_vcf_index
    File sage_plasma_bqr = mergePlasmaBqr.merged_bqr_zip
  }
}

task extractName {
  input {
    String modules
    String refFasta
    String refFai
    File inputBam
    File inputBai
    Int memory = 4
    Int heapRam = 1
    Int timeout = 4
  }

  parameter_meta {
    inputBam: "input .bam file"
    inputBai: "input .bai file"
    refFasta: "Reference FASTA file"
    refFai: "Reference fai index"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    heapRam: "Heap RAM allocation for GATK (GB)"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    if [ -f "~{inputBam}" ]; then
      gatk --java-options "-Xmx~{heapRam}g" GetSampleName -R ~{refFasta} -I ~{inputBam} -O input_name.txt -encode
    fi

    cat input_name.txt
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  meta {
    output_meta: {
      input_name: "name of the input"
    }
  }

  output {
    String input_name = read_string(stdout()) 
  }
}

task splitPonByChromosome {
  input {
    String PON 
    String chromosome
    String modules
    Int memory = 4
    Int timeout = 10
  }
  parameter_meta {
    PON: "Panel of Normal (PON) file, generated for AMBER"
    chromosome: "Chromosome to extract from PON"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail
    
    # Extract just this chromosome from PON
    bcftools view -r ~{chromosome} ~{PON} -O z -o pon.~{chromosome}.vcf.gz
    tabix -p vcf pon.~{chromosome}.vcf.gz
  >>>

  runtime {
    memory: "~{memory} GB"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    File chr_pon = "pon.~{chromosome}.vcf.gz"
    File chr_pon_index = "pon.~{chromosome}.vcf.gz.tbi"
  }
}



task amber {
  input {
    String tumour_name
    File tumour_bam
    File tumour_bai
    String normal_name
    File normal_bam
    File normal_bai
    String? output_prefix
    String PON
    String genomeVersion
    Int min_mapping_quality = 30
    Int min_base_quality = 25
    String modules
    Int threads = 8
    Int memory = 32
    Int heapRam = 32
    Int timeout = 100
    String additionalParameters = ""
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    tumour_bam: "Tumour bam"
    tumour_bai: "Matching bai for Tumour bam"
    normal_name: "Name for Normal sample"
    normal_bam: "Normal bam"
    normal_bai: "Matching bai for Normal bam"
    output_prefix: "prefix for output"
    PON: "Panel of Normal (PON) file, generated for AMBER"
    genomeVersion: "genome version (37 or 38)"
    min_mapping_quality: "Minimum mapping quality for an alignment to be used"
    min_base_quality: "Minimum quality for a base to be considered"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    heapRam: "Heap RAM allocation for AMBER (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
    additionalParameters: "Additional parameters to pass to AMBER"
  }

  String file_prefix = select_first([output_prefix, tumour_name])

  command <<<
    set -euo pipefail

    mkdir ~{file_prefix}.amber

    java -Xmx~{heapRam}G -cp $HMFTOOLS_ROOT/amber.jar com.hartwig.hmftools.amber.AmberApplication \
      -reference ~{normal_name} -reference_bam ~{normal_bam} \
      -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
      -output_dir ~{file_prefix}.amber/ \
      -loci ~{PON} \
      -ref_genome_version ~{genomeVersion} \
      -min_mapping_quality ~{min_mapping_quality} \
      -min_base_quality ~{min_base_quality} \
      ~{additionalParameters}

    zip -r ~{file_prefix}.amber.zip ~{file_prefix}.amber/

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File output_directory = "~{file_prefix}.amber.zip"
  }
  meta {
    output_meta: {
      output_directory: "Zipped AMBER results directory"
    }
  }
}

task mergeAmberChromosomes {
  input {
    Array[File] chr_zips
    String tumour_name
    Int memory = 8
    Int timeout = 10
    String modules
  }

  parameter_meta {
    chr_zips: "Array of chromosome-split AMBER result zips"
    tumour_name: "Name for Tumour sample"
    memory: "Memory allocated for this job (GB)"
    timeout: "Hours before task timeout"
    modules: "Required environment modules"
  }

  command <<<
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
  >>>

  runtime {
    memory: "~{memory} GB"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  output {
    File output_directory = "~{tumour_name}.amber.zip"
  }

  meta {
    output_meta: {
      output_directory: "Zipped merged AMBER results directory"
    }
  }
}

task cobalt {
  input {
    String tumour_name
    File tumour_bam
    File tumour_bai
    String normal_name
    File normal_bam
    File normal_bai
    String gcProfile
    String gamma = "300"
    Int min_mapping_quality = 30
    String modules
    Int threads = 8
    Int memory = 32
    Int heapRam = 8
    Int timeout = 100
    String additionalParameters = ""
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    tumour_bam: "Tumour bam"
    tumour_bai: "Matching bai for Tumour bam"
    normal_name: "Name for Normal sample"
    normal_bam: "Normal bam"
    normal_bai: "Matching bai for Normal bam"
    gcProfile: "GC profile, generated for COBALT"
    gamma: "gamma (penalty) value for segmenting"
    min_mapping_quality: "Minimum mapping quality for an alignment to be used"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    heapRam: "Heap RAM allocation for COBALT (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
    additionalParameters: "Additional parameters to pass to COBALT"
  }

  command <<<
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

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File output_directory = "~{tumour_name}.cobalt.zip"
  }

  meta {
    output_meta: {
      output_directory: "Zipped COBALT results directory"
    }
  }
}

task mergeDirs {
  input {
    File primary_zip
    File plasma_zip
    String output_name
  }
  command <<<
    mkdir -p ~{output_name}
    unzip -o ~{primary_zip} -d ~{output_name}/
    unzip -o ~{plasma_zip} -d ~{output_name}/
    # Flatten if nested
    find ~{output_name} -mindepth 2 -type f -exec mv {} ~{output_name}/ \;
    find ~{output_name} -mindepth 1 -type d -empty -delete
    zip -r ~{output_name}.zip ~{output_name}/
  >>>
  output {
    File output_zip = "~{output_name}.zip"
  }
}

task sage {
  input {
    String chromosome
    String reference_name
    File reference_bam
    File reference_bai
    String refFasta
    String ensemblDir
    String hotspots
    String driverGenePanel
    String highConfBed
    File input_vcf
    File input_vcf_index
    Int min_map_quality = 10
    Int hard_min_tumor_qual = 50
    Int hard_min_tumor_raw_alt_support = 2
    Float hard_min_tumor_vaf = 0.002
    String modules
    Int threads = 8
    Int memory = 40
    Int heapRam = 32
    Int timeout = 24
    String additionalParameters = ""
  }

  parameter_meta {
    chromosome: "Chromosome to process"
    reference_name: "Name for cfDNA sample"
    reference_bam: "cfDNA bam"
    reference_bai: "Matching bai for cfDNA bam"
    refFasta: "Reference genome fasta"
    ensemblDir: "Ensembl data directory"
    hotspots: "Known hotspots VCF"
    driverGenePanel: "Driver gene panel TSV"
    highConfBed: "High confidence regions BED"
    input_vcf: "Input VCF from primary SAGE run"
    input_vcf_index: "Index for input VCF"
    min_map_quality: "Minimum mapping quality"
    hard_min_tumor_qual: "Hard minimum tumor quality"
    hard_min_tumor_raw_alt_support: "Minimum raw alt support"
    hard_min_tumor_vaf: "Minimum tumor VAF"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    heapRam: "Heap RAM allocation for SAGE (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
    additionalParameters: "Additional parameters to pass to SAGE"
  }

  command <<<
    set -euo pipefail

    mkdir -p ~{reference_name}.sage.bqr
    echo $SAGE_ROOT
    ls -la ${SAGE_ROOT}/

    java -Xmx~{heapRam}G -cp $SAGE_ROOT/sage.jar com.hartwig.hmftools.sage.append.SageAppendApplication \
      -reference ~{reference_name} \
      -reference_bam ~{reference_bam} \
      -ref_genome_version 38 \
      -ref_genome ~{refFasta} \
      -input_vcf ~{input_vcf} \
      -specific_chr ~{chromosome} \
      -output_vcf ~{reference_name}.~{chromosome}.sage.vcf.gz \
      -threads ~{threads} \
      -min_map_quality ~{min_map_quality} \
      -hard_min_tumor_qual ~{hard_min_tumor_qual} \
      -hard_min_tumor_raw_alt_support ~{hard_min_tumor_raw_alt_support} \
      -hard_min_tumor_vaf ~{hard_min_tumor_vaf} \
      ~{additionalParameters}

    mv *.sage.bqr.tsv ~{reference_name}.sage.bqr/ 2>/dev/null || true
    zip -r ~{reference_name}.~{chromosome}.sage.bqr.zip ~{reference_name}.sage.bqr/

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File output_bqr_directory = "~{reference_name}.~{chromosome}.sage.bqr.zip"
    File output_vcf = "~{reference_name}.~{chromosome}.sage.vcf.gz"
    File output_vcf_index = "~{reference_name}.~{chromosome}.sage.vcf.gz.tbi"
  }

  meta {
    output_meta: {
      output_bqr_directory: "Zipped SAGE BQR results directory",
      output_vcf: "SAGE output VCF",
      output_vcf_index: "SAGE output VCF index"
    }
  }
}

task mergeVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_indices
    String sample_name
    String modules
    Int memory = 8
    Int timeout = 4
  }
  parameter_meta {
    vcfs: "Array of VCF files to merge"
    vcf_indices: "Array of VCF index files"
    sample_name: "Sample name for output files"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail
    
    # Create file list for bcftools concat
    for vcf in ~{sep=' ' vcfs}; do
      echo "$vcf" >> vcf_list.txt
    done
    
    # Sort by chromosome order
    sort -V vcf_list.txt > vcf_list_sorted.txt
    
    # Concatenate VCFs in order
    bcftools concat \
      --file-list vcf_list_sorted.txt \
      --output-type z \
      --output ~{sample_name}.sage.vcf.gz
    
    # Index the merged VCF
    tabix -p vcf ~{sample_name}.sage.vcf.gz
  >>>

  runtime {
    memory: "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File merged_vcf = "~{sample_name}.sage.vcf.gz"
    File merged_vcf_index = "~{sample_name}.sage.vcf.gz.tbi"
  }
}

task mergeBqrDirs {
  input {
    Array[File] bqr_zips
    String sample_name
    Int memory = 4
    Int timeout = 2
  }
  parameter_meta {
    bqr_zips: "Array of zipped BQR directories to merge"
    sample_name: "Sample name for output files"
    memory: "Memory allocated for this job (GB)"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail
    
    mkdir -p ~{sample_name}.sage.bqr
    
    # Unzip all BQR files into the same directory
    for bqr_zip in ~{sep=' ' bqr_zips}; do
      unzip -o "$bqr_zip" -d temp_bqr/
    done
    
    # Move all .tsv files to final directory
    find temp_bqr -name "*.sage.bqr.tsv" -exec mv {} ~{sample_name}.sage.bqr/ \;
    
    # Zip the merged directory
    zip -r ~{sample_name}.sage.bqr.zip ~{sample_name}.sage.bqr/
  >>>

  runtime {
    memory: "~{memory} GB"
    timeout: "~{timeout}"
  }

  output {
    File merged_bqr_zip = "~{sample_name}.sage.bqr.zip"
  }
}

task runWisp {
  input {
    String donor
    String tumour_name
    String plasma_name
    File? purple_dir
    File amber_dir  
    File cobalt_dir
    File? somatic_vcf
    File? somatic_vcf_index
    File bqr_dir
    String refFasta
    String genomeVersion
    String modules
    String additionalParameters = ""
    Int threads = 4
    Int heapRam = 16
    Int memory = 16
    Int timeout = 50
  }

  parameter_meta {
    donor: "Patient identifier"
    tumour_name: "Primary tumor sample name"
    plasma_name: "Plasma/cfDNA sample name"
    purple_dir: "Zipped PURPLE directory from primary tumor"
    amber_dir: "Zipped merged AMBER directory"
    cobalt_dir: "Zipped merged COBALT directory"
    somatic_vcf: "SAGE VCF with plasma sample appended"
    somatic_vcf_index: "Index for somatic VCF"
    bqr_dir: "Zipped SAGE BQR directory from plasma"
    refFasta: "Reference genome fasta"
    genomeVersion: "Genome version (37 or 38)"
    additionalParameters: "Additional parameters to pass to WISP"
    modules: "Required environment modules"
    heapRam: "Heap RAM allocation for WISP (GB)"
    memory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    # Unzip input directories
    unzip ~{purple_dir}
    unzip ~{amber_dir}
    unzip ~{cobalt_dir}
    unzip ~{bqr_dir}

    # Create output directory
    mkdir -p ~{plasma_name}.wisp

    # Run WISP
    java -Xmx~{heapRam}G -jar $WISP_ROOT/wisp.jar \
      -patient_id ~{donor} \
      -tumor_id ~{tumour_name} \
      -samples ~{plasma_name} \
      -purple_dir ~{tumour_name}.solPrimary.purple/ \
      -amber_dir merged_amber/ \
      -cobalt_dir merged_cobalt/ \
      -somatic_vcf ~{somatic_vcf} \
      -bqr_dir ~{plasma_name}.sage.bqr/ \
      -ref_genome ~{refFasta} \
      -output_dir ~{plasma_name}.wisp/ \
      -threads ~{threads} \
      ~{additionalParameters}

    # Zip output
    zip -r ~{plasma_name}.wisp.zip ~{plasma_name}.wisp/

  >>>

  runtime {
    cpu: "~{threads}"
    memory: "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File wisp_directory = "~{plasma_name}.wisp.zip"
    File wisp_summary = "~{plasma_name}.wisp/~{donor}_~{plasma_name}.wisp.summary.tsv"
    File? wisp_snv_results = "~{plasma_name}.wisp/~{donor}_~{plasma_name}.wisp.snv.tsv"  # Optional - may not exist
  }

  meta {
    output_meta: {
      wisp_directory: "Zipped WISP output directory",
      wisp_summary: "WISP tumor fraction summary",
      wisp_snv_results: "Per-variant SNV results"
    }
  }
}

task purple {
  input {
    String normal_name
    String tumour_name
    String solution_name = "Primary"
    String outfilePrefix = tumour_name + ".sol" + solution_name
    File amber_directory
    File cobalt_directory
    File somatic_vcf
    String ensemblDir
    String refFasta
    String genomeVersion
    String gcProfile
    Int min_diploid_tumor_ratio_count = 60
    String? min_ploidy
    String? max_ploidy
    String? min_purity
    String? max_purity
    String? ploidy_penalty_factor
    String? ploidy_penalty_standard_deviation
    String modules
    Int threads = 8
    Int memory = 32
    Int heapRam = 8
    Int timeout = 100
    String additionalParameters = ""
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    normal_name: "Name for Normal sample"
    solution_name: "Name of solution"
    outfilePrefix: "Prefix of output file"
    amber_directory: "zipped output from AMBER"
    cobalt_directory: "zipped output from COBALT"
    somatic_vcf: "somatic variants vcf"
    ensemblDir: "Directory of Ensembl data for PURPLE"
    refFasta: "fasta of reference genome"
    gcProfile: "GC profile, generated for COBALT"
    min_diploid_tumor_ratio_count: "smooth over contiguous segments"
    min_ploidy: "minimum ploidy to consider"
    max_ploidy: "maximum ploidy to consider"
    min_purity: "minimum purity to consider"
    max_purity: "maximum purity to consider"
    ploidy_penalty_factor: "ploidy penalty factor"
    ploidy_penalty_standard_deviation: "ploidy penalty standard deviation"
    genomeVersion: "genome version for AMBER, default set to V38"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    heapRam: "Heap RAM allocation for PURPLE (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
    additionalParameters: "Additional parameters to pass to PURPLE"
  }

  command <<<
    set -euo pipefail

    unzip ~{amber_directory}
    unzip ~{cobalt_directory}
    mkdir ~{outfilePrefix}.purple

    java -Xmx~{heapRam}G -jar $HMFTOOLS_ROOT/purple.jar \
      -ref_genome_version ~{genomeVersion} \
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

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File? purple_directory = "~{outfilePrefix}.purple.zip"
    File purple_qc = "~{outfilePrefix}.purple/~{tumour_name}.purple.qc"
    File purple_purity = "~{outfilePrefix}.purple/~{tumour_name}.purple.purity.tsv"
    File purple_purity_range = "~{outfilePrefix}.purple/~{tumour_name}.purple.purity.range.tsv"
    File purple_segments = "~{outfilePrefix}.purple/~{tumour_name}.purple.segment.tsv"
    File purple_cnv = "~{outfilePrefix}.purple/~{tumour_name}.purple.cnv.somatic.tsv"
    File purple_cnv_gene = "~{outfilePrefix}.purple/~{tumour_name}.purple.cnv.gene.tsv"
    File purple_somatic_vcf = "~{outfilePrefix}.purple/~{tumour_name}.purple.somatic.vcf.gz"
    File purple_somatic_vcf_index = "~{outfilePrefix}.purple/~{tumour_name}.purple.somatic.vcf.gz.tbi"
  }

  meta {
    output_meta: {
      purple_directory: "Zipped Output PURPLE directory",
      purple_qc: "QC results from PURPLE",
      purple_purity: "Purity estimate from PURPLE",
      purple_purity_range: "Range of Purity estimate from PURPLE",
      purple_segments: "Segments estimated by PURPLE",
      purple_cnv: "Somatic copy number variants from PURPLE",
      purple_cnv_gene: "Gene-level copy number variants from PURPLE"
    }
  }
}


task annotatePlasmaVcfWithPurple {
  input {
    File purple_vcf
    File purple_vcf_index
    File plasma_sage_vcf
    File plasma_sage_vcf_index
    String modules = "bcftools/1.9"
  }
  parameter_meta {
    purple_vcf: "VCF from PURPLE containing copy number and purity annotations"
    purple_vcf_index: "Index for PURPLE VCF"
    plasma_sage_vcf: "SAGE VCF with plasma sample appended"
    plasma_sage_vcf_index: "Index for plasma SAGE VCF"
    modules: "Required environment modules"
  }

  command <<<
    # Extract Purple INFO annotations and merge into plasma VCF
    bcftools annotate \
      -a ~{purple_vcf} \
      -c INFO/SUBCL,INFO/PURPLE_VCN,INFO/PURPLE_AF,INFO/PURPLE_CN \
      ~{plasma_sage_vcf} \
      -Oz -o annotated_plasma.vcf.gz
    
    tabix -p vcf annotated_plasma.vcf.gz
  >>>

   runtime {
    modules: "~{modules}"
  }
  
  output {
    File annotated_vcf = "annotated_plasma.vcf.gz"
    File annotated_vcf_index = "annotated_plasma.vcf.gz.tbi"
  }
}