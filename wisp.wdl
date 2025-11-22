version 1.0

struct GenomeResources {
    String PON
    String modules
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
  }

  parameter_meta {
    tumour_bam: "Input tumor file (bam) of primary sample"
    tumour_bai: "Input tumor file index (bai) of primary sample"
    normal_bam: "Input normal file (bam) of primary sample"
    normal_bai: "Input normal file index (bai) of primary sample"
    plasma_bam: "Input of bam file of plasma from same donor as primary sample"
    plasma_bai: "Input of bam index file of plasma from same donor"
    donor: "The donor of all the samples"
    genomeVersion: "Genome Version, only 38 supported"
  }

Map[String,GenomeResources] resources = {
  "38": {
    "modules": "sage wisp hmftools/1.1 hg38/p12 hmftools-data/53138",
    "gatkModules": "hg38-gridss-index/1.0 gatk/4.1.6.0",
    "refFasta": "$HG38_ROOT/hg38_random.fa",
    "refFai": "$HG38_GRIDSS_INDEX_ROOT/hg38_random.fa.fai",
    "PON" : "$HMFTOOLS_DATA_ROOT/copy_number/GermlineHetPon.38.vcf.gz",
    "ensemblDir": "$HMFTOOLS_DATA_ROOT/ensembl_data",
    "gcProfile": "$HMFTOOLS_DATA_ROOT/copy_number/GC_profile.1000bp.38.cnp",
    "pon_sgl_file": "$HMFTOOLS_DATA_ROOT/sv/sgl_pon.38.bed.gz",
    "pon_sv_file": "$HMFTOOLS_DATA_ROOT/sv/sv_pon.38.bedpe.gz",
    "known_hotspot_file": "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe",
    "repeat_mask_file": "$HMFTOOLS_DATA_ROOT/sv/repeat_mask_data.38.fa.gz",
    "knownfusion": "$HMFTOOLS_DATA_ROOT/sv/known_fusions.38.bedpe",
    "hotspots": "/.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/dev/sage/KnownHotspots.hg38.fixed.vcf.gz",
    "driverGenePanel": "/.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/dev/sage/DriverGenePanel.hg38.tsv",
    "highConfBed": "/.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/dev/sage/NA12878_GIAB_highconf_IllFB-IllGATKHC-CG-Ion-Solid_ALLCHROM_v3.2.2_highconf.bed"
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
    refFasta = resources [ genomeVersion ].refFasta,
    refFai = resources [ genomeVersion ].refFai,
    modules = resources [ genomeVersion ].gatkModules,
    inputBam = plasma_bam,
    inputBai = plasma_bai
  }

  call amber as amberPrimary {
    input:
      tumour_bam = tumour_bam,
      tumour_bai = tumour_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = extractNormalName.input_name,
      tumour_name = extractTumorName.input_name,
      genomeVersion = genomeVersion,
      modules = resources [ genomeVersion ].modules,
      PON = resources [ genomeVersion ].PON
  }

  call cobalt as cobaltPrimary {
    input:
      tumour_bam = tumour_bam,
      tumour_bai = tumour_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = extractNormalName.input_name,
      tumour_name = extractTumorName.input_name,
      modules = resources [ genomeVersion ].modules,
      gcProfile = resources [ genomeVersion ].gcProfile
  }

  call sage as sagePrimary {
    input:
      tumour_bam = tumour_bam,
      tumour_bai = tumour_bai,
      reference_bam = normal_bam,
      reference_bai = normal_bai,
      reference_name = extractNormalName.input_name,
      tumour_name = extractTumorName.input_name,
      append_mode = false,
      modules = resources [ genomeVersion ].modules,
      refFasta = resources [ genomeVersion ].refFasta,
      ensemblDir = resources [ genomeVersion ].ensemblDir,
      hotspots = resources [ genomeVersion ].hotspots,
      driverGenePanel = resources [ genomeVersion ].driverGenePanel,
      highConfBed = resources [ genomeVersion ].highConfBed
  }

  call purple {
    input:
      normal_name = extractNormalName.input_name,
      tumour_name = extractTumorName.input_name,
      amber_directory = amberPrimary.output_directory,
      cobalt_directory = cobaltPrimary.output_directory,
      somatic_vcf = sagePrimary.output_vcf,  
      genomeVersion = genomeVersion,
      modules = resources [ genomeVersion ].modules,
      gcProfile = resources [ genomeVersion ].gcProfile,
      ensemblDir = resources [ genomeVersion ].ensemblDir,
      refFasta = resources [ genomeVersion ].refFasta
  }

  call amber as amberPlasma {
    input:
      tumour_bam = plasma_bam,
      tumour_bai = plasma_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = extractNormalName.input_name,
      tumour_name = extractPlasmaName.input_name,
      genomeVersion = genomeVersion,
      modules = resources [ genomeVersion ].modules,
      PON = resources [ genomeVersion ].PON
  }

  call cobalt as cobaltPlasma {
    input:
      tumour_bam = plasma_bam,
      tumour_bai = plasma_bai,
      normal_bam = normal_bam,
      normal_bai = normal_bai,
      normal_name = extractNormalName.input_name,
      tumour_name = extractPlasmaName.input_name,
      modules = resources [ genomeVersion ].modules,
      gcProfile = resources [ genomeVersion ].gcProfile
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

  call sage as sagePlasma {
    input:
      reference_bam = plasma_bam,
      reference_bai = plasma_bai,
      reference_name = extractPlasmaName.input_name,
      refFasta = resources[genomeVersion].refFasta,
      append_mode = true,
      input_vcf = sagePrimary.output_vcf,  #  run in append mode, use VCF from primary
      modules = resources[genomeVersion].modules,
      ensemblDir = resources [ genomeVersion ].ensemblDir,
      hotspots = resources [ genomeVersion ].hotspots,
      driverGenePanel = resources [ genomeVersion ].driverGenePanel,
      highConfBed = resources [ genomeVersion ].highConfBed
  }

  call runWisp {
    input:
      donor = donor,  
      tumour_name = extractTumorName.input_name,
      plasma_name = extractPlasmaName.input_name,
      purple_dir = purple.purple_directory,
      amber_dir = mergeAmber.output_zip,
      cobalt_dir = mergeCobalt.output_zip,
      somatic_vcf = sagePlasma.output_vcf,
      somatic_vcf_index = sagePlasma.output_vcf_index,
      bqr_dir = sagePlasma.output_bqr_directory,
      refFasta = resources[genomeVersion].refFasta,
      genomeVersion = genomeVersion,
      modules = resources[genomeVersion].modules
  }

  meta {
    author: "Gavin Peng"
    email: "gpeng@oicr.on.ca"
    description: "The WISP (Whole-genome Inference of Somatic Plasma) workflow estimates circulating tumor DNA (ctDNA) fraction in plasma samples by leveraging somatic variants and copy number profiles derived from matched primary tumor sequencing data."
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
            name: "bcftools",
            url: "https://www.htslib.org/doc/1.9/bcftools.html"
        },
        {
            name: "wisp",
            url: "https://github.com/hartwigmedical/hmftools/tree/master/wisp"
        },
        {
            name: "sage",
            url: "https://github.com/hartwigmedical/hmftools/tree/master/sage"
        }
    ]
    output_meta: {
      wisp_directory: "Zipped WISP output directory",
      wisp_summary: "WISP tumor fraction summary",
      wisp_snv_results: "Per-variant SNV results"
    }
  }
  output {
    File wisp_directory = runWisp.wisp_directory
    File wisp_summary = runWisp.wisp_summary
    File wisp_snv_results = runWisp.wisp_snv_results
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
    Int timeout = 4
  }

  parameter_meta {
    inputBam: "input .bam file"
    inputBai: "input .bai file"
    refFasta: "Reference FASTA file"
    refFai: "Reference fai index"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    timeout: "Hours before task timeout"
  }

  command <<<
    set -euo pipefail

    if [ -f "~{inputBam}" ]; then
      gatk --java-options "-Xmx1g" GetSampleName -R ~{refFasta} -I ~{inputBam} -O input_name.txt -encode
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

task amber {
  input {
    String tumour_name
    File tumour_bam
    File tumour_bai
    String normal_name
    File normal_bam
    File normal_bai
    String amberScript = "java -Xmx32G -cp $HMFTOOLS_ROOT/amber.jar com.hartwig.hmftools.amber.AmberApplication"
    String PON
    String genomeVersion
    Int min_mapping_quality = 30
    Int min_base_quality = 25
    String modules
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    tumour_bam: "Tumour bam"
    tumour_bai: "Matching bai for Tumour bam"
    normal_name: "Name for Normal sample"
    normal_bam: "Normal bam"
    normal_bai: "Matching bai for Normal bam"
    amberScript: "location of AMBER script"
    PON: "Panel of Normal (PON) file, generated for AMBER"
    genomeVersion: "genome version (37 or 38)"
    min_mapping_quality: "Minimum mapping quality for an alignment to be used"
    min_base_quality: "Minimum quality for a base to be considered"
		modules: "Required environment modules"
		memory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

  command <<<
    set -euo pipefail

    mkdir ~{tumour_name}.amber  

    ~{amberScript} \
      -reference ~{normal_name} -reference_bam ~{normal_bam} \
      -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
      -output_dir ~{tumour_name}.amber/ \
      -loci ~{PON} \
      -ref_genome_version ~{genomeVersion} \
      -min_mapping_quality ~{min_mapping_quality} \
      -min_base_quality ~{min_base_quality} 

    zip -r ~{tumour_name}.amber.zip ~{tumour_name}.amber/

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File output_directory = "~{tumour_name}.amber.zip"
  }
  meta {
		output_meta: {
			output_directory: "Zipped AMBER results directory"
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
    String colbaltScript = "java -Xmx8G -cp $HMFTOOLS_ROOT/cobalt.jar com.hartwig.hmftools.cobalt.CobaltApplication"
    String gcProfile
    String gamma = 300
    Int min_mapping_quality = 30
    String modules
    Int threads = 8
    Int memory = 32
    Int timeout = 100
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    tumour_bam: "Tumour bam"
    tumour_bai: "Matching bai for Tumour bam"
    normal_name: "Name for Normal sample"
    normal_bam: "Normal bam"
    normal_bai: "Matching bai for Normal bam"
    colbaltScript: "location of COBALT script"
    gcProfile: "GC profile, generated for COBALT"
    gamma: "gamma (penalty) value for segmenting"
    min_mapping_quality: "Minimum mapping quality for an alignment to be used"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
	}

  command <<<
    set -euo pipefail

    mkdir ~{tumour_name}.cobalt 

    ~{colbaltScript} \
      -reference ~{normal_name} -reference_bam ~{normal_bam} \
      -tumor ~{tumour_name} -tumor_bam ~{tumour_bam} \
      -output_dir ~{tumour_name}.cobalt/ \
      -gc_profile ~{gcProfile} \
      -pcf_gamma ~{gamma} \
      -min_quality ~{min_mapping_quality}

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
    String? tumour_name
    File? tumour_bam
    File? tumour_bai
    String reference_name
    File reference_bam
    File reference_bai
    String refFasta
    String ensemblDir
    String hotspots
    String driverGenePanel
    String highConfBed
    Boolean append_mode = false
    File? input_vcf  # For append mode
    Int min_map_quality = 10 
    Int hard_min_tumor_qual = 50
    Int hard_min_tumor_raw_alt_support = 2 
    Float hard_min_tumor_vaf = 0.002 
    String modules
    Int threads = 8
    Int memory = 40
    Int timeout = 100
  }

  parameter_meta {
    tumour_name: "Name for Tumour sample"
    tumour_bam: "Tumour bam"
    tumour_bai: "Matching bai for Tumour bam"
    reference_name: "Name for reference sample, in append mode this is cfDNA sample"
    reference_bam: "reference bam"
    reference_bai: "Matching bai for reference bam"
    refFasta: "Reference genome fasta"
    ensemblDir: "Ensembl data directory"
    hotspots: "Known hotspots VCF"
    driverGenePanel: "Driver gene panel TSV"
    highConfBed: "High confidence regions BED"
    append_mode: "whether run sage in append mode"
    input_vcf: "Input VCF for append mode (optional)"
    min_map_quality: "Minimum mapping quality"
    hard_min_tumor_qual: "Hard minimum tumor quality"
    hard_min_tumor_raw_alt_support: "Minimum raw alt support" 
    hard_min_tumor_vaf: "Minimum tumor VAF"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
  }
  

  command <<<
    set -euo pipefail
    mkdir -p ~{tumour_name}.sage.bqr 

    if ~{append_mode}; then
        sageClass="com.hartwig.hmftools.sage.append.SageAppendApplication"
    else
        sageClass="com.hartwig.hmftools.sage.SageApplication"
    fi
    sage_jar="/.mounts/labs/gsiprojects/gsi/gsiusers/gpeng/dev/sage/sage_v3.4.4.jar"
    
  
    java -Xmx32G -cp ${sage_jar} ${sageClass}   \
      ~{if append_mode then "" else "-tumor " + tumour_name} \
      ~{if append_mode then "" else "-tumor_bam " + tumour_bam} \
      -reference ~{reference_name} \
      -reference_bam ~{reference_bam} \
      -ref_genome_version 38 \
      -ref_genome ~{refFasta} \
      ~{"-input_vcf " + input_vcf} \
      ~{if !defined(input_vcf) then "-ensembl_data_dir " + ensemblDir else ""} \
      ~{if !defined(input_vcf) then "-high_confidence_bed " + highConfBed else ""} \
      ~{if !defined(input_vcf) then "-hotspots " + hotspots else ""} \
      -output_vcf ~{tumour_name}.sage.vcf.gz \
      -threads ~{threads} \
      -min_map_quality ~{min_map_quality} \
      -hard_min_tumor_qual ~{hard_min_tumor_qual} \
      -hard_min_tumor_raw_alt_support ~{hard_min_tumor_raw_alt_support} \
      -hard_min_tumor_vaf ~{hard_min_tumor_vaf}

    mv *.sage.bqr.tsv ~{tumour_name}.sage.bqr/ 2>/dev/null || true
    zip -r ~{tumour_name}.sage.bqr.zip ~{tumour_name}.sage.bqr/

  >>>

  runtime {
    cpu: "~{threads}"
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    File output_bqr_directory = "~{tumour_name}.sage.bqr.zip"
    File output_vcf = "~{tumour_name}.sage.vcf.gz"
    File output_vcf_index = "~{tumour_name}.sage.vcf.gz.tbi"
  }

  meta {
    output_meta: {
      output_bqr_directory: "Zipped SAGE BQR results directory",
      output_vcf: "SAGE output VCF",
      output_vcf_index: "SAGE output VCF index"
    }
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
    Int threads = 4
    Int memory = 16
    Int timeout = 50
  }

  parameter_meta {
    donor: "Patient identifier"
    tumour_name: "Primary tumor sample name"
    plasma_name: "Plasma/cfDNA sample name"
    purple_dir: "Zipped PURPLE directory from primary tumor"
    amber_dir: "Zipped AMBER directory from plasma"
    cobalt_dir: "Zipped COBALT directory from plasma"
    somatic_vcf: "SAGE VCF with plasma sample appended"
    somatic_vcf_index: "Index for somatic VCF"
    bqr_dir: "Zipped SAGE BQR directory from plasma"
    refFasta: "Reference genome fasta"
    genomeVersion: "Genome version (37 or 38)"
    modules: "Required environment modules"
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
    java -Xmx16G -jar $WISP_ROOT/wisp.jar \
      -patient_id ~{donor} \
      -tumor_id ~{tumour_name} \
      -samples ~{plasma_name} \
      -purple_dir ~{tumour_name}.solPrimary.purple/ \
      -amber_dir "merged_amber"/ \
      -cobalt_dir "merged_cobalt"/ \
      -somatic_vcf ~{somatic_vcf} \
      -bqr_dir ~{plasma_name}.sage.bqr/ \
      -ref_genome ~{refFasta} \
      -output_dir ~{plasma_name}.wisp/ \
      -threads ~{threads}

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
    File wisp_summary = "~{plasma_name}.wisp/~{plasma_name}.wisp.summary.tsv"
    File wisp_snv_results = "~{plasma_name}.wisp/~{plasma_name}.wisp.snv.tsv"
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
    String purpleScript = "java -Xmx8G -jar $HMFTOOLS_ROOT/purple.jar"
    String? min_ploidy
    String? max_ploidy
    String? min_purity
    String? max_purity
    String? ploidy_penalty_factor
    String? ploidy_penalty_standard_deviation
    String modules
    Int threads = 8
    Int memory = 32
    Int timeout = 100
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
    min_diploid_tumor_ratio_count: "smooth over contiguous segments which are fewer than this number of depth windows long and which have no SV support on either side and which are bounded on both sides by copy number regions which could be smoothed together using our normal smoothing rules."
    genomeVersion: "genome version for AMBER, default set to V38"
    purpleScript: "location of PURPLE script"
    min_ploidy: "minimum ploidy"
    max_ploidy: "max ploidy"
    min_purity: "mininimum purity"
    max_purity: "max purity"
    ploidy_penalty_factor: "multiplies aggregate event penalty by this factor"
    ploidy_penalty_standard_deviation: "not entirely sure what this does"
    modules: "Required environment modules"
    memory: "Memory allocated for this job (GB)"
    threads: "Requested CPU threads"
    timeout: "Hours before task timeout"
	}

  command <<<
    set -euo pipefail

    unzip ~{amber_directory} 
    unzip ~{cobalt_directory} 
    mkdir ~{outfilePrefix}.purple 

    ~{purpleScript} \
      -ref_genome_version ~{genomeVersion} \
      -ref_genome ~{refFasta}  \
      -gc_profile ~{gcProfile} \
      -ensembl_data_dir ~{ensemblDir}  \
      -reference ~{normal_name} -tumor ~{tumour_name}  \
      -amber ~{tumour_name}.amber -cobalt ~{tumour_name}.cobalt \
      ~{"-ploidy_penalty_factor" + ploidy_penalty_factor} \
      ~{"-ploidy_penalty_standard_deviation" + ploidy_penalty_standard_deviation} \
      ~{"-somatic_vcf " + somatic_vcf} \
      ~{"-min_ploidy " + min_ploidy} \
      ~{"-max_ploidy " + max_ploidy} \
      ~{"-min_purity " + min_purity} \
      ~{"-max_purity " + max_purity} \
      -no_charts \
      -min_diploid_tumor_ratio_count ~{min_diploid_tumor_ratio_count} \
      -output_dir ~{outfilePrefix}.purple 

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
  }

  meta {
    output_meta: {
      purple_directory: "Zipped Output PURPLE directory",
      purple_qc: "QC results from PURPLE",
      purple_purity: "tab seperated Purity estimate from PURPLE",
      purple_purity_range: "tab seperated range of Purity estimate from PURPLE",
      purple_segments: "tab seperated segments estimated by PURPLE",
      purple_cnv: "tab seperated somatic copy number variants from PURPLE",
      purple_cnv_gene: "tab seperated somatic gene-level copy number variants from PURPLE"
		}
	}
}
