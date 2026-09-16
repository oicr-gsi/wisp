version 1.0

# Alignment file plus its index. BAM or CRAM: REDUX names the argument -input_bam but reads
# either, given the reference, and its own output is a BAM whichever it was handed. A CRAM
# needs its index beside it, which the tasks arrange.
struct Alignment {
    File aln
    File idx
}

# wisp - SNV-based MRD detection, implemented as discrete tasks.
#
# User-facing docs (modes, valid input combinations, what each mode produces) live in
# meta.description at the bottom of this file, which is what generates README.md. Keep
# them there, not here. Notes below are for maintainers.
#
# Every task runs one tool from WiGiTS, the Hartwig Medical Foundation toolkit also known as
# hmftools, inside its own container image, with the arguments the tool needs and nothing
# else. The command in each task is derived from the corresponding
# process in nf-core/oncoanalyser 3.0.0, including the conditional arguments; when the
# upstream process changes, the task must be re-derived from the conditionals rather than
# from a single observed invocation.
#
# Copy-number and LOH evidence are out of scope. WISP is asked for SOMATIC_VARIANT and,
# when use_copy_number is set, COPY_NUMBER; AMBER_LOH needs the primary AMBER directory
# together with the primary normal alignment, which this workflow does not carry into the
# longitudinal stage.
#
# Resource locations are declarations holding literal environment-variable text, so the
# workflow carries no absolute paths and a site supplies them through a module or through
# the inputs file.

workflow wisp {
    input {
        String mode = "WG_PE"
        String outputFileNamePrefix
        String donor_id
        Array[Alignment]? tumor_alignments
        Array[Alignment]? normal_alignments
        Array[Alignment]? longitudinal_alignments
        String? tumor_redux_dir
        String? normal_redux_dir
        String? longitudinal_redux_dir
        File? primary_tarball
        String? tumor_sample_id
        String? normal_sample_id
        String? longitudinal_sample_id
        String? sequencing_platform
        String? longitudinal_sequencing_platform
        Boolean use_copy_number = true
        String hmftools_log_level = "INFO"
        Array[String] container_binds = []
        String images_dir = "$WISP_IMAGES_DIR"
        String ref_data_dir = "$WISP_REF_DATA_DIR"
        String genome_fasta = "$WISP_GENOME_FASTA"
    }

    parameter_meta {
        mode:                    "Which stages to run: WG (primary only, producing the tarball a later PE run consumes), PE (longitudinal sample against an existing primary tarball) or WG_PE (both in sequence)"
        outputFileNamePrefix:    "Prefix for the longitudinal-stage provisioned files, so runs of different samples do not provision the same name. Primary-stage files are named from the primary tumour sample id instead: one run can produce both, and naming a primary call set after the longitudinal sample reads as the wrong sample. The names the tools use among themselves always come from the sample ids"
        donor_id:                "Donor the samples came from, passed to WISP as patient_id and reported as a column of its summary. Groups a primary with every longitudinal sample drawn from the same donor, so it identifies the donor rather than the run"
        tumor_alignments:        "Primary tumour alignments and indexes, BAM or CRAM, which REDUX merges and processes. Supply these or tumor_redux_dir, not both; one of the two is MANDATORY for WG and WG_PE"
        normal_alignments:       "Matched normal alignments and indexes, BAM or CRAM. One of these or normal_redux_dir is MANDATORY for WG and WG_PE, with no override. Without a matched normal SAGE cannot subtract germline variants, the somatic call set fills with germline sites, and WISP measures those in the patient's own cfDNA and reports a large spurious tumour fraction"
        longitudinal_alignments: "Longitudinal (plasma) alignments and indexes, BAM or CRAM. Supply these or longitudinal_redux_dir, not both; one of the two is MANDATORY for PE and WG_PE"
        tumor_redux_dir:         "An existing REDUX output directory for the primary tumour, holding {sample_id}.redux.bam, its index and the recalibration, jitter and microsatellite tables. Supplied instead of tumor_alignments, so REDUX does not run again"
        normal_redux_dir:        "An existing REDUX output directory for the matched normal, supplied instead of normal_alignments"
        longitudinal_redux_dir:  "An existing REDUX output directory for the longitudinal sample, supplied instead of longitudinal_alignments"
        primary_tarball:         "Primary-stage output from an earlier WG run. MANDATORY for PE, and must not be supplied for WG or WG_PE"
        tumor_sample_id:         "Overrides the primary tumour sample id, which is otherwise read from the alignment's read-group SM tag"
        normal_sample_id:        "Overrides the matched normal sample id, which is otherwise read from the alignment's read-group SM tag"
        longitudinal_sample_id:  "Overrides the longitudinal sample id, which is otherwise read from the alignment's read-group SM tag"
        sequencing_platform:     "Platform of the primary pair: ILLUMINA, ULTIMA or SBX. Read from the read-group PL tag when not set. The tumour and the normal must agree, because AMBER and SAGE are each given both in one call"
        longitudinal_sequencing_platform: "Platform of the longitudinal sample, which may differ from the primary's. Read from its read-group PL tag when not set"
        use_copy_number:         "Whether WISP is also asked for COPY_NUMBER. When false, COBALT does not run on the longitudinal sample and WISP reports somatic-variant evidence alone"
        hmftools_log_level:      "Log level passed to every tool: ERROR, WARN, INFO, DEBUG or TRACE"
        container_binds:         "Extra host paths to bind into every container, each reduced to its filesystem root. Rarely needed, and empty is the normal case: the task directory, the reference data and wherever the alignments really live are all discovered and bound automatically"
        images_dir:              "Directory holding the container images, normally the literal $WISP_IMAGES_DIR"
        ref_data_dir:            "Root of the extracted HMF resource bundle, normally the literal $WISP_REF_DATA_DIR"
        genome_fasta:            "Reference genome FASTA, normally the literal $WISP_GENOME_FASTA. Its .fai and .dict must sit beside it"
    }

    # Only GRCh38 is supported: the resource filenames below carry the build in their
    # names and the 37 bundle names them differently.
    String genome_version = "38"

    call resolve_resources {
        input: images_dir = images_dir, ref_data_dir = ref_data_dir, genome_fasta = genome_fasta
    }

    # What the caller asked to bind, plus wherever the resources turned out to be.
    Array[String] all_binds = flatten([container_binds, resolve_resources.binds])

    String images = resolve_resources.images
    String ref_data = resolve_resources.ref_data
    String genome = resolve_resources.genome

    String driver_gene_panel     = "~{ref_data}/common/DriverGenePanel.38.tsv"
    String ensembl_data_dir      = "~{ref_data}/common/ensembl_data"
    String unmap_regions         = "~{ref_data}/common/unmap_regions.38.tsv"
    String heterozygous_sites    = "~{ref_data}/dna/copy_number/AmberGermlineSites.38.tsv.gz"
    String diploid_bed           = "~{ref_data}/dna/copy_number/DiploidRegions.38.bed.gz"
    String gc_profile            = "~{ref_data}/dna/copy_number/GC_profile.1000bp.38.cnp"
    String germline_amp_del_freq = "~{ref_data}/dna/copy_number/cohort_germline_amp_del_freq.38.csv"
    String hotspots_somatic      = "~{ref_data}/dna/variants/KnownHotspots.somatic.38.vcf.gz"
    String hotspots_germline     = "~{ref_data}/dna/variants/KnownHotspots.germline.38.vcf.gz"
    String clinvar_vcf           = "~{ref_data}/dna/variants/clinvar.38.vcf.gz"
    String gnomad_dir            = "~{ref_data}/dna/variants/gnomad"
    String mappability_bed       = "~{ref_data}/dna/variants/mappability_150.38.bed.gz"
    String msi_jitter_sites      = "~{ref_data}/dna/variants/msi_jitter_sites.38.tsv.gz"
    String high_confidence_bed   = "~{ref_data}/dna/variants/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed.gz"

    Boolean run_primary      = mode == "WG" || mode == "WG_PE"
    Boolean run_longitudinal = mode == "PE" || mode == "WG_PE"

    # ---------------------------------------------------------------------------------
    # Read what the alignments themselves report, so the sample ids, the platform and the
    # header contig order are checked rather than assumed.
    # ---------------------------------------------------------------------------------

    if (run_primary) {
        # Either raw alignments to run through REDUX, or a directory REDUX has already
        # written. Whichever it is, the alignment is read the same way afterwards.
        if (defined(tumor_redux_dir)) {
            call stage_redux_dir as stage_tumor {
                input: redux_dir = select_first([tumor_redux_dir]), role = "tumor",
                       sample_id_override = tumor_sample_id
            }
        }
        scatter (a in select_first([tumor_alignments, []])) {
            File tumor_raw_bam = a.aln
            File tumor_raw_bai = a.idx
        }
        call probe_alignments as probe_tumor {
            input: alignments = select_first([stage_tumor.alignments, tumor_raw_bam]),
                   indexes = select_first([stage_tumor.indexes, tumor_raw_bai]),
                   role = "tumor", sample_id_override = tumor_sample_id,
                   requires_mate_cigar = !defined(tumor_redux_dir)
        }

        if (defined(normal_redux_dir)) {
            call stage_redux_dir as stage_normal {
                input: redux_dir = select_first([normal_redux_dir]), role = "normal",
                       sample_id_override = normal_sample_id
            }
        }
        scatter (a in select_first([normal_alignments, []])) {
            File normal_raw_bam = a.aln
            File normal_raw_bai = a.idx
        }
        call probe_alignments as probe_normal {
            input: alignments = select_first([stage_normal.alignments, normal_raw_bam]),
                   indexes = select_first([stage_normal.indexes, normal_raw_bai]),
                   role = "normal", sample_id_override = normal_sample_id,
                   requires_mate_cigar = !defined(normal_redux_dir)
        }
    }

    if (run_longitudinal) {
        # Ahead of the preflight, so what the archive says about the primary can be checked
        # alongside everything else rather than discovered once WISP is running.
        if (defined(primary_tarball)) {
            call extract_primary {
                input: tarball = select_first([primary_tarball])
            }
        }

        if (defined(longitudinal_redux_dir)) {
            call stage_redux_dir as stage_longitudinal {
                input: redux_dir = select_first([longitudinal_redux_dir]),
                       role = "longitudinal", sample_id_override = longitudinal_sample_id
            }
        }
        scatter (a in select_first([longitudinal_alignments, []])) {
            File longitudinal_raw_bam = a.aln
            File longitudinal_raw_bai = a.idx
        }
        call probe_alignments as probe_longitudinal {
            input: alignments = select_first([stage_longitudinal.alignments, longitudinal_raw_bam]),
                   indexes = select_first([stage_longitudinal.indexes, longitudinal_raw_bai]),
                   role = "longitudinal", sample_id_override = longitudinal_sample_id,
                   requires_mate_cigar = !defined(longitudinal_redux_dir)
        }
    }

    call validate_inputs {
        input:
            mode                    = mode,
            genome_version          = genome_version,
            roles                   = select_all([probe_tumor.sample_role, probe_normal.sample_role, probe_longitudinal.sample_role]),
            sample_ids              = select_all([probe_tumor.sample_id, probe_normal.sample_id, probe_longitudinal.sample_id]),
            platforms               = select_all([probe_tumor.platform, probe_normal.platform, probe_longitudinal.platform]),
            mate_cigar_present      = select_all([probe_tumor.has_mate_cigar, probe_normal.has_mate_cigar, probe_longitudinal.has_mate_cigar]),
            contig_lists            = select_all([probe_tumor.contig_names, probe_normal.contig_names, probe_longitudinal.contig_names]),
            tumor_count             = length(select_first([tumor_alignments, []])),
            normal_count            = length(select_first([normal_alignments, []])),
            longitudinal_count      = length(select_first([longitudinal_alignments, []])),
            has_tumor_redux_dir     = defined(tumor_redux_dir),
            has_normal_redux_dir    = defined(normal_redux_dir),
            has_longitudinal_redux_dir = defined(longitudinal_redux_dir),
            has_primary_tarball     = defined(primary_tarball),
            primary_contigs         = extract_primary.contig_names,
            primary_purple_version  = extract_primary.purple_version,
            primary_tumor_id        = extract_primary.tumor_id,
            outputFileNamePrefix    = outputFileNamePrefix,
            sequencing_platform     = sequencing_platform,
            longitudinal_sequencing_platform = longitudinal_sequencing_platform,
            use_copy_number         = use_copy_number,
            genome_fasta            = genome,
            images_dir              = images,
            reference_files         = [driver_gene_panel, ensembl_data_dir, unmap_regions,
                                       heterozygous_sites, gc_profile, germline_amp_del_freq,
                                       hotspots_somatic, hotspots_germline, clinvar_vcf,
                                       gnomad_dir, mappability_bed, msi_jitter_sites,
                                       high_confidence_bed]
    }

    String primary_platform = validate_inputs.primary_platform
    String longitudinal_platform = validate_inputs.longitudinal_platform

    # SAGE and PAVE only ever see the primary, so its platform selects the panel of normals.
    # SAGE append takes no panel of normals, so the longitudinal platform needs none.
    String sage_pon = if primary_platform == "ULTIMA"
                      then "~{ref_data}/dna/variants/hmf_wgs_sage_pon.ill_ult.38.tsv.gz"
                      else if primary_platform == "SBX"
                           then "~{ref_data}/dna/variants/hmf_wgs_sage_pon.ill_sbx.38.tsv.gz"
                           else "~{ref_data}/dna/variants/hmf_wgs_sage_pon_1000.38.tsv.gz"

    # ---------------------------------------------------------------------------------
    # Primary stage: call somatic variants in the tumour against its matched normal and
    # fit them with PURPLE, producing the call set the longitudinal stage measures.
    # ---------------------------------------------------------------------------------

    if (run_primary) {
        if (! defined(tumor_redux_dir)) {
            call redux as redux_tumor {
            input:
                sample_id       = select_first([probe_tumor.sample_id]),
                alignments      = select_first([tumor_raw_bam]),
                indexes         = select_first([tumor_raw_bai]),
                platform        = primary_platform,
                genome_fasta    = genome,
                genome_version  = genome_version,
                msi_jitter_sites = msi_jitter_sites,
                unmap_regions   = unmap_regions,
                log_level       = hmftools_log_level,
                images_dir      = images,
                container_binds = all_binds,
                checked         = validate_inputs.checked
            }
        }
        File tumor_bam = select_first([redux_tumor.redux_bam, stage_tumor.bam])
        File tumor_bai = select_first([redux_tumor.redux_bai, stage_tumor.bai])
        Array[File] tumor_tsvs = select_first([redux_tumor.redux_tsvs, stage_tumor.tsvs])

        if (! defined(normal_redux_dir)) {
            call redux as redux_normal {
            input:
                sample_id       = select_first([probe_normal.sample_id]),
                alignments      = select_first([normal_raw_bam]),
                indexes         = select_first([normal_raw_bai]),
                platform        = primary_platform,
                genome_fasta    = genome,
                genome_version  = genome_version,
                msi_jitter_sites = msi_jitter_sites,
                unmap_regions   = unmap_regions,
                log_level       = hmftools_log_level,
                images_dir      = images,
                container_binds = all_binds,
                checked         = validate_inputs.checked
            }
        }
        File normal_bam = select_first([redux_normal.redux_bam, stage_normal.bam])
        File normal_bai = select_first([redux_normal.redux_bai, stage_normal.bai])
        Array[File] normal_tsvs = select_first([redux_normal.redux_tsvs, stage_normal.tsvs])

        call amber as amber_primary {
            input:
                tumor_id       = select_first([probe_tumor.sample_id]),
                tumor_bam      = tumor_bam,
                tumor_bai      = tumor_bai,
                normal_id      = probe_normal.sample_id,
                normal_bam     = normal_bam,
                normal_bai     = normal_bai,
                platform       = primary_platform,
                genome_fasta   = genome,
                genome_version = genome_version,
                heterozygous_sites = heterozygous_sites,
                checked        = validate_inputs.checked,
                log_level      = hmftools_log_level,
                images_dir     = images,
                container_binds = all_binds
        }

        call cobalt as cobalt_primary {
            input:
                tumor_id       = select_first([probe_tumor.sample_id]),
                tumor_bam      = tumor_bam,
                tumor_bai      = tumor_bai,
                normal_id      = probe_normal.sample_id,
                normal_bam     = normal_bam,
                normal_bai     = normal_bai,
                genome_fasta   = genome,
                genome_version = genome_version,
                gc_profile     = gc_profile,
                checked        = validate_inputs.checked,
                log_level      = hmftools_log_level,
                images_dir     = images,
                container_binds = all_binds
        }

        call sage_somatic {
            input:
                tumor_id       = select_first([probe_tumor.sample_id]),
                tumor_bam      = tumor_bam,
                tumor_bai      = tumor_bai,
                tumor_tsvs     = tumor_tsvs,
                normal_id      = select_first([probe_normal.sample_id]),
                normal_bam     = normal_bam,
                normal_bai     = normal_bai,
                normal_tsvs    = normal_tsvs,
                platform       = primary_platform,
                genome_fasta   = genome,
                genome_version = genome_version,
                hotspots       = hotspots_somatic,
                high_confidence_bed = high_confidence_bed,
                driver_gene_panel = driver_gene_panel,
                ensembl_data_dir = ensembl_data_dir,
                sage_pon       = sage_pon,
                gnomad_dir     = gnomad_dir,
                checked        = validate_inputs.checked,
                log_level      = hmftools_log_level,
                images_dir     = images,
                container_binds = all_binds
        }

        call pave_somatic {
            input:
                tumor_id       = select_first([probe_tumor.sample_id]),
                sage_vcf       = sage_somatic.somatic_vcf,
                sage_tbi       = sage_somatic.somatic_tbi,
                platform       = primary_platform,
                genome_fasta   = genome,
                genome_version = genome_version,
                sage_pon       = sage_pon,
                gnomad_dir     = gnomad_dir,
                clinvar_vcf    = clinvar_vcf,
                driver_gene_panel = driver_gene_panel,
                mappability_bed = mappability_bed,
                ensembl_data_dir = ensembl_data_dir,
                log_level      = hmftools_log_level,
                images_dir     = images,
                container_binds = all_binds
        }

        call purple {
            input:
                tumor_id       = select_first([probe_tumor.sample_id]),
                normal_id      = probe_normal.sample_id,
                amber_files    = amber_primary.amber_files,
                cobalt_files   = cobalt_primary.cobalt_files,
                pave_vcf       = pave_somatic.pave_vcf,
                pave_tbi       = pave_somatic.pave_tbi,
                redux_tumor_tsvs = tumor_tsvs,
                genome_fasta   = genome,
                genome_version = genome_version,
                gc_profile     = gc_profile,
                hotspots_somatic = hotspots_somatic,
                hotspots_germline = hotspots_germline,
                driver_gene_panel = driver_gene_panel,
                ensembl_data_dir = ensembl_data_dir,
                germline_amp_del_freq = germline_amp_del_freq,
                log_level      = hmftools_log_level,
                images_dir     = images,
                container_binds = all_binds
        }

        call pack_primary {
            input:
                tumor_id     = select_first([probe_tumor.sample_id]),
                purple_files = purple.purple_files,
                amber_files  = amber_primary.amber_files,
                purple_plots = purple.purple_plots,
                amber_plots  = amber_primary.amber_plots,
                cobalt_plots = cobalt_primary.cobalt_plots
        }
    }

    # ---------------------------------------------------------------------------------
    # Longitudinal stage: force-call the primary's somatic sites in the plasma sample and
    # estimate the tumour fraction they support.
    # ---------------------------------------------------------------------------------

    if (run_longitudinal) {
        String primary_id = select_first([probe_tumor.sample_id, extract_primary.tumor_id])
        Array[File] primary_purple = select_first([purple.purple_files, extract_primary.purple_files])

        if (! defined(longitudinal_redux_dir)) {
            call redux as redux_longitudinal {
            input:
                sample_id       = select_first([probe_longitudinal.sample_id]),
                alignments      = select_first([longitudinal_raw_bam]),
                indexes         = select_first([longitudinal_raw_bai]),
                platform        = longitudinal_platform,
                genome_fasta    = genome,
                genome_version  = genome_version,
                msi_jitter_sites = msi_jitter_sites,
                unmap_regions   = unmap_regions,
                log_level       = hmftools_log_level,
                images_dir      = images,
                container_binds = all_binds,
                checked         = validate_inputs.checked
            }
        }
        File longitudinal_bam = select_first([redux_longitudinal.redux_bam, stage_longitudinal.bam])
        File longitudinal_bai = select_first([redux_longitudinal.redux_bai, stage_longitudinal.bai])
        Array[File] longitudinal_tsvs = select_first([redux_longitudinal.redux_tsvs, stage_longitudinal.tsvs])

        # Tumour-only, so COBALT normalises against the diploid regions BED rather than a
        # matched normal.
        if (use_copy_number) {
            call cobalt as cobalt_longitudinal {
                input:
                    tumor_id       = select_first([probe_longitudinal.sample_id]),
                    tumor_bam      = longitudinal_bam,
                    tumor_bai      = longitudinal_bai,
                    genome_fasta   = genome,
                    genome_version = genome_version,
                    gc_profile     = gc_profile,
                    diploid_bed    = diploid_bed,
                    checked        = validate_inputs.checked,
                    log_level      = hmftools_log_level,
                    images_dir     = images,
                    container_binds = all_binds
            }
        }

        call sage_append {
            input:
                primary_id     = primary_id,
                longitudinal_id = select_first([probe_longitudinal.sample_id]),
                purple_files   = primary_purple,
                longitudinal_bam = longitudinal_bam,
                longitudinal_bai = longitudinal_bai,
                longitudinal_tsvs = longitudinal_tsvs,
                platform       = longitudinal_platform,
                genome_fasta   = genome,
                genome_version = genome_version,
                outputFileNamePrefix = outputFileNamePrefix,
                log_level      = hmftools_log_level,
                images_dir     = images,
                container_binds = all_binds
        }

        call wisp_purity {
            input:
                donor_id      = donor_id,
                primary_id      = primary_id,
                longitudinal_id = select_first([probe_longitudinal.sample_id]),
                purple_files    = primary_purple,
                append_vcf      = sage_append.append_vcf,
                append_tbi      = sage_append.append_tbi,
                longitudinal_tsvs = longitudinal_tsvs,
                cobalt_files    = cobalt_longitudinal.cobalt_files,
                use_copy_number = use_copy_number,
                genome_fasta    = genome,
                outputFileNamePrefix = outputFileNamePrefix,
                log_level       = hmftools_log_level,
                images_dir      = images,
                container_binds = all_binds
        }
    }

    output {
        File validation_log = validate_inputs.log
        File? primary_output = pack_primary.tarball
        File? primary_somatic_vcf = purple.somatic_vcf
        File? primary_purity = purple.purity_tsv
        File? longitudinal_append_vcf = sage_append.append_vcf
        File? wisp_summary = wisp_purity.summary
        File? wisp_output = wisp_purity.tarball
    }

    meta {
        author: "Gavin Peng"
        description: "SNV-based MRD detection with the Hartwig WiGiTS tools, run as discrete tasks rather than through a pipeline engine.\n\n![wisp workflow](docs/wisp.flow.svg)\n\nThe chart is a declaration-level view: every box is one Cromwell task running one tool in its own container, and the dashed clusters are the mode conditionals. Tasks that only prepare or check inputs are hidden -- `resolve_resources`, `probe_alignments`, `validate_inputs`, `stage_redux_dir`, `extract_primary`, `pack_primary`. Diagram source is Graphviz, in docs/.\n\nA primary tumour is called against its matched normal and fitted with PURPLE, producing a somatic call set. A longitudinal (plasma) sample is then force-called at exactly those sites and WISP estimates the tumour fraction they support. Copy-number and LOH evidence are deliberately out of scope: WISP is asked for SOMATIC_VARIANT and, optionally, COPY_NUMBER.\n\n## Modes\n\n| mode | inputs | produces |\n|---|---|---|\n| `WG` | `tumor_alignments`, `normal_alignments` | `primary_output` tarball for later `PE` runs |\n| `PE` | `longitudinal_alignments`, `primary_tarball` | `wisp_summary`, `wisp_output` |\n| `WG_PE` | all three alignment sets | both, in sequence |\n\n`normal_alignments` is mandatory for `WG` and `WG_PE`, with no override. Without a matched normal SAGE cannot subtract germline variants, the somatic call set fills with germline sites, and WISP measures those in the patient's own cfDNA and reports a large spurious tumour fraction.\n\nRun `WG` once per primary and `PE` once per timepoint. `WG_PE` suits a single-timepoint case, where re-running the primary costs nothing extra.\n\n## Inputs the alignments must satisfy\n\nSeveral properties are read from the alignments and checked before any expensive task runs, because getting them wrong produces a plausible but wrong answer rather than a failure:\n\n- **The tumour and the normal must share a sequencing platform.** AMBER and SAGE are each given both in one call and take a single platform. The longitudinal sample is only ever processed on its own, so it may differ -- an Illumina primary with an Ultima plasma is a valid run, and each sample is then processed with its own error model while the site list still comes from the primary. Set `sequencing_platform` and `longitudinal_sequencing_platform` to override what the read-group PL tags report.\n- **Mate CIGAR (MC) tags must be present**, wherever REDUX is going to mark duplicates. Without them it marks them wrong and reports nothing unusual. Alignments produced by bwa-mem2 carry them. Ultima is exempt, its reads being single-ended, and so is a sample supplied as a REDUX directory.\n- **The header must order the called contigs the same way the reference does.** The tools address a contig by its position in the alignment header, so a header sorted differently -- chr10 before chr2, as an alphabetically sorted reference produces -- makes them read the wrong contig and discard the evidence without failing. Only the relative order of the contigs the two have in common is checked, so a header carrying fewer decoys than the reference is fine; and a disagreement confined to the decoys that follow chr1..chrM is reported as a note rather than refused, because no variant is called there.\n\nAlignments may be BAM or CRAM.\n\n## Skipping REDUX\n\nA sample that has already been through REDUX is supplied as `tumor_redux_dir`, `normal_redux_dir` or `longitudinal_redux_dir` instead of its alignments, and REDUX does not run for it. The directory must hold `{sample_id}.redux.bam` with its index and the recalibration, jitter and microsatellite tables, all named by the same prefix; that prefix has to equal what the alignment reports as its read-group SM tag, because the tools resolve the files by prefix but read the sample from the header. Each sample takes one route or the other, never both.\n\n## Copy number\n\n`use_copy_number` adds COPY_NUMBER to the purity methods WISP applies, at the cost of running COBALT on the longitudinal sample. Somatic-variant evidence alone is what the assay reports, so the flag exists to let the two be compared."
        dependencies: [
            {name: "wisp/3.0.0", url: "https://github.com/hartwigmedical/hmftools/tree/master/wisp"},
            {name: "apptainer/1.4.5", url: "https://apptainer.org/"},
            {name: "hmftools-redux/2.0.5", url: "https://github.com/hartwigmedical/hmftools/tree/master/redux"},
            {name: "hmftools-amber/4.3", url: "https://github.com/hartwigmedical/hmftools/tree/master/amber"},
            {name: "hmftools-cobalt/3.0", url: "https://github.com/hartwigmedical/hmftools/tree/master/cobalt"},
            {name: "hmftools-sage/5.0.2", url: "https://github.com/hartwigmedical/hmftools/tree/master/sage"},
            {name: "hmftools-pave/1.9", url: "https://github.com/hartwigmedical/hmftools/tree/master/pave"},
            {name: "hmftools-purple/4.4", url: "https://github.com/hartwigmedical/hmftools/tree/master/purple"},
            {name: "hmftools-wisp/1.3.1", url: "https://github.com/hartwigmedical/hmftools/tree/master/wisp"}
        ]
        output_meta: {
            validation_log: {
                description: "What the preflight checks read from the alignments and decided, including the sample ids and the platform in use.",
                vidarr_label: "validation_log"
            },
            primary_output: {
                description: "Primary-stage PURPLE and AMBER output with the tool plots, as the tarball a later PE run consumes. Named from the primary tumour sample id. WG and WG_PE only.",
                vidarr_label: "primary_output"
            },
            primary_somatic_vcf: {
                description: "PURPLE somatic small-variant VCF for the primary tumour, the call set the longitudinal stage measures. Named from the primary tumour sample id. WG and WG_PE only.",
                vidarr_label: "primary_somatic_vcf"
            },
            primary_purity: {
                description: "PURPLE purity and ploidy fit for the primary tumour. WG and WG_PE only.",
                vidarr_label: "primary_purity"
            },
            longitudinal_append_vcf: {
                description: "The primary's somatic sites force-called in the longitudinal sample. PE and WG_PE only.",
                vidarr_label: "longitudinal_append_vcf"
            },
            wisp_summary: {
                description: "WISP purity estimate for the longitudinal sample, one row per purity method. PE and WG_PE only.",
                vidarr_label: "wisp_summary"
            },
            wisp_output: {
                description: "Full WISP output directory as a tarball, including the per-variant table and plots. PE and WG_PE only.",
                vidarr_label: "wisp_output"
            }
        }
    }
}

# Takes a directory REDUX has already written and presents its contents the way the redux
# task does, so the rest of the workflow cannot tell which produced them. The sample id comes
# from the filenames, because that is how the tools resolve these files: everything is named
# {sample_id}.redux.*, and a prefix that disagrees with the alignment's own read-group SM tag
# is caught by validate_inputs rather than here.
# The resource locations are given as environment-variable text so that the workflow carries
# no absolute paths. A shell has to expand them, and only a task runs a shell, so they are
# expanded once here and every later use is a real path. Without this the values survive into
# places no shell reads them: a list written with write_lines, or a script handed to a
# container.
task resolve_resources {
    input {
        String images_dir
        String ref_data_dir
        String genome_fasta
        Int jobMemory = 1
        Int cores = 1
        Int timeout = 1
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        images_dir:   "Directory holding the container images, as a path or as environment-variable text"
        ref_data_dir: "Root of the reference data, as a path or as environment-variable text"
        genome_fasta: "Reference genome FASTA, as a path or as environment-variable text"
        jobMemory:    "Memory allocated to the job, in GB"
        cores:        "Number of CPUs allocated to the job"
        timeout:      "Maximum run time, in hours"
        modules:      "Environment modules to load"
    }

    command <<<
        set -euo pipefail

        printf '%s' "~{images_dir}"   > images_dir.txt
        printf '%s' "~{ref_data_dir}" > ref_data_dir.txt
        printf '%s' "~{genome_fasta}" > genome_fasta.txt

        # Where the resources really live, as well as where they are named. A resource root
        # reached through a symlink has to be bound at both paths: the tools open the name they
        # were given, and the kernel resolves it to the target.
        {
            for d in "$(cat images_dir.txt)" "$(cat ref_data_dir.txt)" "$(dirname "$(cat genome_fasta.txt)")"; do
                printf '%s\n' "${d}"
                realpath "${d}" 2>/dev/null || true
            done
        } | sort -u | grep -v '^$' > binds.txt

        for f in images_dir.txt ref_data_dir.txt genome_fasta.txt; do
            value=$(cat "${f}")
            if [ -z "${value}" ]; then
                echo "ERROR: ${f%.txt} resolved to nothing; the module that sets it is probably not loaded" >&2
                exit 1
            fi
            echo "${f%.txt}: ${value}" >&2
        done
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        String images = read_string("images_dir.txt")
        String ref_data = read_string("ref_data_dir.txt")
        String genome = read_string("genome_fasta.txt")
        Array[String] binds = read_lines("binds.txt")
    }
}


task stage_redux_dir {
    input {
        String redux_dir
        String role
        String? sample_id_override
        Int jobMemory = 2
        Int cores = 1
        Int timeout = 1
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        redux_dir:          "Existing REDUX output directory, holding either a BAM with its .bai or a CRAM with its .crai"
        role:               "Which sample this is: tumor, normal or longitudinal. Used to name the sample in an error"
        sample_id_override: "Selects which sample to take when the directory holds more than one. Otherwise the directory must hold exactly one"
        jobMemory:          "Memory allocated to the job, in GB"
        cores:              "Number of CPUs allocated to the job"
        timeout:            "Maximum run time, in hours"
        modules:            "Environment modules to load"
    }

    command <<<
        set -euo pipefail

        dir="~{redux_dir}"
        [ -d "${dir}" ] || { echo "ERROR: ~{role} redux_dir is not a directory: ${dir}" >&2; exit 1; }

        # CRAM first, as the resolver upstream does, then BAM.
        wanted="~{default="" sample_id_override}"
        if [ -n "${wanted}" ]; then
            for ext in cram bam; do
                [ -f "${dir}/${wanted}.redux.${ext}" ] && aln="${dir}/${wanted}.redux.${ext}" && break
            done
            [ -n "${aln:-}" ] || {
                echo "ERROR: no ${wanted}.redux.cram or ${wanted}.redux.bam in ${dir}" >&2; exit 1; }
        else
            count=$(find "${dir}" -maxdepth 1 \( -name '*.redux.cram' -o -name '*.redux.bam' \) | wc -l)
            if [ "${count}" -ne 1 ]; then
                echo "ERROR: ~{role} redux_dir holds ${count} *.redux.cram or *.redux.bam files," \
                     "expected one; name the one to use with the sample id input: ${dir}" >&2
                exit 1
            fi
            aln=$(find "${dir}" -maxdepth 1 \( -name '*.redux.cram' -o -name '*.redux.bam' \))
        fi

        aln_ext="${aln##*.}"
        case "${aln_ext}" in
            cram) idx_ext=crai ;;
            bam)  idx_ext=bai ;;
        esac

        sample_id=$(basename "${aln}" ".redux.${aln_ext}")
        echo "${sample_id}" > sample_id.txt

        # The tools resolve every file in the directory by this exact prefix, but they read
        # the sample out of the alignment header, so the two disagreeing produces a run that
        # looks for files that are not there. Refused here rather than several hours in.
        sm=$(samtools view -H "${aln}" | awk -F'\t' '$1 == "@RG" {
                 for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) { print substr($i, 4); }
             }' | sort -u)
        if [ "${sm}" != "${sample_id}" ]; then
            echo "ERROR: ~{role} redux_dir file prefix '${sample_id}' does not match the" \
                 "read-group SM tag '${sm}'; regenerate the directory with -sample ${sm}," \
                 "or rename its files to match" >&2
            exit 1
        fi

        # The three tables SAGE reads. Their absence would silently switch it to skipping
        # recalibration and jitter fitting, so a directory without them is refused.
        mkdir -p redux
        ln -s "${aln}" "redux/${sample_id}.redux.${aln_ext}"
        for suffix in "redux.${aln_ext}.${idx_ext}" redux.bqr.tsv redux.jitter_params.tsv redux.ms_table.tsv.gz; do
            src="${dir}/${sample_id}.${suffix}"
            [ -f "${src}" ] || { echo "ERROR: ${sample_id}.${suffix} is missing from ${dir}" >&2; exit 1; }
            ln -s "${src}" "redux/${sample_id}.${suffix}"
        done

        # Whatever else REDUX wrote is carried along, since PURPLE reads the microsatellite
        # tables and none of it is large.
        for suffix in redux.duplicate_freq.tsv redux.msi_prediction.tsv redux.repeat.tsv.gz; do
            src="${dir}/${sample_id}.${suffix}"
            [ -f "${src}" ] && ln -s "${src}" "redux/${sample_id}.${suffix}"
        done

        echo "~{role}: taking ${sample_id} from ${dir}, REDUX will not run for it" >&2
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        String sample_id = read_string("sample_id.txt")
        # Matched rather than named: the sample id is only known once the task has run, and a
        # path built from a function call is one the engine cannot evaluate when it works out
        # what the job produces. Exactly one alignment is linked in, so the match is unique.
        # The patterns take both a BAM and its .bai and a CRAM and its .crai.
        File bam = glob("redux/*.redux.*am")[0]
        File bai = glob("redux/*.redux.*ai")[0]
        Array[File] alignments = glob("redux/*.redux.*am")
        Array[File] indexes = glob("redux/*.redux.*ai")
        Array[File] tsvs = flatten([glob("redux/*.tsv"), glob("redux/*.tsv.gz")])
    }
}


# Reads what the alignments report about themselves: the sample id, the sequencing
# platform, whether mate CIGAR tags are present, and the header contig order. All four are
# checked in validate_inputs rather than assumed, because each one is wrong silently.
task probe_alignments {
    input {
        Array[File] alignments
        Array[File] indexes
        String role
        String? sample_id_override
        Boolean requires_mate_cigar = true
        Int records = 10000
        Int jobMemory = 4
        Int cores = 1
        Int timeout = 2
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        alignments:        "Alignments for one sample. Every entry must report the same read-group SM tag"
        indexes:           "Indexes for the alignments, in the same order"
        role:              "Which sample this is: tumor, normal or longitudinal. Reported back so validate_inputs can name the sample in an error"
        sample_id_override: "Used as the sample id instead of the read-group SM tag. Still checked against what the alignments report, and a disagreement is reported"
        requires_mate_cigar: "Whether mate CIGAR tags have to be present. False for an alignment REDUX has already processed, which has consumed them"
        records:           "How many records to read when looking for mate CIGAR tags"
        jobMemory:         "Memory allocated to the job, in GB"
        cores:             "Number of CPUs allocated to the job"
        timeout:           "Maximum run time, in hours"
        modules:           "Environment modules to load"
    }

    command <<<
        set -euo pipefail

        # Sample id and platform come from the read groups, which must agree across the
        # inputs: they are merged into one sample downstream.
        sms=""
        pls=""
        while IFS= read -r aln; do
            sm=$(samtools view -H "${aln}" | awk -F'\t' '$1 == "@RG" {
                     for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) { print substr($i, 4); }
                 }' | sort -u)
            pl=$(samtools view -H "${aln}" | awk -F'\t' '$1 == "@RG" {
                     for (i = 2; i <= NF; i++) if ($i ~ /^PL:/) { print toupper(substr($i, 4)); }
                 }' | sort -u)
            sms="${sms}${sm}"$'\n'
            pls="${pls}${pl}"$'\n'
        done < ~{write_lines(alignments)}

        sm_unique=$(printf '%s' "${sms}" | grep -v '^$' | sort -u)
        pl_unique=$(printf '%s' "${pls}" | grep -v '^$' | sort -u)

        if [ "$(printf '%s\n' "${sm_unique}" | grep -c .)" -ne 1 ]; then
            echo "ERROR: ~{role} alignments do not agree on one read-group SM tag:" >&2
            printf '  %s\n' ${sm_unique} >&2
            exit 1
        fi

        # Tested as a string rather than through a file: `echo` of an empty override still
        # writes a newline, which reads back as a file with content and would override the
        # sample id with nothing.
        override="~{default="" sample_id_override}"
        if [ -n "${override}" ]; then
            if [ "${override}" != "${sm_unique}" ]; then
                echo "NOTE: ~{role} sample id overridden: alignments report ${sm_unique}," \
                     "using ${override}" >&2
            fi
            printf '%s' "${override}" > sample_id.txt
        else
            printf '%s' "${sm_unique}" > sample_id.txt
        fi

        # An empty or unrecognised PL is left for validate_inputs to report, together with
        # whatever the other samples said, rather than failed here one sample at a time.
        case "${pl_unique}" in
            ILLUMINA) echo ILLUMINA ;;
            ULTIMA)   echo ULTIMA ;;
            SBX|ELEMENT) echo SBX ;;
            *)        echo "UNKNOWN:${pl_unique}" ;;
        esac | head -1 > platform.txt

        # The tools address contigs by their position in the alignment header, so the
        # header order is what has to match the reference, not merely the contig names.
        samtools view -H "$(head -1 ~{write_lines(alignments)})" \
            | awk -F'\t' '$1 == "@SQ" {
                  for (i = 2; i <= NF; i++) if ($i ~ /^SN:/) { print substr($i, 4); }
              }' > contig_names.txt

        # REDUX marks duplicates from the mate CIGAR tag. Reading records rather than the
        # header, because the tag is per record. An alignment REDUX has already processed is
        # reported as satisfying the requirement without being read.
        if ~{if requires_mate_cigar then "true" else "false"}; then
            found=$(samtools view "$(head -1 ~{write_lines(alignments)})" \
                    | head -~{records} | grep -c 'MC:Z:' || true)
            if [ "${found}" -gt 0 ]; then echo true > has_mate_cigar.txt; else echo false > has_mate_cigar.txt; fi
        else
            echo true > has_mate_cigar.txt
        fi

        echo "~{role}: sample_id=$(cat sample_id.txt) platform=$(cat platform.txt)" \
             "contigs=$(wc -l < contig_names.txt) mate_cigar=$(cat has_mate_cigar.txt)" >&2
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        String sample_role = "~{role}"
        String sample_id = read_string("sample_id.txt")
        String platform = read_string("platform.txt")
        Boolean has_mate_cigar = read_boolean("has_mate_cigar.txt")
        File contig_names = "contig_names.txt"
    }
}

# Checks everything that would otherwise be discovered hours in, or not at all. Runs before
# any tool, and the primary REDUX tasks take its output so nothing expensive starts until it
# passes.
task validate_inputs {
    input {
        String mode
        String genome_version
        Array[String] roles
        Array[String] sample_ids
        Array[String] platforms
        Array[Boolean] mate_cigar_present
        Array[File] contig_lists
        Int tumor_count
        Int normal_count
        Int longitudinal_count
        Boolean has_tumor_redux_dir
        Boolean has_normal_redux_dir
        Boolean has_longitudinal_redux_dir
        Boolean has_primary_tarball
        File? primary_contigs
        String? primary_purple_version
        String? primary_tumor_id
        String outputFileNamePrefix
        String? sequencing_platform
        String? longitudinal_sequencing_platform
        Boolean use_copy_number
        String genome_fasta
        String images_dir
        Array[String] reference_files
        Int jobMemory = 2
        Int cores = 1
        Int timeout = 1
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        mode:                "Which stages the run asked for"
        genome_version:      "Reference build the resource filenames were built for. Only 38 is supported"
        roles:               "Sample roles reported by the probes, used to name a sample in an error"
        sample_ids:          "Sample ids reported by the probes, in the same order as roles"
        platforms:           "Sequencing platforms reported by the probes, in the same order as roles"
        mate_cigar_present:  "Whether each sample's alignments carry mate CIGAR tags, in the same order as roles"
        contig_lists:        "Header contig names for each sample, in header order, in the same order as roles"
        tumor_count:         "How many primary tumour alignments were supplied"
        normal_count:        "How many matched normal alignments were supplied"
        longitudinal_count:  "How many longitudinal alignments were supplied"
        has_tumor_redux_dir: "Whether an existing REDUX directory was supplied for the primary tumour"
        has_normal_redux_dir: "Whether an existing REDUX directory was supplied for the matched normal"
        has_longitudinal_redux_dir: "Whether an existing REDUX directory was supplied for the longitudinal sample"
        has_primary_tarball: "Whether a primary tarball was supplied"
        primary_contigs:     "Contigs the primary call set was made against, in the order its caller used. In PE mode there are no primary alignments to read, so this is the only evidence of what the primary was called with"
        primary_purple_version: "PURPLE version that produced the primary call set, recorded in the log"
        primary_tumor_id:    "Primary tumour sample id read out of the archive, recorded in the log"
        outputFileNamePrefix: "Prefix for the provisioned log"
        sequencing_platform: "Primary platform override, which wins over what the tumour and normal read groups report"
        longitudinal_sequencing_platform: "Longitudinal platform override, which wins over what its read groups report"
        use_copy_number:     "Whether COPY_NUMBER was requested, which needs the diploid regions BED"
        genome_fasta:        "Reference genome FASTA, whose .fai and .dict must sit beside it"
        images_dir:          "Directory that must hold every container image the run needs"
        reference_files:     "Reference data paths that must be readable"
        jobMemory:           "Memory allocated to the job, in GB"
        cores:               "Number of CPUs allocated to the job"
        timeout:             "Maximum run time, in hours"
        modules:             "Environment modules to load"
    }

    Array[String] required_images = ["hmftools-redux-2.0.5--hdfd78af_0.img",
                                    "hmftools-amber-4.3--hdfd78af_0.img",
                                    "hmftools-cobalt-3.0--hdfd78af_0.img",
                                    "hmftools-sage-5.0.2--hdfd78af_0.img",
                                    "hmftools-pave-1.9--hdfd78af_0.img",
                                    "hmftools-purple-4.4--hdfd78af_0.img",
                                    "hmftools-wisp-1.3.1--hdfd78af_0.img"]

    command <<<
        set -euo pipefail
        errors=()

        # chr1..chr22, chrX, chrY, chrM: the contigs variants are called on, and the size of
        # the default sequence dictionary htsjdk builds when a tool supplies none.
        MAIN_CONTIGS=25

        mode="~{mode}"
        case "${mode}" in
            WG|PE|WG_PE) ;;
            *) errors+=("mode must be WG, PE or WG_PE, got '${mode}'") ;;
        esac

        if [ "~{genome_version}" != "38" ]; then
            errors+=("only GRCh38 is supported, got genome_version '~{genome_version}'")
        fi

        # Input combinations. A tumour-only primary is refused rather than degraded: with no
        # matched normal the somatic call set fills with germline sites and the reported
        # tumour fraction is large and wrong.
        # Each sample arrives either as alignments REDUX has still to process, or as a
        # directory REDUX has already written. Exactly one of the two, so that a run cannot
        # silently ignore half of what was supplied.
        check_source() {
            local role="$1" count="$2" has_dir="$3" required="$4"
            if [ "${count}" -gt 0 ] && [ "${has_dir}" = "true" ]; then
                errors+=("${role}: supply either ${role}_alignments or ${role}_redux_dir, not both")
            elif [ "${count}" -eq 0 ] && [ "${has_dir}" != "true" ] && [ "${required}" = "true" ]; then
                errors+=("${mode} needs ${role}_alignments or ${role}_redux_dir")
            elif [ "${count}" -gt 0 ] || [ "${has_dir}" = "true" ]; then
                [ "${required}" = "true" ] || errors+=("${mode} does not use ${role}_alignments or ${role}_redux_dir")
            fi
        }

        primary_required=false
        longitudinal_required=false
        case "${mode}" in
            WG)    primary_required=true ;;
            PE)    longitudinal_required=true ;;
            WG_PE) primary_required=true; longitudinal_required=true ;;
        esac

        check_source tumor "~{tumor_count}" "~{if has_tumor_redux_dir then "true" else "false"}" "${primary_required}"
        check_source normal "~{normal_count}" "~{if has_normal_redux_dir then "true" else "false"}" "${primary_required}"
        check_source longitudinal "~{longitudinal_count}" "~{if has_longitudinal_redux_dir then "true" else "false"}" "${longitudinal_required}"

        if [ "${primary_required}" = "true" ] && [ "~{normal_count}" -eq 0 ] \
           && [ "~{if has_normal_redux_dir then "true" else "false"}" != "true" ]; then
            errors+=("a tumour-only primary inverts the MRD result rather than degrading it, so a matched normal is required")
        fi

        if [ "${mode}" = "PE" ]; then
            ~{if has_primary_tarball then "true" else "false"} || errors+=("PE needs primary_tarball")
        elif ~{if has_primary_tarball then "true" else "false"}; then
            errors+=("${mode} computes the primary itself, so primary_tarball must not be supplied")
        fi

        # Sample ids have to be distinct: the tools address samples by id and several write
        # into a shared directory keyed by it.
        dupes=$(sort "~{write_lines(sample_ids)}" | uniq -d)
        if [ -n "${dupes}" ]; then
            errors+=("sample ids must be distinct, these repeat: $(echo ${dupes} | tr '\n' ' ')")
        fi
        while IFS= read -r sid; do
            [ -n "${sid}" ] || errors+=("a sample reported an empty sample id")
        done < "~{write_lines(sample_ids)}"

        # The primary pair has to agree, because AMBER and SAGE are each given the tumour
        # and the normal in one call and take a single -sequencing_type. The longitudinal
        # sample is only ever processed on its own, so it is free to differ.
        check_platform() {
            local what="$1" value="$2"
            case "${value}" in
                ILLUMINA|ULTIMA|SBX) printf '%s' "${value}" ;;
                '') errors+=("cannot tell the ${what} sequencing platform from the read-group PL tag; set the platform input") ;;
                *)  errors+=("${what} sequencing platform must be ILLUMINA, ULTIMA or SBX, got '${value}'") ;;
            esac
        }

        roles_platforms=$(paste "~{write_lines(roles)}" "~{write_lines(platforms)}")

        primary_reported=$(printf '%s\n' "${roles_platforms}" \
            | awk -F'\t' '$1 == "tumor" || $1 == "normal" { print $2 }' | sort -u)
        long_reported=$(printf '%s\n' "${roles_platforms}" \
            | awk -F'\t' '$1 == "longitudinal" { print $2 }' | sort -u)

        if [ "$(printf '%s\n' "${primary_reported}" | grep -c . || true)" -gt 1 ]; then
            errors+=("the tumour and the normal report different sequencing platforms ($(echo ${primary_reported} | tr '\n' ' ')); they are called together, so they must agree -- set sequencing_platform")
            primary_reported=""
        fi

        primary_override="~{default="" sequencing_platform}"
        long_override="~{default="" longitudinal_sequencing_platform}"

        primary_platform=""
        if [ -n "${primary_override}" ]; then
            primary_platform=$(check_platform "primary" "$(echo "${primary_override}" | tr '[:lower:]' '[:upper:]')")
        elif [ -n "${primary_reported}" ]; then
            primary_platform=$(check_platform "primary" "${primary_reported}")
        fi

        longitudinal_platform=""
        if [ -n "${long_override}" ]; then
            longitudinal_platform=$(check_platform "longitudinal" "$(echo "${long_override}" | tr '[:lower:]' '[:upper:]')")
        elif [ -n "${long_reported}" ]; then
            longitudinal_platform=$(check_platform "longitudinal" "${long_reported}")
        fi

        if [ -n "${primary_platform}" ] && [ -n "${longitudinal_platform}" ] \
           && [ "${primary_platform}" != "${longitudinal_platform}" ]; then
            echo "NOTE: primary is ${primary_platform} and longitudinal is ${longitudinal_platform}." \
                 "Each is processed with its own error model; the site list still comes from the primary." >&2
        fi

        echo "${primary_platform}" > primary_platform.txt
        echo "${longitudinal_platform}" > longitudinal_platform.txt

        # REDUX marks duplicates from the mate CIGAR tag, so missing tags would produce a
        # plausible BAM with wrong duplicate flags. The requirement applies only where REDUX
        # will actually mark duplicates: not to a sample already processed, and not to Ultima,
        # whose reads are single-ended and therefore have no mates to describe.
        paste "~{write_lines(roles)}" "~{write_lines(platforms)}" "~{write_lines(mate_cigar_present)}" \
            | while IFS=$'\t' read -r role role_platform present; do
                  [ "${present}" != "true" ] || continue
                  [ "${role_platform}" != "ULTIMA" ] || continue
                  case "${role}" in
                      tumor)        [ "~{if has_tumor_redux_dir then "true" else "false"}" = "true" ] && continue ;;
                      normal)       [ "~{if has_normal_redux_dir then "true" else "false"}" = "true" ] && continue ;;
                      longitudinal) [ "~{if has_longitudinal_redux_dir then "true" else "false"}" = "true" ] && continue ;;
                  esac
                  echo "${role}"
              done > missing_mc.txt
        if [ -s missing_mc.txt ]; then
            errors+=("no mate CIGAR (MC) tags in the alignments for: $(tr '\n' ' ' < missing_mc.txt)-- REDUX needs them to mark duplicates correctly")
        fi

        # Contig order. The tools address a contig by its position in the alignment header, so
        # a header that orders them differently from the reference makes them read the wrong
        # contig and discard the evidence without failing. What matters is the relative order
        # of the contigs the two have in common: a header carrying fewer decoys than the
        # reference is ordinary and harmless, while one sorted differently -- chr10 before
        # chr2, as an alphabetically sorted reference produces -- is not.
        fai="~{genome_fasta}.fai"
        if [ -r "${fai}" ]; then
            cut -f1 "${fai}" > reference_contigs.txt

            # Reports the first contig that appears earlier than one already seen. A
            # disagreement among the contigs variants are called on is a different matter from
            # one among the decoys that follow them, so the two are labelled apart.
            check_order() {
                awk -v mains="${MAIN_CONTIGS}" '
                     NR == FNR { idx[$1] = FNR; next }
                     {
                         if (!($1 in idx)) { absent = absent " " $1; next }
                         if (idx[$1] < prev_idx) {
                             printf "%s %s after %s (reference positions %d and %d)\n", \
                                    (idx[$1] <= mains || prev_idx <= mains ? "MAIN:" : "DECOY:"), \
                                    $1, prev_name, idx[$1], prev_idx
                             exit
                         }
                         prev_idx = idx[$1]; prev_name = $1
                     }
                     END { if (absent != "") printf "ABSENT:%s\n", absent }' \
                    reference_contigs.txt "$1"
            }

            paste "~{write_lines(roles)}" "~{write_lines(contig_lists)}" \
                | while IFS=$'\t' read -r role contigs; do
                      check_order "${contigs}" | while IFS= read -r line; do
                          echo "${role}|${line}"
                      done
                  done > contig_report.txt

            while IFS='|' read -r role detail; do
                [ -n "${role}" ] || continue
                case "${detail}" in
                    ABSENT:*)
                        echo "NOTE: the ${role} alignment header lists contigs the reference does" \
                             "not have:${detail#ABSENT:}. Reads on them cannot be called." >&2 ;;
                    DECOY:*)
                        echo "NOTE: the ${role} alignment header orders decoy contigs differently" \
                             "from the reference:${detail#DECOY:}. The called contigs agree, so" \
                             "this does not affect the variants reported." >&2 ;;
                    MAIN:*)
                        errors+=("the ${role} alignment header orders the called contigs differently from the reference:${detail#MAIN:}; the tools address a contig by its position in the header, so read evidence would be silently discarded") ;;
                esac
            done < contig_report.txt

            # The primary arrives already called, so its alignments cannot be checked. Its
            # call set carries the dictionary its caller used, which shows the same defect.
            primary_contigs="~{default="" primary_contigs}"
            if [ -n "${primary_contigs}" ] && [ -s "${primary_contigs}" ]; then
                check_order "${primary_contigs}" > primary_report.txt
                while IFS= read -r detail; do
                    [ -n "${detail}" ] || continue
                    case "${detail}" in
                        ABSENT:*|DECOY:*) ;;
                        MAIN:*) errors+=("the primary call set in primary_tarball was made against a reference that orders the called contigs differently from this run's:${detail#MAIN:}; its variant calls would have been made with read evidence silently discarded") ;;
                    esac
                done < primary_report.txt
            fi
        else
            errors+=("reference index not readable: ${fai}")
        fi

        # Reference genome and its companions. The tools resolve the sequence dictionary
        # beside the FASTA, with the extension replaced rather than appended.
        fasta="~{genome_fasta}"
        for f in "${fasta}" "${fasta}.fai" "${fasta%.*}.dict"; do
            [ -r "${f}" ] || errors+=("reference file not readable: ${f}")
        done

        # Reference data. A directory is checked for readability, a file for both.
        while IFS= read -r f; do
            [ -n "${f}" ] || continue
            if [ -d "${f}" ]; then
                [ -r "${f}" ] || errors+=("reference directory not readable: ${f}")
            else
                [ -f "${f}" ] && [ -r "${f}" ] || errors+=("reference file not readable: ${f}")
            fi
        done < "~{write_lines(reference_files)}"

        if ~{if use_copy_number then "true" else "false"}; then
            echo "copy-number evidence requested; COBALT will run on the longitudinal sample" >&2
        fi

        # Container images. Cheap to check, and a missing one otherwise surfaces only when
        # that task starts.
        while IFS= read -r img; do
            [ -n "${img}" ] || continue
            [ -r "~{images_dir}/${img}" ] || errors+=("container image not readable: ~{images_dir}/${img}")
        done < "~{write_lines(required_images)}"

        if [ "${#errors[@]}" -gt 0 ]; then
            echo "ERROR: inputs rejected" >&2
            printf '  - %s\n' "${errors[@]}" >&2
            exit 1
        fi

        {
            echo "mode: ${mode}"
            echo "primary_platform: ${primary_platform}"
            echo "longitudinal_platform: ${longitudinal_platform}"
            echo "copy_number: ~{if use_copy_number then "true" else "false"}"
            paste "~{write_lines(roles)}" "~{write_lines(sample_ids)}" \
                | while IFS=$'\t' read -r role sid; do echo "${role}: ${sid}"; done
            if ~{if has_primary_tarball then "true" else "false"}; then
                echo "primary from tarball: ~{default="unknown" primary_tumor_id}"
                echo "primary purple version: ~{default="unknown" primary_purple_version}"
                echo "primary call set contigs: $(grep -c . "${primary_contigs:-/dev/null}" 2>/dev/null || echo 0)"
            fi
            echo "reference: ~{genome_fasta}"
            echo "contigs: $(wc -l < reference_contigs.txt)"
        } | tee ~{outputFileNamePrefix}.validation.log
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        String checked = read_string(stdout())
        String primary_platform = read_string("primary_platform.txt")
        String longitudinal_platform = read_string("longitudinal_platform.txt")
        File log = "~{outputFileNamePrefix}.validation.log"
    }
}

# REDUX marks duplicates, builds consensus reads and writes the base-quality recalibration
# and microsatellite jitter tables every downstream caller reads. Several alignments are
# passed as one list, which REDUX merges.
task redux {
    input {
        String sample_id
        Array[File] alignments
        Array[File] indexes
        String platform
        String genome_fasta
        String genome_version
        String msi_jitter_sites
        String unmap_regions
        String log_level
        String images_dir
        Array[String] container_binds
        String? checked
        String image = "hmftools-redux-2.0.5--hdfd78af_0.img"
        Float heapFraction = 0.75
        Int jobMemory = 48
        Int cores = 8
        Int timeout = 48
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        sample_id:        "Sample id, used for every output filename. Downstream tasks resolve the alignment by this exact prefix"
        alignments:       "Alignments for this sample. Several are merged by REDUX itself"
        indexes:          "Indexes for the alignments, in the same order"
        platform:         "ILLUMINA, ULTIMA or SBX. Selects whether duplicates are marked and consensus reads are formed"
        genome_fasta:     "Reference genome FASTA"
        genome_version:   "Reference build, 37 or 38"
        msi_jitter_sites: "Microsatellite sites REDUX fits the jitter model on"
        unmap_regions:    "Regions whose reads are unmapped before duplicate marking"
        log_level:        "Log level passed to the tool"
        images_dir:       "Directory holding the container images"
        container_binds: "Host paths to bind into the container. Each is reduced to its filesystem root, so naming a directory below one already bound is harmless"
        checked:          "Preflight result, taken only so that no expensive task starts before validate_inputs passes"
        image:            "Container image filename within images_dir"
        heapFraction:     "Fraction of jobMemory given to the JVM heap. The remainder covers the helper processes the tool forks, which are charged to the same allocation"
        jobMemory:        "Memory allocated to the job, in GB"
        cores:            "Number of CPUs allocated to the job"
        timeout:          "Maximum run time, in hours"
        modules:          "Environment modules to load"
    }

    Int heapMb = floor(jobMemory * 1024 * heapFraction)

    command <<<
        set -euo pipefail

        # Everything the task will open is bound by its filesystem root: the task directory
        # itself, whatever the caller asked for, and wherever the engine's input links
        # actually resolve to. Roots rather than the directories themselves, because binding a
        # path nested inside another bind makes the container runtime skip it, and one root
        # covers every file beneath it anyway.
        bind_args=(--pwd "$(pwd)")
        bind_seen=()
        add_bind_root() {
            local p real root
            p="${1:-}"
            [ -n "${p}" ] || return 0
            real=$(readlink -f "${p}") || return 0
            [ -e "${real}" ] || return 0
            root="/$(printf '%s' "${real#/}" | cut -d/ -f1)"
            [ "${root}" != "/" ] || return 0
            case " ${bind_seen[*]:-} " in *" ${root} "*) return 0 ;; esac
            bind_seen+=("${root}")
            bind_args+=(-B "${root}")
        }

        add_bind_root "$(pwd)"
        while IFS= read -r b; do add_bind_root "${b}"; done < ~{write_lines(container_binds)}

        # The engine localizes an alignment and its index into separate directories, but the
        # tools resolve an index by its position beside the alignment, so both are linked into
        # one directory here under their own names.
        mkdir -p aln redux
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -sf "${f}" "aln/$(basename "${f}")"; done < ~{write_lines(alignments)}
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -sf "${f}" "aln/$(basename "${f}")"; done < ~{write_lines(indexes)}
        input_bams=$(while IFS= read -r f; do
                         printf '%s\n' "aln/$(basename "${f}")"
                     done < ~{write_lines(alignments)} | paste -sd,)

        # Ultima reads are single-ended, so there are no duplicates to mark and no mates to
        # build a consensus from; every other platform gets both.
        cat > redux.sh <<'COMMAND'
        set -euo pipefail

        for f in aln/*; do
            [ -r "${f}" ] && continue
            echo "ERROR: ${f} cannot be read inside the container; it resolves to" \
                 "$(readlink -f "${f}")" >&2
            exit 1
        done

        redux \
            -Xmx~{heapMb}m \
            -sample ~{sample_id} \
            -input_bam INPUT_BAMS \
            -ref_genome ~{genome_fasta} \
            -ref_genome_version ~{genome_version} \
            -ref_genome_msi_file ~{msi_jitter_sites} \
            -unmap_regions ~{unmap_regions} \
            -bamtool $(which samtools) \
            -sequencing_type ~{platform} \
            -bqr_write_plot \
            PLATFORM_ARGS \
            -threads ~{cores} \
            -log_level ~{log_level} \
            -output_bam redux/~{sample_id}.redux.bam \
            -output_dir redux/
COMMAND

        if [ "~{platform}" = "ULTIMA" ]; then
            platform_args="-skip_duplicate_marking"
        else
            platform_args="-form_consensus"
        fi

        sed -i "s|INPUT_BAMS|${input_bams}|; s|PLATFORM_ARGS|${platform_args}|" redux.sh

        while IFS= read -r link; do add_bind_root "${link}"; done < <(find . -type l)

        apptainer exec "${bind_args[@]}" \
            "~{images_dir}/~{image}" bash redux.sh

        [ -f "redux/~{sample_id}.redux.bam" ]
        [ -f "redux/~{sample_id}.redux.bam.bai" ]
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        File redux_bam = "redux/~{sample_id}.redux.bam"
        File redux_bai = "redux/~{sample_id}.redux.bam.bai"
        Array[File] redux_tsvs = flatten([glob("redux/*.tsv"), glob("redux/*.tsv.gz")])
        Array[File] redux_plots = glob("redux/*.png")
    }
}

# AMBER measures B-allele frequencies at germline heterozygous sites, which PURPLE fits
# alongside the read ratios COBALT produces.
task amber {
    input {
        String tumor_id
        File tumor_bam
        File tumor_bai
        String? normal_id
        File? normal_bam
        File? normal_bai
        String platform
        String genome_fasta
        String genome_version
        String heterozygous_sites
        String? checked
        String log_level
        String images_dir
        Array[String] container_binds
        Int? tumor_min_depth
        String image = "hmftools-amber-4.3--hdfd78af_0.img"
        Float heapFraction = 0.75
        Int jobMemory = 32
        Int cores = 8
        Int timeout = 24
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        tumor_id:           "Tumour sample id"
        tumor_bam:          "Tumour REDUX alignment"
        tumor_bai:          "Index for the tumour alignment"
        normal_id:          "Matched normal sample id. Omitted for a tumour-only run"
        normal_bam:         "Matched normal REDUX alignment"
        normal_bai:         "Index for the matched normal alignment"
        platform:           "ILLUMINA, ULTIMA or SBX"
        genome_fasta:       "Reference genome FASTA"
        genome_version:     "Reference build, 37 or 38"
        heterozygous_sites: "Germline heterozygous sites to measure"
        checked:            "Preflight result, taken only so that no expensive task starts before validate_inputs passes"
        log_level:          "Log level passed to the tool"
        images_dir:         "Directory holding the container images"
        container_binds:  "Host paths to bind into the container"
        tumor_min_depth:    "Minimum tumour depth for a site to be used. Left unset for a primary, where the default applies"
        image:              "Container image filename within images_dir"
        heapFraction:       "Fraction of jobMemory given to the JVM heap"
        jobMemory:          "Memory allocated to the job, in GB"
        cores:              "Number of CPUs allocated to the job"
        timeout:            "Maximum run time, in hours"
        modules:            "Environment modules to load"
    }

    Int heapMb = floor(jobMemory * 1024 * heapFraction)
    String tumor_ext = sub(basename(tumor_bam), "^.*\\.", "")
    String tumor_idx_ext = sub(basename(tumor_bai), "^.*\\.", "")
    String normal_ext = sub(basename(select_first([normal_bam, tumor_bam])), "^.*\\.", "")
    String normal_idx_ext = sub(basename(select_first([normal_bai, tumor_bai])), "^.*\\.", "")

    command <<<
        set -euo pipefail

        # Everything the task will open is bound by its filesystem root: the task directory
        # itself, whatever the caller asked for, and wherever the engine's input links
        # actually resolve to. Roots rather than the directories themselves, because binding a
        # path nested inside another bind makes the container runtime skip it, and one root
        # covers every file beneath it anyway.
        bind_args=(--pwd "$(pwd)")
        bind_seen=()
        add_bind_root() {
            local p real root
            p="${1:-}"
            [ -n "${p}" ] || return 0
            real=$(readlink -f "${p}") || return 0
            [ -e "${real}" ] || return 0
            root="/$(printf '%s' "${real#/}" | cut -d/ -f1)"
            [ "${root}" != "/" ] || return 0
            case " ${bind_seen[*]:-} " in *" ${root} "*) return 0 ;; esac
            bind_seen+=("${root}")
            bind_args+=(-B "${root}")
        }

        add_bind_root "$(pwd)"
        while IFS= read -r b; do add_bind_root "${b}"; done < ~{write_lines(container_binds)}

        # The tools find an index by its position beside the alignment, so both are linked
        # into the task directory under their expected names.
        ln -s "~{tumor_bam}" "~{tumor_id}.redux.~{tumor_ext}"
        ln -s "~{tumor_bai}" "~{tumor_id}.redux.~{tumor_ext}.~{tumor_idx_ext}"
        ~{if defined(normal_bam) then "ln -s \"" + normal_bam + "\" \"" + normal_id + ".redux." + normal_ext + "\"" else ""}
        ~{if defined(normal_bai) then "ln -s \"" + normal_bai + "\" \"" + normal_id + ".redux." + normal_ext + "." + normal_idx_ext + "\"" else ""}

        mkdir -p amber

        cat > amber.sh <<'COMMAND'
        set -euo pipefail
        amber \
            -Xmx~{heapMb}m \
            -tumor ~{tumor_id} \
            -tumor_bam ~{tumor_id}.redux.~{tumor_ext} \
            ~{if defined(normal_id) then "-reference " + normal_id else ""} \
            ~{if defined(normal_bam) then "-reference_bam " + normal_id + ".redux." + normal_ext else ""} \
            -ref_genome ~{genome_fasta} \
            -ref_genome_version ~{genome_version} \
            -sequencing_type ~{platform} \
            -loci ~{heterozygous_sites} \
            ~{"-tumor_min_depth " + tumor_min_depth} \
            -log_level ~{log_level} \
            -threads ~{cores} \
            -output_dir amber/
COMMAND

        while IFS= read -r link; do add_bind_root "${link}"; done < <(find . -type l)

        apptainer exec "${bind_args[@]}" \
            "~{images_dir}/~{image}" bash amber.sh

        # Cromwell delocalizes files, not directories, so any subdirectory the tool wrote is
        # archived rather than listed as an output.
        mkdir -p amber_extra
        find amber -mindepth 1 -maxdepth 1 -type d -exec mv {} amber_extra/ \;
        tar -czf ~{tumor_id}.amber_plots.tar.gz -C amber_extra .
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        Array[File] amber_files = glob("amber/*")
        File amber_plots = "~{tumor_id}.amber_plots.tar.gz"
    }
}

# COBALT measures read-depth ratios in 1kb windows. With a matched normal the ratio is
# tumour over normal; without one it is normalised against the diploid regions BED, which is
# how the longitudinal sample is profiled.
task cobalt {
    input {
        String tumor_id
        File tumor_bam
        File tumor_bai
        String? normal_id
        File? normal_bam
        File? normal_bai
        String genome_fasta
        String genome_version
        String gc_profile
        String? diploid_bed
        String? checked
        String log_level
        String images_dir
        Array[String] container_binds
        String image = "hmftools-cobalt-3.0--hdfd78af_0.img"
        Float heapFraction = 0.75
        Int jobMemory = 32
        Int cores = 8
        Int timeout = 24
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        tumor_id:          "Sample being profiled. The longitudinal sample is the tumour in the longitudinal stage"
        tumor_bam:         "REDUX alignment for the sample being profiled"
        tumor_bai:         "Index for that alignment"
        normal_id:         "Matched normal sample id. Omitted for a tumour-only run"
        normal_bam:        "Matched normal REDUX alignment"
        normal_bai:        "Index for the matched normal alignment"
        genome_fasta:      "Reference genome FASTA"
        genome_version:    "Reference build, 37 or 38"
        gc_profile:        "GC content per window, used to correct the ratios"
        diploid_bed:       "Diploid regions to normalise against. MANDATORY for a tumour-only run and unused otherwise"
        checked:           "Preflight result, taken only so that no expensive task starts before validate_inputs passes"
        log_level:         "Log level passed to the tool"
        images_dir:        "Directory holding the container images"
        container_binds: "Host paths to bind into the container. Each is reduced to its filesystem root, so naming a directory below one already bound is harmless"
        image:             "Container image filename within images_dir"
        heapFraction:      "Fraction of jobMemory given to the JVM heap"
        jobMemory:         "Memory allocated to the job, in GB"
        cores:             "Number of CPUs allocated to the job"
        timeout:           "Maximum run time, in hours"
        modules:           "Environment modules to load"
    }

    Int heapMb = floor(jobMemory * 1024 * heapFraction)
    String tumor_ext = sub(basename(tumor_bam), "^.*\\.", "")
    String tumor_idx_ext = sub(basename(tumor_bai), "^.*\\.", "")
    String normal_ext = sub(basename(select_first([normal_bam, tumor_bam])), "^.*\\.", "")
    String normal_idx_ext = sub(basename(select_first([normal_bai, tumor_bai])), "^.*\\.", "")

    command <<<
        set -euo pipefail

        # Everything the task will open is bound by its filesystem root: the task directory
        # itself, whatever the caller asked for, and wherever the engine's input links
        # actually resolve to. Roots rather than the directories themselves, because binding a
        # path nested inside another bind makes the container runtime skip it, and one root
        # covers every file beneath it anyway.
        bind_args=(--pwd "$(pwd)")
        bind_seen=()
        add_bind_root() {
            local p real root
            p="${1:-}"
            [ -n "${p}" ] || return 0
            real=$(readlink -f "${p}") || return 0
            [ -e "${real}" ] || return 0
            root="/$(printf '%s' "${real#/}" | cut -d/ -f1)"
            [ "${root}" != "/" ] || return 0
            case " ${bind_seen[*]:-} " in *" ${root} "*) return 0 ;; esac
            bind_seen+=("${root}")
            bind_args+=(-B "${root}")
        }

        add_bind_root "$(pwd)"
        while IFS= read -r b; do add_bind_root "${b}"; done < ~{write_lines(container_binds)}

        ln -s "~{tumor_bam}" "~{tumor_id}.redux.~{tumor_ext}"
        ln -s "~{tumor_bai}" "~{tumor_id}.redux.~{tumor_ext}.~{tumor_idx_ext}"
        ~{if defined(normal_bam) then "ln -s \"" + normal_bam + "\" \"" + normal_id + ".redux." + normal_ext + "\"" else ""}
        ~{if defined(normal_bai) then "ln -s \"" + normal_bai + "\" \"" + normal_id + ".redux." + normal_ext + "." + normal_idx_ext + "\"" else ""}

        mkdir -p cobalt

        cat > cobalt.sh <<'COMMAND'
        set -euo pipefail
        cobalt \
            -Xmx~{heapMb}m \
            -tumor ~{tumor_id} \
            -tumor_bam ~{tumor_id}.redux.~{tumor_ext} \
            ~{if defined(normal_id) then "-reference " + normal_id else ""} \
            ~{if defined(normal_bam) then "-reference_bam " + normal_id + ".redux." + normal_ext else ""} \
            -ref_genome ~{genome_fasta} \
            -ref_genome_version ~{genome_version} \
            -gc_profile ~{gc_profile} \
            ~{if !defined(normal_bam) then "-tumor_only_diploid_bed " + diploid_bed else ""} \
            -log_level ~{log_level} \
            -threads ~{cores} \
            -output_dir cobalt/
COMMAND

        while IFS= read -r link; do add_bind_root "${link}"; done < <(find . -type l)

        apptainer exec "${bind_args[@]}" \
            "~{images_dir}/~{image}" bash cobalt.sh

        # Cromwell delocalizes files, not directories, so any subdirectory the tool wrote is
        # archived rather than listed as an output.
        mkdir -p cobalt_extra
        find cobalt -mindepth 1 -maxdepth 1 -type d -exec mv {} cobalt_extra/ \;
        tar -czf ~{tumor_id}.cobalt_plots.tar.gz -C cobalt_extra .
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        Array[File] cobalt_files = glob("cobalt/*")
        File cobalt_plots = "~{tumor_id}.cobalt_plots.tar.gz"
    }
}

# SAGE discovers somatic small variants in the tumour, using the matched normal to subtract
# germline sites. This is the call set the whole assay rests on.
task sage_somatic {
    input {
        String tumor_id
        File tumor_bam
        File tumor_bai
        Array[File] tumor_tsvs
        String normal_id
        File normal_bam
        File normal_bai
        Array[File] normal_tsvs
        String platform
        String genome_fasta
        String genome_version
        String hotspots
        String high_confidence_bed
        String driver_gene_panel
        String ensembl_data_dir
        String sage_pon
        String gnomad_dir
        String? checked
        String log_level
        String images_dir
        Array[String] container_binds
        String image = "hmftools-sage-5.0.2--hdfd78af_0.img"
        Float heapFraction = 0.75
        Int jobMemory = 80
        Int cores = 12
        Int timeout = 72
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        tumor_id:            "Tumour sample id, which also prefixes the output VCF"
        tumor_bam:           "Tumour REDUX alignment"
        tumor_bai:           "Index for the tumour alignment"
        tumor_tsvs:          "Tumour REDUX tables. SAGE reads the recalibration and jitter tables from the alignment's own directory, so these are placed beside it"
        normal_id:           "Matched normal sample id"
        normal_bam:          "Matched normal REDUX alignment"
        normal_bai:          "Index for the matched normal alignment"
        normal_tsvs:         "Matched normal REDUX tables, placed beside its alignment for the same reason"
        platform:            "ILLUMINA, ULTIMA or SBX. Also decides whether the tumour-in-normal check runs"
        genome_fasta:        "Reference genome FASTA"
        genome_version:      "Reference build, 37 or 38"
        hotspots:            "Known somatic hotspots, called with relaxed thresholds"
        high_confidence_bed: "High-confidence regions, where stricter germline filters apply"
        driver_gene_panel:   "Driver gene panel"
        ensembl_data_dir:    "Ensembl gene and transcript tables"
        sage_pon:            "Panel of normals for the sequencing platform in use"
        gnomad_dir:          "Population frequencies, read by the tumour-in-normal check"
        checked:             "Preflight result, taken only so that no expensive task starts before validate_inputs passes"
        log_level:           "Log level passed to the tool"
        images_dir:          "Directory holding the container images"
        container_binds:   "Host paths to bind into the container"
        image:               "Container image filename within images_dir"
        heapFraction:        "Fraction of jobMemory given to the JVM heap"
        jobMemory:           "Memory allocated to the job, in GB"
        cores:               "Number of CPUs allocated to the job"
        timeout:             "Maximum run time, in hours"
        modules:             "Environment modules to load"
    }

    Int heapMb = floor(jobMemory * 1024 * heapFraction)
    String tumor_ext = sub(basename(tumor_bam), "^.*\\.", "")
    String tumor_idx_ext = sub(basename(tumor_bai), "^.*\\.", "")
    String normal_ext = sub(basename(normal_bam), "^.*\\.", "")
    String normal_idx_ext = sub(basename(normal_bai), "^.*\\.", "")

    command <<<
        set -euo pipefail

        # Everything the task will open is bound by its filesystem root: the task directory
        # itself, whatever the caller asked for, and wherever the engine's input links
        # actually resolve to. Roots rather than the directories themselves, because binding a
        # path nested inside another bind makes the container runtime skip it, and one root
        # covers every file beneath it anyway.
        bind_args=(--pwd "$(pwd)")
        bind_seen=()
        add_bind_root() {
            local p real root
            p="${1:-}"
            [ -n "${p}" ] || return 0
            real=$(readlink -f "${p}") || return 0
            [ -e "${real}" ] || return 0
            root="/$(printf '%s' "${real#/}" | cut -d/ -f1)"
            [ "${root}" != "/" ] || return 0
            case " ${bind_seen[*]:-} " in *" ${root} "*) return 0 ;; esac
            bind_seen+=("${root}")
            bind_args+=(-B "${root}")
        }

        add_bind_root "$(pwd)"
        while IFS= read -r b; do add_bind_root "${b}"; done < ~{write_lines(container_binds)}

        # SAGE takes the recalibration and jitter tables from the directory holding the
        # alignment, so alignments, indexes and tables all go into one directory. The names
        # are sample-prefixed, so the two samples' tables do not collide.
        ln -s "~{tumor_bam}" "~{tumor_id}.redux.~{tumor_ext}"
        ln -s "~{tumor_bai}" "~{tumor_id}.redux.~{tumor_ext}.~{tumor_idx_ext}"
        ln -s "~{normal_bam}" "~{normal_id}.redux.~{normal_ext}"
        ln -s "~{normal_bai}" "~{normal_id}.redux.~{normal_ext}.~{normal_idx_ext}"
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -sf "${f}" .; done < ~{write_lines(tumor_tsvs)}
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -sf "${f}" .; done < ~{write_lines(normal_tsvs)}

        mkdir -p somatic

        # The tumour-in-normal check is Illumina-only, and pulls in the panel of normals and
        # the population frequencies with it.
        if [ "~{platform}" = "ILLUMINA" ]; then
            tinc_args="-run_tinc -write_fit_variants -pon_file ~{sage_pon} -gnomad_freq_dir ~{gnomad_dir}"
        else
            tinc_args=""
        fi

        cat > sage.sh <<'COMMAND'
        set -euo pipefail
        sage \
            -Xmx~{heapMb}m \
            -reference ~{normal_id} \
            -reference_bam ~{normal_id}.redux.~{normal_ext} \
            -ref_sample_count 1 \
            -tumor ~{tumor_id} \
            -tumor_bam ~{tumor_id}.redux.~{tumor_ext} \
            -ref_genome ~{genome_fasta} \
            -ref_genome_version ~{genome_version} \
            -hotspots ~{hotspots} \
            -driver_gene_panel ~{driver_gene_panel} \
            -high_confidence_bed ~{high_confidence_bed} \
            -ensembl_data_dir ~{ensembl_data_dir} \
            -sequencing_type ~{platform} \
            -include_mt \
            TINC_ARGS \
            -threads ~{cores} \
            -log_level ~{log_level} \
            -output_vcf somatic/~{tumor_id}.sage.somatic.vcf.gz
COMMAND

        sed -i "s|TINC_ARGS|${tinc_args}|" sage.sh

        while IFS= read -r link; do add_bind_root "${link}"; done < <(find . -type l)

        apptainer exec "${bind_args[@]}" \
            "~{images_dir}/~{image}" bash sage.sh
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        File somatic_vcf = "somatic/~{tumor_id}.sage.somatic.vcf.gz"
        File somatic_tbi = "somatic/~{tumor_id}.sage.somatic.vcf.gz.tbi"
        Array[File] sage_files = glob("somatic/*")
    }
}

# PAVE annotates the SAGE calls with gene consequence, population frequency and panel-of-
# normals membership. PURPLE reads its output rather than SAGE's.
task pave_somatic {
    input {
        String tumor_id
        File sage_vcf
        File sage_tbi
        String platform
        String genome_fasta
        String genome_version
        String sage_pon
        String gnomad_dir
        String clinvar_vcf
        String driver_gene_panel
        String mappability_bed
        String ensembl_data_dir
        String log_level
        String images_dir
        Array[String] container_binds
        String image = "hmftools-pave-1.9--hdfd78af_0.img"
        Float heapFraction = 0.75
        Int jobMemory = 32
        Int cores = 6
        Int timeout = 12
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        tumor_id:          "Tumour sample id, which prefixes the output VCF"
        sage_vcf:          "Somatic VCF from SAGE"
        sage_tbi:          "Index for the SAGE VCF, which the tool resolves beside the VCF"
        platform:          "ILLUMINA, ULTIMA or SBX"
        genome_fasta:      "Reference genome FASTA"
        genome_version:    "Reference build, 37 or 38"
        sage_pon:          "Panel of normals for the sequencing platform in use"
        gnomad_dir:        "Population frequencies"
        clinvar_vcf:       "ClinVar annotations"
        driver_gene_panel: "Driver gene panel"
        mappability_bed:   "Per-region mappability, one of the filters that decides which sites can carry MRD signal"
        ensembl_data_dir:  "Ensembl gene and transcript tables"
        log_level:         "Log level passed to the tool"
        images_dir:        "Directory holding the container images"
        container_binds: "Host paths to bind into the container. Each is reduced to its filesystem root, so naming a directory below one already bound is harmless"
        image:             "Container image filename within images_dir"
        heapFraction:      "Fraction of jobMemory given to the JVM heap"
        jobMemory:         "Memory allocated to the job, in GB"
        cores:             "Number of CPUs allocated to the job"
        timeout:           "Maximum run time, in hours"
        modules:           "Environment modules to load"
    }

    Int heapMb = floor(jobMemory * 1024 * heapFraction)

    command <<<
        set -euo pipefail

        # Everything the task will open is bound by its filesystem root: the task directory
        # itself, whatever the caller asked for, and wherever the engine's input links
        # actually resolve to. Roots rather than the directories themselves, because binding a
        # path nested inside another bind makes the container runtime skip it, and one root
        # covers every file beneath it anyway.
        bind_args=(--pwd "$(pwd)")
        bind_seen=()
        add_bind_root() {
            local p real root
            p="${1:-}"
            [ -n "${p}" ] || return 0
            real=$(readlink -f "${p}") || return 0
            [ -e "${real}" ] || return 0
            root="/$(printf '%s' "${real#/}" | cut -d/ -f1)"
            [ "${root}" != "/" ] || return 0
            case " ${bind_seen[*]:-} " in *" ${root} "*) return 0 ;; esac
            bind_seen+=("${root}")
            bind_args+=(-B "${root}")
        }

        add_bind_root "$(pwd)"
        while IFS= read -r b; do add_bind_root "${b}"; done < ~{write_lines(container_binds)}

        # The index is resolved beside the VCF, so both go into one directory of their own
        # rather than into the task directory, where the localized copy already sits.
        mkdir -p sage_input pave_somatic
        ln -s "~{sage_vcf}" sage_input/
        ln -s "~{sage_tbi}" sage_input/

        cat > pave.sh <<'COMMAND'
        set -euo pipefail
        pave \
            -Xmx~{heapMb}m \
            -sample ~{tumor_id} \
            -input_vcf sage_input/$(basename "~{sage_vcf}") \
            -ref_genome ~{genome_fasta} \
            -ref_genome_version ~{genome_version} \
            -pon_file ~{sage_pon} \
            -gnomad_freq_dir ~{gnomad_dir} \
            -clinvar_vcf ~{clinvar_vcf} \
            -driver_gene_panel ~{driver_gene_panel} \
            -mappability_bed ~{mappability_bed} \
            -ensembl_data_dir ~{ensembl_data_dir} \
            -sequencing_type ~{platform} \
            -threads ~{cores} \
            -log_level ~{log_level} \
            -output_dir pave_somatic/ \
            -output_vcf pave_somatic/~{tumor_id}.pave.somatic.vcf.gz
COMMAND

        while IFS= read -r link; do add_bind_root "${link}"; done < <(find . -type l)

        apptainer exec "${bind_args[@]}" \
            "~{images_dir}/~{image}" bash pave.sh
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        File pave_vcf = "pave_somatic/~{tumor_id}.pave.somatic.vcf.gz"
        File pave_tbi = "pave_somatic/~{tumor_id}.pave.somatic.vcf.gz.tbi"
    }
}

# PURPLE fits purity and ploidy from the AMBER and COBALT signals and rewrites the variant
# calls with copy number and clonality. WISP reads its purity fit and its somatic VCF, so
# this task cannot be skipped even though copy number is out of scope.
task purple {
    input {
        String tumor_id
        String? normal_id
        Array[File] amber_files
        Array[File] cobalt_files
        File pave_vcf
        File pave_tbi
        Array[File] redux_tumor_tsvs
        String genome_fasta
        String genome_version
        String gc_profile
        String hotspots_somatic
        String hotspots_germline
        String driver_gene_panel
        String ensembl_data_dir
        String germline_amp_del_freq
        String log_level
        String images_dir
        Array[String] container_binds
        String image = "hmftools-purple-4.4--hdfd78af_0.img"
        Float heapFraction = 0.75
        Int jobMemory = 32
        Int cores = 6
        Int timeout = 24
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        tumor_id:              "Tumour sample id, which prefixes every output file"
        normal_id:             "Matched normal sample id"
        amber_files:           "AMBER output, reassembled into a directory"
        cobalt_files:          "COBALT output, reassembled into a directory"
        pave_vcf:              "Annotated somatic VCF from PAVE"
        pave_tbi:              "Index for the PAVE VCF"
        redux_tumor_tsvs:      "Tumour REDUX tables, from which PURPLE reads the microsatellite jitter fit"
        genome_fasta:          "Reference genome FASTA"
        genome_version:        "Reference build, 37 or 38"
        gc_profile:            "GC content per window"
        hotspots_somatic:      "Known somatic hotspots"
        hotspots_germline:     "Known germline hotspots"
        driver_gene_panel:     "Driver gene panel"
        ensembl_data_dir:      "Ensembl gene and transcript tables"
        germline_amp_del_freq: "Cohort germline amplification and deletion frequencies"
        log_level:             "Log level passed to the tool"
        images_dir:            "Directory holding the container images"
        container_binds:     "Host paths to bind into the container"
        image:                 "Container image filename within images_dir"
        heapFraction:          "Fraction of jobMemory given to the JVM heap"
        jobMemory:             "Memory allocated to the job, in GB"
        cores:                 "Number of CPUs allocated to the job"
        timeout:               "Maximum run time, in hours"
        modules:               "Environment modules to load"
    }

    Int heapMb = floor(jobMemory * 1024 * heapFraction)

    command <<<
        set -euo pipefail

        # Everything the task will open is bound by its filesystem root: the task directory
        # itself, whatever the caller asked for, and wherever the engine's input links
        # actually resolve to. Roots rather than the directories themselves, because binding a
        # path nested inside another bind makes the container runtime skip it, and one root
        # covers every file beneath it anyway.
        bind_args=(--pwd "$(pwd)")
        bind_seen=()
        add_bind_root() {
            local p real root
            p="${1:-}"
            [ -n "${p}" ] || return 0
            real=$(readlink -f "${p}") || return 0
            [ -e "${real}" ] || return 0
            root="/$(printf '%s' "${real#/}" | cut -d/ -f1)"
            [ "${root}" != "/" ] || return 0
            case " ${bind_seen[*]:-} " in *" ${root} "*) return 0 ;; esac
            bind_seen+=("${root}")
            bind_args+=(-B "${root}")
        }

        add_bind_root "$(pwd)"
        while IFS= read -r b; do add_bind_root "${b}"; done < ~{write_lines(container_binds)}

        # PURPLE takes directories, so the upstream outputs are relinked into the shapes it
        # expects.
        mkdir -p amber cobalt pave_somatic redux_tumor purple
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -s "${f}" amber/; done < ~{write_lines(amber_files)}
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -s "${f}" cobalt/; done < ~{write_lines(cobalt_files)}
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -s "${f}" redux_tumor/; done < ~{write_lines(redux_tumor_tsvs)}
        ln -s "~{pave_vcf}" pave_somatic/
        ln -s "~{pave_tbi}" pave_somatic/

        # No structural variants are called, so PURPLE is not given an ESVEE directory and
        # fits copy number from the read ratios alone.
        cat > purple.sh <<'COMMAND'
        set -euo pipefail
        purple \
            -Xmx~{heapMb}m \
            -tumor ~{tumor_id} \
            ~{if defined(normal_id) then "-reference " + normal_id else ""} \
            -amber amber/ \
            -cobalt cobalt/ \
            -pave_somatic_dir pave_somatic/ \
            -redux_tumor_dir redux_tumor/ \
            -ref_genome ~{genome_fasta} \
            -ref_genome_version ~{genome_version} \
            -driver_gene_panel ~{driver_gene_panel} \
            -ensembl_data_dir ~{ensembl_data_dir} \
            -somatic_hotspots ~{hotspots_somatic} \
            -germline_hotspots ~{hotspots_germline} \
            -germline_amp_del_freq_file ~{germline_amp_del_freq} \
            -gc_profile ~{gc_profile} \
            -circos $(which circos) \
            -threads ~{cores} \
            -log_level ~{log_level} \
            -output_dir purple/
COMMAND

        while IFS= read -r link; do add_bind_root "${link}"; done < <(find . -type l)

        apptainer exec "${bind_args[@]}" \
            "~{images_dir}/~{image}" bash purple.sh

        [ -f "purple/~{tumor_id}.purple.somatic.vcf.gz" ]
        [ -f "purple/~{tumor_id}.purple.purity.tsv" ]

        # Cromwell delocalizes files, not directories, so any subdirectory the tool wrote is
        # archived rather than listed as an output.
        mkdir -p purple_extra
        find purple -mindepth 1 -maxdepth 1 -type d -exec mv {} purple_extra/ \;
        tar -czf ~{tumor_id}.purple_plots.tar.gz -C purple_extra .
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        File somatic_vcf = "purple/~{tumor_id}.purple.somatic.vcf.gz"
        File somatic_tbi = "purple/~{tumor_id}.purple.somatic.vcf.gz.tbi"
        File purity_tsv = "purple/~{tumor_id}.purple.purity.tsv"
        Array[File] purple_files = glob("purple/*")
        File purple_plots = "~{tumor_id}.purple_plots.tar.gz"
    }
}

# Packs the primary-stage results a later longitudinal run needs, so that one primary can be
# computed once and measured against many timepoints.
task pack_primary {
    input {
        String tumor_id
        Array[File] purple_files
        Array[File] amber_files
        File purple_plots
        File amber_plots
        File cobalt_plots
        Int jobMemory = 8
        Int cores = 1
        Int timeout = 4
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        tumor_id:     "Primary tumour sample id, which names the archive"
        purple_plots: "PURPLE plot directories, archived by the task that made them"
        amber_plots:  "AMBER plot directories, archived by the task that made them"
        cobalt_plots: "COBALT plot directories, archived by the task that made them"
        purple_files: "PURPLE output"
        amber_files:  "AMBER output, carried so that a later run can add LOH evidence without recomputing the primary"
        jobMemory:    "Memory allocated to the job, in GB"
        cores:        "Number of CPUs allocated to the job"
        timeout:      "Maximum run time, in hours"
        modules:      "Environment modules to load"
    }

    command <<<
        set -euo pipefail

        mkdir -p primary/purple primary/amber primary/plots
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -s "${f}" primary/purple/; done < ~{write_lines(purple_files)}
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -s "${f}" primary/amber/; done < ~{write_lines(amber_files)}

        # The tools write plots into subdirectories, which the tasks archive so that their
        # outputs are files. Carried here so the plots are not simply discarded.
        for plots in "~{purple_plots}" "~{amber_plots}" "~{cobalt_plots}"; do
            [ -s "${plots}" ] || continue
            ln -s "${plots}" primary/plots/
        done

        # Dereferenced, so the archive holds the files rather than links into a task
        # directory that will be cleaned up.
        tar -czhf ~{tumor_id}.primary.tar.gz primary
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        File tarball = "~{tumor_id}.primary.tar.gz"
    }
}

# Unpacks a primary tarball and reads the primary tumour sample id back out of the PURPLE
# filenames, so a longitudinal run needs no separate declaration of it.
task extract_primary {
    input {
        File tarball
        Int jobMemory = 8
        Int cores = 1
        Int timeout = 4
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        tarball:   "Primary-stage archive from an earlier run"
        jobMemory: "Memory allocated to the job, in GB"
        cores:     "Number of CPUs allocated to the job"
        timeout:   "Maximum run time, in hours"
        modules:   "Environment modules to load"
    }

    command <<<
        set -euo pipefail

        mkdir -p primary
        tar -xzf "~{tarball}" -C primary --strip-components=1

        purity=$(find primary/purple -name '*.purple.purity.tsv' | head -1)
        if [ -z "${purity}" ]; then
            echo "ERROR: no *.purple.purity.tsv in the primary archive; it does not look like" \
                 "a primary-stage output" >&2
            exit 1
        fi
        basename "${purity}" .purple.purity.tsv > tumor_id.txt

        # What the primary was called with. The somatic VCF carries the sequence dictionary
        # the caller used, which is the one thing in the archive that can show the primary was
        # called against a different reference from the one this run uses.
        sed -n 's/^version=//p' primary/purple/purple.version 2>/dev/null | head -1 > purple_version.txt
        [ -s purple_version.txt ] || echo unknown > purple_version.txt

        vcf="primary/purple/$(cat tumor_id.txt).purple.somatic.vcf.gz"
        : > primary_contigs.txt
        if [ -f "${vcf}" ]; then
            zcat "${vcf}" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' > primary_contigs.txt
        fi

        echo "primary tumour sample id: $(cat tumor_id.txt)," \
             "called by PURPLE $(cat purple_version.txt)," \
             "$(grep -c . primary_contigs.txt) contigs in its call set" >&2
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        String tumor_id = read_string("tumor_id.txt")
        String purple_version = read_string("purple_version.txt")
        File contig_names = "primary_contigs.txt"
        Array[File] purple_files = glob("primary/purple/*")
        Array[File] amber_files = glob("primary/amber/*")
    }
}

# SAGE append force-calls the primary's somatic sites in the longitudinal sample. No
# discovery happens: the site list is fixed by the primary's VCF, which is what makes the
# depths comparable across timepoints.
task sage_append {
    input {
        String primary_id
        String longitudinal_id
        Array[File] purple_files
        File longitudinal_bam
        File longitudinal_bai
        Array[File] longitudinal_tsvs
        String platform
        String genome_fasta
        String genome_version
        String outputFileNamePrefix
        String log_level
        String images_dir
        Array[String] container_binds
        String image = "hmftools-sage-5.0.2--hdfd78af_0.img"
        Float heapFraction = 0.75
        Int jobMemory = 48
        Int cores = 8
        Int timeout = 48
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        primary_id:         "Primary tumour sample id, which names the PURPLE VCF the sites come from"
        longitudinal_id:    "Longitudinal sample id, which prefixes the output VCF"
        purple_files:       "Primary PURPLE output, from which the somatic VCF and its index are taken"
        longitudinal_bam:   "Longitudinal REDUX alignment"
        longitudinal_bai:   "Index for the longitudinal alignment"
        longitudinal_tsvs:  "Longitudinal REDUX tables. Their absence would silently switch the tool to skipping recalibration and jitter fitting, changing the depths it reports, so they are required here"
        platform:           "ILLUMINA, ULTIMA or SBX"
        genome_fasta:       "Reference genome FASTA"
        genome_version:     "Reference build, 37 or 38"
        outputFileNamePrefix: "Prefix for the output VCF, which is provisioned"
        log_level:          "Log level passed to the tool"
        images_dir:         "Directory holding the container images"
        container_binds:  "Host paths to bind into the container"
        image:              "Container image filename within images_dir"
        heapFraction:       "Fraction of jobMemory given to the JVM heap"
        jobMemory:          "Memory allocated to the job, in GB"
        cores:              "Number of CPUs allocated to the job"
        timeout:            "Maximum run time, in hours"
        modules:            "Environment modules to load"
    }

    Int heapMb = floor(jobMemory * 1024 * heapFraction)
    String long_ext = sub(basename(longitudinal_bam), "^.*\\.", "")
    String long_idx_ext = sub(basename(longitudinal_bai), "^.*\\.", "")

    command <<<
        set -euo pipefail

        # Everything the task will open is bound by its filesystem root: the task directory
        # itself, whatever the caller asked for, and wherever the engine's input links
        # actually resolve to. Roots rather than the directories themselves, because binding a
        # path nested inside another bind makes the container runtime skip it, and one root
        # covers every file beneath it anyway.
        bind_args=(--pwd "$(pwd)")
        bind_seen=()
        add_bind_root() {
            local p real root
            p="${1:-}"
            [ -n "${p}" ] || return 0
            real=$(readlink -f "${p}") || return 0
            [ -e "${real}" ] || return 0
            root="/$(printf '%s' "${real#/}" | cut -d/ -f1)"
            [ "${root}" != "/" ] || return 0
            case " ${bind_seen[*]:-} " in *" ${root} "*) return 0 ;; esac
            bind_seen+=("${root}")
            bind_args+=(-B "${root}")
        }

        add_bind_root "$(pwd)"
        while IFS= read -r b; do add_bind_root "${b}"; done < ~{write_lines(container_binds)}

        # The index has to sit beside the VCF, and the recalibration and jitter tables beside
        # the alignment.
        mkdir -p purple_primary
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -s "${f}" purple_primary/; done < ~{write_lines(purple_files)}
        ln -s "~{longitudinal_bam}" "~{longitudinal_id}.redux.~{long_ext}"
        ln -s "~{longitudinal_bai}" "~{longitudinal_id}.redux.~{long_ext}.~{long_idx_ext}"
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -sf "${f}" .; done < ~{write_lines(longitudinal_tsvs)}

        for required in redux.bqr.tsv redux.jitter_params.tsv redux.ms_table.tsv.gz; do
            if [ ! -e "~{longitudinal_id}.${required}" ]; then
                echo "ERROR: ~{longitudinal_id}.${required} is missing; without it the tool" \
                     "skips recalibration and jitter fitting and reports different depths" >&2
                exit 1
            fi
        done

        mkdir -p sage_append

        cat > sage_append.sh <<'COMMAND'
        set -euo pipefail
        sage \
            -Xmx~{heapMb}m \
            com.hartwig.hmftools.sage.append.SageAppendApplication \
            -input_vcf purple_primary/~{primary_id}.purple.somatic.vcf.gz \
            -max_read_depth 100000 \
            -reference ~{longitudinal_id} \
            -reference_bam ~{longitudinal_id}.redux.~{long_ext} \
            -ref_genome ~{genome_fasta} \
            -ref_genome_version ~{genome_version} \
            -sequencing_type ~{platform} \
            -write_frag_lengths \
            -threads ~{cores} \
            -log_level ~{log_level} \
            -output_vcf sage_append/~{outputFileNamePrefix}.sage.append.vcf.gz
COMMAND

        while IFS= read -r link; do add_bind_root "${link}"; done < <(find . -type l)

        apptainer exec "${bind_args[@]}" \
            "~{images_dir}/~{image}" bash sage_append.sh
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        File append_vcf = "sage_append/~{outputFileNamePrefix}.sage.append.vcf.gz"
        File append_tbi = "sage_append/~{outputFileNamePrefix}.sage.append.vcf.gz.tbi"
        Array[File] sage_append_files = glob("sage_append/*")
    }
}

# WISP estimates what fraction of the longitudinal sample is tumour, from the depths SAGE
# append measured at the primary's somatic sites.
task wisp_purity {
    input {
        String donor_id
        String primary_id
        String longitudinal_id
        Array[File] purple_files
        File append_vcf
        File append_tbi
        Array[File] longitudinal_tsvs
        Array[File]? cobalt_files
        Boolean use_copy_number
        String genome_fasta
        String outputFileNamePrefix
        String log_level
        String images_dir
        Array[String] container_binds
        String image = "hmftools-wisp-1.3.1--hdfd78af_0.img"
        Float heapFraction = 0.75
        Int jobMemory = 32
        Int cores = 2
        Int timeout = 12
        String modules = "wisp/3.0.0"
    }

    parameter_meta {
        donor_id:          "Donor the samples came from, reported as patient_id"
        primary_id:        "Primary tumour sample id, whose PURPLE fit supplies the expected variant copy numbers"
        longitudinal_id:   "Longitudinal sample id, the sample being measured"
        purple_files:      "Primary PURPLE output, reassembled into a directory"
        append_vcf:        "The primary's somatic sites as force-called in the longitudinal sample"
        append_tbi:        "Index for that VCF"
        longitudinal_tsvs: "Longitudinal REDUX tables, from which the per-base error rates are read"
        cobalt_files:      "Longitudinal COBALT output. MANDATORY when use_copy_number is set and unused otherwise"
        use_copy_number:   "Whether COPY_NUMBER joins SOMATIC_VARIANT in the purity methods applied"
        genome_fasta:      "Reference genome FASTA"
        outputFileNamePrefix: "Prefix for the provisioned summary and archive"
        log_level:         "Log level passed to the tool"
        images_dir:        "Directory holding the container images"
        container_binds: "Host paths to bind into the container. Each is reduced to its filesystem root, so naming a directory below one already bound is harmless"
        image:             "Container image filename within images_dir"
        heapFraction:      "Fraction of jobMemory given to the JVM heap"
        jobMemory:         "Memory allocated to the job, in GB"
        cores:             "Number of CPUs allocated to the job"
        timeout:           "Maximum run time, in hours"
        modules:           "Environment modules to load"
    }

    Int heapMb = floor(jobMemory * 1024 * heapFraction)

    command <<<
        set -euo pipefail

        # Everything the task will open is bound by its filesystem root: the task directory
        # itself, whatever the caller asked for, and wherever the engine's input links
        # actually resolve to. Roots rather than the directories themselves, because binding a
        # path nested inside another bind makes the container runtime skip it, and one root
        # covers every file beneath it anyway.
        bind_args=(--pwd "$(pwd)")
        bind_seen=()
        add_bind_root() {
            local p real root
            p="${1:-}"
            [ -n "${p}" ] || return 0
            real=$(readlink -f "${p}") || return 0
            [ -e "${real}" ] || return 0
            root="/$(printf '%s' "${real#/}" | cut -d/ -f1)"
            [ "${root}" != "/" ] || return 0
            case " ${bind_seen[*]:-} " in *" ${root} "*) return 0 ;; esac
            bind_seen+=("${root}")
            bind_args+=(-B "${root}")
        }

        add_bind_root "$(pwd)"
        while IFS= read -r b; do add_bind_root "${b}"; done < ~{write_lines(container_binds)}

        mkdir -p purple_primary sage_append_longitudinal redux_longitudinal cobalt_longitudinal wisp
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -s "${f}" purple_primary/; done < ~{write_lines(purple_files)}
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -s "${f}" redux_longitudinal/; done < ~{write_lines(longitudinal_tsvs)}
        ln -s "~{append_vcf}" sage_append_longitudinal/
        ln -s "~{append_tbi}" sage_append_longitudinal/
        while IFS= read -r f; do [ -n "${f}" ] || continue; ln -s "${f}" cobalt_longitudinal/; done < ~{write_lines(select_first([cobalt_files, []]))}

        # AMBER_LOH is not offered: it needs the primary AMBER directory together with the
        # primary normal alignment, which the longitudinal stage does not carry.
        if ~{if use_copy_number then "true" else "false"}; then
            methods='SOMATIC_VARIANT;COPY_NUMBER'
            extra_args="-cobalt_dir cobalt_longitudinal/ -write_types ALL"
        else
            methods='SOMATIC_VARIANT'
            extra_args="-write_types 'SOMATIC_DATA;SOMATIC_PLOT'"
        fi

        cat > wisp.sh <<'COMMAND'
        set -euo pipefail
        wisp \
            -Xmx~{heapMb}m \
            com.hartwig.hmftools.wisp.purity.PurityEstimator \
            -patient_id ~{donor_id} \
            -tumor_id ~{primary_id} \
            -samples ~{longitudinal_id} \
            -purity_methods 'METHODS' \
            -somatic_vcf sage_append_longitudinal/~{basename(append_vcf)} \
            -purple_dir purple_primary/ \
            -bqr_dir redux_longitudinal/ \
            EXTRA_ARGS \
            -ref_genome ~{genome_fasta} \
            -log_level ~{log_level} \
            -output_dir wisp/
COMMAND

        sed -i "s|METHODS|${methods}|; s|EXTRA_ARGS|${extra_args}|" wisp.sh

        while IFS= read -r link; do add_bind_root "${link}"; done < <(find . -type l)

        apptainer exec "${bind_args[@]}" \
            "~{images_dir}/~{image}" bash wisp.sh

        cp wisp/~{donor_id}_~{longitudinal_id}.wisp.summary.tsv \
           ~{outputFileNamePrefix}.wisp.summary.tsv
        tar -czhf ~{outputFileNamePrefix}.wisp.tar.gz wisp
    >>>

    runtime {
        memory: "~{jobMemory} GB"
        cpu: "~{cores}"
        timeout: "~{timeout}"
        modules: "~{modules}"
    }

    output {
        File summary = "~{outputFileNamePrefix}.wisp.summary.tsv"
        File tarball = "~{outputFileNamePrefix}.wisp.tar.gz"
        Array[File] wisp_files = glob("wisp/*")
    }
}
