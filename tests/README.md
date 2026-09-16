# Test data

The regression cases read from `/.mounts/labs/gsi/testdata/wisp/test_data/`.

Carried over from the previous workflow, already present:

    test_01_Nk_P_WG_Novaseq_X_Validation_v1-2.downsampled.bam   primary tumour
    test_01_Ly_R_WG_Novaseq_X_Validation_v1-2.downsampled.bam   matched normal
    test_01_Ct_T_PG_Novaseq_X_validation_v1_2.downsampled.bam   longitudinal

Two more have to be staged by hand, because a case cannot depend on another case and the
`PE` cases need a primary that already exists:

    test_01_Nk_P.primary.tar.gz   the primary_output of a WG run over the tumour and normal
    test_01_Ct_T.redux/           a REDUX output directory for the longitudinal sample

Both come out of one `WG_PE` run over the three alignments above. Take `primary_output` for
the tarball, and the `redux/` directory from that run's `call-redux_longitudinal/execution`
for the REDUX directory. The REDUX directory must hold `{sample_id}.redux.bam` with its
index and the bqr, jitter_params and ms_table tables, all sharing a prefix equal to the
alignment's read-group SM tag.

## What each case covers

| case | exercises |
|---|---|
| `test_01_wg` | the primary stage alone, and the tarball a later run consumes |
| `test_02_pe` | `extract_primary`, and a longitudinal sample from alignments |
| `test_03_pe_redux_dir_snv_only` | the REDUX-directory route, and `use_copy_number = false` |
| `test_04_wg_pe` | both stages in one run |

The cases differ in output shape rather than in input values, so each one catches something
the others cannot. `test_03` is the exception: it produces the same file list as `test_02`
but reaches it by a different route and with a shorter WISP summary, which
`tests/calculate.sh` records.

## Before the first run

Every alignment must list contigs in the same order as the reference the module ships, or
`validate_inputs` rejects it. Check with:

    diff <(samtools view -H <aln> | awk '$1=="@SQ"{for(i=2;i<=NF;i++) if($i~/^SN:/) print substr($i,4)}') \
         <(cut -f1 "$WISP_GENOME_FASTA.fai" | head -n "$(samtools view -H <aln> | grep -c '^@SQ')")

Empty output means it passes.
