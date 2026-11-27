# oncoanalyser with bwa-mem3 and the Rust ports of hmftools

This fork of nf-core/oncoanalyser (branch `rust-ports`) runs the pipeline with:

- **bwa-mem3** for DNA read alignment (RNA: STAR by default, see "RNA alignment" below). The branch
  starts from upstream's `bwa-mem3` branch (Hartwig's development line for oncoanalyser 3.1), which
  already replaces BWA-MEM2 with bwa-mem3. Upstream runs it with `--compat=bwa-mem2`, which writes the
  same records as BWA-MEM2. This fork runs bwa-mem3's own mode by default (`--bwamem3_compat none`;
  `bwa-mem2` gives upstream's records), which bwa-mem3's authors recommend for new pipelines: on a
  whole-genome lane it gave the same alignment as compat mode in every record, in the same time, and
  adds the `MQ` (mate mapping quality) and `HN` tags (Alignment, below).
- **No fastp.** fastp only split FASTQ files into chunks for parallel alignment (or extracted UMIs);
  it changed no reads. `max_fastq_records` defaults to `-1`, so FASTQs are aligned whole and fastp
  does not run unless UMIs are enabled.
- **Rust ports of SAGE, REDUX, ESVEE, AMBER, COBALT, the BamTools metrics, TEAL, LILAC, CIDER, PAVE,
  LINX, SIGS, Virus Interpreter, PEACH, CHORD, NEO, CUPPA and qsee** (and PURPLE, see below) in place of the Java tools, chosen with a
  profile. The ports live in the hmftools fork
  (https://github.com/slowkow/hmftools, branch `rust-ports`): `sage-rs/`, `redux-rs/`, `esvee-rs/`,
  `amber-rs/`, `cobalt-rs/`, `bam-tools-rs/`, `teal-rs/`, `lilac-rs/`, `cider-rs/`, `pave-rs/`, `purple-rs/`, `linx-rs/`, `sigs-rs/`, `virus-interpreter-rs/`, `peach-rs/`, `chord-rs/`, `neo-rs/`, `cuppa-rs/`, `qsee-rs/`. Each
  takes the Java
  tool's options, writes the same files, and on every test input gives Java's output exactly when asked
  to (see below).

## Running

```bash
nextflow run slowkow/oncoanalyser -r rust-ports \
    -profile docker,rust_ports \
    --mode wgts --genome GRCh38_hmf \
    --input samplesheet.csv --outdir output/
```

| Profile | Effect |
|---|---|
| `rust_ports` | SAGE, REDUX, ESVEE, AMBER, COBALT, BamTools, TEAL, LILAC, CIDER, PAVE, PURPLE (with its plots), LINX (LINX_SOMATIC, LINX_GERMLINE, LINX_VISUALISER), SIGS, Virus Interpreter, PEACH, CHORD, NEO (NEO_FINDER, NEO_SCORER), CUPPA, qsee and ORANGE all run as Rust ports |
| `rust_sage`, `rust_redux`, `rust_esvee`, `rust_amber`, `rust_cobalt`, `rust_bamtools`, `rust_teal`, `rust_lilac`, `rust_cider`, `rust_pave`, `rust_purple`, `rust_linx`, `rust_sigs`, `rust_virusinterpreter`, `rust_peach`, `rust_chord`, `rust_neo`, `rust_cuppa`, `rust_qsee`, `rust_orange` | one tool as a Rust port, the others in Java; combine as needed |

| Parameter | Default | Meaning |
|---|---|---|
| `--rust_ports_container` | `hmftools-rust:dev` | image with the ports (built locally, see below) |
| `--rust_ports_bin` | not set | folder with the ports on the host; these steps then run without a container (see below) |
| `--rust_java_compat` | `false` | reproduce the Java tools' output exactly, including the faults the ports correct by default |
| `--rust_low_base_qual_threshold` | not set (26) | SAGE and ESVEE: bases at or below this quality are low quality |
| `--rust_high_base_qual_threshold` | not set (30) | SAGE and ESVEE: Illumina bases at or above this quality are high quality |
| `--publish_alignments` | `true` | `false` leaves REDUX's BAMs and side files (`alignments/`) out of the output folder: for a rerun whose REDUX tasks come from the `-resume` cache (the pipeline publishes through its `output` block, so a process `publishDir` setting does not stop it) |

Give `--rust_ports_bin` and `--rust_java_compat` on the command line or in a `-params-file`: they are
read while the configuration is loaded, so a value set in a `-c` config file is not seen.
`--rust_ports_container` can be set in either place (also in a `-c` file). It is read as each task
starts, and Nextflow then warns once at start-up, before the pipeline's parameters are loaded:
`Access to undefined parameter rust_ports_container`. That warning is harmless.

By default the ports correct the faults found in the Java code (each is listed with an id in the
port's `REPORT.md`). `--rust_java_compat true` adds `-java_compat` to the command line of each ported
step, which makes each port behave as Java does. Until 2026-11-27 it set `<TOOL>_JAVA_COMPAT` through
Nextflow's `env` scope instead, which Google Batch does not pass to tasks, so on Batch it had no effect.

The two base-quality parameters are options of the ports, not of the Java tools, so they apply only to
the steps that run as ports. They are for data whose base qualities were rewritten by recalibration and
stop below Q30: Java SAGE then passes no SNV, since it counts only bases at Q30 or above towards
`minTumorSupport`, and Java ESVEE treats most of the tumor's bases as low quality, so the tumor
supports few SVs. Each port's README says what the bands do. When neither is set the task commands are
unchanged (so is the `-resume` cache); they can be given anywhere, including a `-c` file. Hartwig
tuned its filters on data with the usual qualities, so calls made with moved bands need checking.

The profiles also lower the memory requested for these steps (SAGE, REDUX, AMBER and COBALT 8 GB,
BamTools, TEAL and LILAC 4 GB, CIDER 8 GB (for the pages of the BWA index image it maps), PAVE and LINX 2 GB, SIGS, Virus Interpreter,
PEACH, CHORD and qsee 1 GB, NEO and CUPPA 2 GB, ESVEE 24 GB,
each doubled on retry), since the ports need a few GB where the Java tools are given 24 to 72 GB. More
samples then fit on one machine at once.

## RNA alignment: STAR (default) or bwa-mem3 + TARS

The RNA route of upstream's `bwa-mem3` branch aligns RNA with bwa-mem3 to the genome with transcript
contigs appended, lifts the alignments back to genome coordinates with TARS and runs REDUX on them. It
needs reference data Hartwig has not published: that genome and its bwa-mem2 index
(`--ref_data_genome_fasta_rna`, `--ref_data_genome_bwamem2_index_rna`) and the resource files
`contigs_mapping_rna` and `unmap_regions_rna`. So RNA FASTQ input could not run. `--rna_aligner`
chooses the route:

| `--rna_aligner` | Steps | Needs |
|---|---|---|
| `star` (default) | STAR, samtools sort, Picard MarkDuplicates, as nf-core/oncoanalyser's `dev` branch does; REDUX does not run on RNA | `--ref_data_genome_star_index`, a folder or `.tar.gz`; set for `GRCh37_hmf` and `GRCh38_hmf` to Hartwig's published STAR 2.7.3a indexes. The pipeline does not build one |
| `bwa-mem3` | bwa-mem3, TARS, REDUX (the route before this option) | the RNA reference data above |

With `star`, the MarkDuplicates folder (`<sample>.md.bam`, `.md.bam.bai`, `.md.metrics`) takes the place of the RNA
REDUX folder: Isofox, BamTools, SAGE append, LILAC, CIDER and NEO read `<sample>.md.bam`. It is published in
`<group>/alignments/<sample>/` (left out with `--publish_alignments false`, as the REDUX BAMs are), and STAR's
`Log.final.out` in `<group>/star/<sample>.Log.final.out`. STAR aligns all lanes of a sample in one task, so
`max_fastq_records` does not split RNA FASTQs. STAR asks for `process_high` (12 CPUs, 64 GB), samtools sort and
MarkDuplicates for `process_medium` (6 CPUs, 36 GB). A later run cannot start from the published `md.bam` (the
samplesheet's REDUX inputs need REDUX's TSVs), so give RNA as FASTQ.

## Resources

What each step asks the executor for, and what the ports were measured to use on a whole-genome
Illumina tumor/normal pair (from the ports' REPORT.md files; Java's peak depends on the heap it is
given, which oncoanalyser sets to 95% of the task's memory):

| Step | Requested (CPUs, memory) | Rust, measured | Java, measured |
|---|---|---|---|
| SAGE somatic | 12 CPUs, 8 GB | 316 s, 1.7 GB at 16 threads | 3,407 s, 26.9 GB at 16 threads; 70 GB in an oncoanalyser run |
| REDUX, per sample | 12 CPUs, 8 GB | 26 min, 1.6 GB at 16 threads (Linux) | 50 min, 35.7 GB at 16 threads; 51 GB in an oncoanalyser run |
| ESVEE | 12 CPUs, 24 GB | 14 min at 8 threads; peak 7.2 GB (assembly), steps run one at a time | 121 GB in an oncoanalyser run |

CPUs and memory are doubled on each retry (`task.attempt`). The CPU requests are upstream's
(`process_high`); the memory requests are this fork's (upstream asks 72 GB). The one-call ESVEE mode,
which oncoanalyser runs, peaked at 12.9 and 14.3 GB in two runs of a 35x whole genome on Google Batch (section "Whole
genome" below), within the 24 GB requested.

On Google Batch every task gets its own VM, chosen from the CPU and memory requested, so the requests
set the machine type and the cost. On AWS Batch tasks share instances, and the requests decide how
many fit on one. Either way, keep the requests close to measured use.

The trace of each run (`pipeline_info/execution_trace_*.txt`) records, per task, what was requested
(`cpus`, `memory`, `attempt`) and what was used (`peak_rss`, `realtime`, `%cpu`). Summarise one or more
runs per step with:

```bash
tools/trace_summary.py output/pipeline_info/execution_trace_*.txt
```

It prints, per process, the requests, the largest and median peak memory, the largest share of its
request a task used, run times and CPU hours, then lists tasks that came near their memory request
and tasks that failed (exit 137 is usually out of memory).

## Without a container

The ports are single programs with no dependencies beyond the operating system, so they do not need
a container. With `--rust_ports_bin <folder>`, the ported steps run on the host from that folder,
while every other step keeps its container (or conda environment). This is the better choice on Apple
Silicon, where Docker runs the pipeline's x86-64 images under emulation in a virtual machine: the ports
then run as native arm64 programs. Build the folder from the hmftools fork:

```bash
git clone -b rust-ports https://github.com/slowkow/hmftools.git
hmftools/rust-docker/install.sh ~/rust-ports
nextflow run slowkow/oncoanalyser -r rust-ports -profile docker,rust_ports \
    --rust_ports_bin ~/rust-ports ...
```

`samtools` and `sambamba` must also be on the host's `PATH`: oncoanalyser passes them to REDUX and
ESVEE as `-bamtool`, and the ports check that the path exists (as the Java tools do), though they sort
and index BAMs themselves.

Executors that run every task in a container (Google Cloud Batch, AWS Batch, Kubernetes) need the
image instead.

## The container

The image holds `sage-rs`, `redux-rs`, `esvee-rs`, `amber-rs`, `cobalt-rs`, `bam-tools-rs`, `teal-rs`,
`lilac-rs`, `cider-rs`, `pave-rs`, `purple-rs`, `linx-rs`, `sigs-rs`, `virus-interpreter-rs`, `peach-rs`, `chord-rs`,
`neo-rs`, `cuppa-rs` and `qsee-rs`, wrapper commands `sage`, `redux`, `esvee`, `amber`, `cobalt`, `bamtools`, `teal`, `lilac`, `cider`, `pave`, `purple`,
`linx`, `sigs`, `virusinterpreter`, `peach`, `chord`, `neo`, `cuppa` and `qsee` (and a `python` that only runs `python -m
cuppa.predict`) that accept the Java tools' command lines (JVM options are dropped, Java main classes become
subcommands), samtools and sambamba (passed to REDUX and ESVEE as `-bamtool`), and `ps`, which
Nextflow needs to record each task's memory and CPU use. It is built from the hmftools fork, for
linux/amd64 (the machines of Google Batch and AWS Batch). On Google Cloud, Cloud Build builds it and
pushes it to Artifact Registry, with no Docker needed locally:

```bash
git clone -b rust-ports https://github.com/slowkow/hmftools.git
cd hmftools
rust-docker/cloud-build.sh <region>-docker.pkg.dev/<project>/<repository>/hmftools-rust
```

The image is tagged with the hmftools commit id, and the build runs the ports' unit tests on Linux.
For AWS Batch, build the same Dockerfile (`rust-docker/Dockerfile`) and push it to ECR. Then point
the pipeline at it, for example in a private config file given with `-c`:

```groovy
params.rust_ports_container = '<region>-docker.pkg.dev/<project>/<repository>/hmftools-rust:<commit id>'
```

The Batch machines need permission to pull from that registry (on Google Cloud, the Artifact Registry
Reader role for the service account the jobs run as).

## Versions

The Java images in this branch (`conf/hmf_containers.config`) are Hartwig's pre-release builds:
SAGE 5.1-beta.3 and ESVEE 2.1-beta.4; REDUX is the bioconda 2.0.7. The ports reproduce SAGE 5.1 and
REDUX 2.1 as of hmftools commit a224a2cd8e, and ESVEE 2.1 as of commit ca04e6c39d. So a run with the
Java tools and a run with `rust_ports --rust_java_compat true` are not expected to be identical
where those versions differ, REDUX in particular.

AMBER, COBALT and BamTools are closer: the branch's images `amber:4.4-beta.1`, `cobalt:3.1-beta.1` and
`bam-tools:1.7-beta.3` were built from hmftools tags whose sources for these tools (and the hmf-common
code they use) are the same as at the commit the ports follow, so their Java runs and the ports in
Java mode should agree.

TEAL is 1.4 in both the module's image and the port. On an earlier pipeline run's telbams, the port in
Java mode wrote the same tellength and breakend files as that run's Java TEAL. LILAC's image is 2.0, and
the port follows 2.1. On the HG008 pair, LILAC 2.1 gives the same rows in all three files as the
pipeline's LILAC 2.0, in another gene group order. CIDER is 1.2 in both, but the two builds differ:
the port follows the hmftools 1.2 sources, which write two more columns in `cider.vdj.tsv.gz`
(`vImgtRefUnaligned`, `jImgtRefUnaligned`) than the module's image; the other columns and files were the
same on HG008.

## Known gaps

- The modules also record the Java version in `versions.yml` (`java --version`), and COBALT's the R
  and R package versions. The image has no Java or R, so those entries are empty; the step is not
  affected.
- On macOS (`--rust_ports_bin`), the modules read tool versions with GNU `sed` syntax, which the
  system `sed` rejects, so the versions in `versions.yml` are empty unless GNU sed is first on `PATH`
  (for example Homebrew's `gnu-sed`). The steps themselves are not affected.
- redux-rs does not draw the base-quality recalibration plot (`-bqr_write_plot`); no later step
  reads it.
- Run end to end on Google Batch on a 2M-read-pair test and on two whole-genome pairs (below); not
  yet on AWS Batch. AMBER, COBALT and BamTools, ported on 2026-11-26, have run inside the pipeline on
  the HG008 pair only (section "AMBER, COBALT and BamTools in the pipeline"). TEAL, LILAC and CIDER,
  ported on 2026-11-27, have run inside the pipeline on HG008 too (section "TEAL, LILAC and CIDER in the
  pipeline").
- CUPPA and qsee have run inside the pipeline on HG008 (section "CUPPA and qsee in the pipeline"). CUPPA's
  numbers differ from the Linux pipeline's pycuppa in the last digits; neither draws its plot without R.
- ORANGE has run inside the pipeline on HG008 (section "ORANGE in the pipeline"): the JSON only, no PDF; its
  virus key is `virusInterpreter` where the pipeline's newer image writes `viridian`.
- SIGS, Virus Interpreter, PEACH, CHORD and NEO have run inside the pipeline on HG008 (sections "SIGS,
  Virus Interpreter, PEACH and CHORD in the pipeline" and "NEO in the pipeline").
- LINX, PAVE and PURPLE have run inside the pipeline on HG008 (sections "LINX in the pipeline" and
  "PAVE and PURPLE in the pipeline"). PURPLE draws its plots natively since 2026-11-27 and is part of
  `rust_ports` (section "PURPLE's plots and ORANGE in the pipeline").
  LINX_VISUALISER keeps the Java image and reads the port's visualiser files.
- TEAL's prep took 7.5 minutes on Google Batch against Java's 13.8 (HG008), after reading each BAM
  ahead with several handles; before that it was slower than Java (15.5 minutes).

## Google Batch

A 2M-read-pair tumor/normal test from FASTQ (DNA only, `-profile rust_ports`, image built with
`rust-docker/cloud-build.sh`) completed on Google Batch on 2026-11-23: bwa-mem3, the three Rust ports
and every later step. ORANGE did not run: this branch runs it only when PURPLE has drawn its Circos
plot, which PURPLE does not do when it finds no tumor, as at this depth. What the first attempts ran
into, all fixed or worked around:

- **Command files among the outputs.** The modules of this branch list `.command.*` among their
  outputs. On Google Batch each task runs in a local scratch folder in which `.command.run` and
  `.command.sh` are links to the copies in the work folder, so copying the outputs back moves each
  link onto the file it points to: GNU `mv` and `cp` refuse ("are the same file", the task fails), and
  busybox `mv` replaces the file with a link to itself. Configs that shadow `mv` with `cp -fR` (a
  workaround for older oncoanalyser releases, which moved a `dir/` output before `dir/file` outputs;
  no module here has such a pair) fail on every task whose image has GNU coreutils. A `beforeScript`
  that skips such links works:

  ```groovy
  process.beforeScript = 'mv() { local dst="${@: -1}" src; for src in "${@:1:$#-1}"; do case "$src" in -*) continue;; esac; if [ -L "$src" ] && [ "$(readlink -f "$src")" = "$(readlink -f "$dst")/$(basename "$src")" ]; then continue; fi; command mv -f "$src" "$dst" || return; done; }; export -f mv'
  ```
- **Resource files.** This branch names 15 files that no published resource bundle has yet (newest:
  `hmf_pipeline_resources.38_v3.0.0--8`), and checks that each named file exists. For a DNA-only
  Illumina run, set in `params.hmf_data_paths['38']`: `unmap_regions_dna` to
  `common/unmap_regions.38.tsv`, and `[]` for `saga_germline_variants` (and its `_fai`, `_dict`,
  `_img`), `msi_model_error_rates_illumina` (REDUX then skips its MSI prediction; PURPLE still reports
  MSI), the SBX and Ultima MSI files, `sage_pon_ultima`, and the RNA files (`isofox_alt_sj_distribution`,
  `isofox_fusion_cohort`, `isofox_gene_distribution`, `unmap_regions_rna`, `contigs_mapping_rna`).
- **Fixed in the ports** (hmftools `rust-ports`): the REDUX index pointed some reads at the end of a
  compressed block, which the Java tools' htsjdk refuses (BamTools failed); SAGE looked for REDUX's
  BQR file at `/` when the BAM was given as a bare file name, as oncoanalyser does; ESVEE read the
  mapped genome image page by page through the bucket mount (28 minutes at 1% CPU, now 55 seconds).
- **Fixed here:** a missing space in the bwa-mem3 modules' `versions.yml`; the image parameter now
  also counts when set in a `-c` file; `purple:4.5-beta.4` was never published (merged upstream's move
  to beta.5).

Peak memory on the test, for scale only (2M read pairs is about 0.1x): bwa-mem3 17.7 GB (of the
96 GB requested), REDUX 1.3 GB, SAGE 0.4 GB, ESVEE 1.9 GB; of the Java steps, PAVE germline 9.6 GB,
LILAC 9.6 GB and AMBER 6.9 GB used the most.

### Whole genome

A 35x Illumina tumor/normal pair from FASTQ (7 lanes, 3 and 4 per sample; DNA only, `-profile
rust_ports`) was run end to end three times, each time from the FASTQs: on 2026-11-24; on 2026-11-25
with the gcsfuse file cache (below) and ESVEE's correction of finding A5; and on 2026-11-26 with
bwa-mem3's own mode (this fork's default, `--bwamem3_compat none`), the FASTQs left out of the cache,
and the corrections found by the review of 2026-11-25 (hmftools ed4c31e8c6: REDUX R41-R47, SAGE
S49-S51, ESVEE A37-A41 and C18). Every step completed each time. Alignment is bwa-mem3 0.15.0, one task
per lane on n2-standard-64 machines, five lanes at a time; the first two runs used `--compat=bwa-mem2`.
Time is inside the task, peak memory is the task's; costs are each Batch job's machine and boot disk
from the creation of its VM to the end of the job, at us-central1 on-demand N2 prices (storage and
network not included):

| Step | Machine | First run | Second run (cache) | Third run (own mode, FASTQs not cached) | Cost, third run |
|---|---|---|---|---|---|
| Alignment, bwa-mem3 0.15.0, per lane | 64 vCPU | 19-39 min, 62-66 GB | 21-41 min, 62-66 GB | 19-39 min, 62-66 GB | |
| Alignment, all 7 lanes | | 53 min, 3.6 VM hours, $11.85 | 64 min, 4.0 VM hours, $13.03 | 52 min, 3.6 VM hours | $11.58 |
| REDUX, per sample (its lane BAMs as `-input_bam a,b,...`) | 16 vCPU | 39.5 min, 2.4 GB | 38-39 min, 2.3 GB | 33-36 min, 2.4 GB | $1.41 (both) |
| SAGE somatic | 16 vCPU | 16.5 min, 1.9 GB | 14.5 min, 1.8 GB | 16.6 min, 2.0 GB | $0.18 |
| SAGE germline | 16 vCPU | 1.8 min, 0.8 GB | 1.4 min, 0.8 GB | 1.5 min, 0.8 GB | $0.02 |
| ESVEE (one call: prep, assembly, depth, caller) | 16 vCPU | 63.5 min, 12.9 GB | 21 min, 14.3 GB | 16.5 min, 13.2 GB | $0.25 |
| Every other step (Java tools) | | | | | $1.81 |
| Whole run | | 3 h 18 min, $17.21 | 2 h 42 min, $16.78 | 2 h 27 min | $15.25 |

Alignment is most of the cost (69%, 78% and 76%). After alignment the runs cost $5.36, $3.75 and
$3.67.

- **The bucket mount was the bottleneck.** Task inputs are links into the bucket mounted with
  gcsfuse, and without a cache every small read of a file is a network read. ESVEE's assembly reads
  the reference FASTA in many small pieces: 47 of its 63 minutes. From the same prep output, assembly
  took 2.1 minutes with the inputs copied to local disk first, and 2.5 minutes read through the mount
  with the gcsfuse file cache on. A file is cached only if it fits, so a 12 GB cache holds the FASTA
  (3.2 GB), the BWA image (5.6 GB) and the resource files, while the BAMs (18 GB and up) are still
  streamed:

  ```groovy
  google.batch.gcsfuseOptions = ['-o rw', '-implicit-dirs', '--cache-dir=/tmp', '--file-cache-max-size-mb=12288',
                                 '--file-cache-cache-file-for-range-read=true', '--file-cache-enable-parallel-downloads=true']
  ```

  Over the whole second run: ESVEE 63.5 to 21 minutes, CIDER (Java) 54 minutes to 31 seconds, AMBER
  19 to 16, SAGE somatic 16.5 to 14.5; REDUX, which reads only the BAMs, did not change. Alignment was
  slower: by bwa-mem3's own profile it spent 1.5 to 2 times as long reading the FASTQs of the normal
  lanes, so all lanes took 64 minutes instead of 53. The cause is the cache: each FASTQ of the normal
  (12.2 to 12.6 GB) is just under the cache size, so gcsfuse copies the whole file to the boot disk as
  it is read, in parallel downloads, on the disk that the sort writes to; a lane's two FASTQs do not fit
  together. The tumor's FASTQs (16 to 17 GB) are too large to be cached, and those lanes were not slower
  on average. Leaving FASTQs out of the cache fixes it (Alignment, below):

  ```groovy
  google.batch.gcsfuseOptions = ['-o rw', '-implicit-dirs', '--cache-dir=/tmp', '--file-cache-max-size-mb=12288',
                                 '--file-cache-cache-file-for-range-read=true', '--file-cache-enable-parallel-downloads=true',
                                 '--file-cache-exclude-regex=fastq']
  ```

  The pattern is unanchored, so it leaves out every path that contains `fastq`. Google Batch refuses
  mount options with `*`, `(`, `|` or `$` in them ("potential shell injection"), so a narrower regex
  cannot be given.
- **Same inputs, same outputs.** The seven lane BAMs of the two runs are byte-identical (bwa-mem3 is
  deterministic), and so is the normal's REDUX BAM. The tumor's REDUX BAM differed in 7,751 of 1.59
  billion records (60 duplicate flags, the rest bases and qualities of consensus reads), because REDUX
  was given the tumor's three lane BAMs in another order: Nextflow listed them in the order the lanes
  finished, and REDUX's result depends on the order of its inputs, Java's too (finding R41 in
  `redux-rs/REPORT.md`: a tie between two bases of a consensus is given to the first read). The lane
  BAMs are now listed by file name (`read_alignment_dna`, `read_alignment_rna`), and from hmftools
  commit e23fce9d03 redux-rs corrects R41 by default: it reads its inputs in the order of their headers,
  whatever the order given. With the same three lane BAMs given in the two orders of these runs, its
  BAM, index and side files were byte-identical. Downstream, the two runs have the same purity, ploidy, copy number and driver
  catalog, the same 3,824 PASS small variants, and the same 298 PASS somatic structural variant
  breakends.
- **Third run.** With the FASTQs out of the cache, the normal lanes took 19-21 minutes again (bwa-mem3
  read their FASTQs in 401-433 s, against 614-797 s in the second run) and all lanes 52 minutes. Purity,
  ploidy, the somatic and germline driver catalogs and the reported fusions (none) are the same as in the
  second run; copy numbers differ in the fourth decimal and PURPLE has 1,189 segments instead of 1,200.
  PASS small variants 3,824 to 3,823, PASS somatic structural variant breakends 298 to 289, germline
  27,428 to 27,202. To separate the new alignments from the new corrections, the second run's SAGE and
  ESVEE tasks were run again with the third run's image ($0.52): on the second run's BAMs the new image
  gives the third run's PASS calls exactly (SAGE somatic, ESVEE somatic), so every PASS difference comes
  from the corrections and none from bwa-mem3's own mode or the REDUX corrections (whose BAMs mark 16
  fewer duplicates in each sample of 1.2 to 1.3 billion reads, and whose base quality recalibration moves by
  0.03 or less). By correction:
  - SAGE S49 (deletions of 6 or more units of a microsatellite took the MSI quality path): 17 SAGE PASS
    calls become `minTumorQual` and 2 become PASS, all homopolymer deletions of 6 to 23 bases; one of the
    17 had reached PURPLE's PASS calls.
  - ESVEE A40 (fuzzy local match of an assembly's extension): one ESVEE run per new correction, with that
    correction alone taken back to Java's behaviour ($1.58 for six), puts all 9 lost somatic breakends on
    A40: with Java's match they come back, and nothing else does. They are 5 events, all in satellite or
    simple-repeat sequence near centromeres, with quality 32 to 50; Java leaves each assembly unlinked and
    calls a breakend to another satellite copy (or a single breakend); the corrected match links it to the
    local reference within 500 bases, as a duplication or deletion of 51 to 407 bases, which the somatic
    filters then drop. Seven of the 9 are also in the Java run. A40 also accounts for 191 of the 226 fewer
    germline breakends (A37 18, A41 16; P13 and A38 nothing).
- **ESVEE output depended on thread timing** on this pair, whose tumor and normal come from different
  sequencing runs: assembly shortened read names as Java ESVEE does (finding A5 in the port's
  REPORT.md), and two 16-thread runs from the same prep output differed. Corrected by default from
  hmftools commit 1b2682730e (the image of the second run); on this pair the assembly output is then
  identical at 16, 16 and 1 threads. Runs with an earlier image are affected only when the read names
  of tumor and normal (or of the lanes of one sample) come from runs with prefixes of different
  lengths.

### A second whole-genome pair: GIAB HG008

GIAB HG008 is a public tumor/normal pair: a pancreatic ductal adenocarcinoma cell line, HG008-T, and
normal duodenal tissue from the same donor, HG008-N-D. Its NYGC Illumina FASTQs were run on 2026-11-26
with the same profile and settings as the third run above. The image was hmftools aef317426f, which has
the ports on the shared crate `hmf-common-rs` and gives the same results as ed4c31e8c6.

The tumor has 8 lanes from two flowcells, the normal 4. The lanes are named `<flowcell>_<lane>` in the
samplesheet, because the tumor's two flowcells reuse lane numbers. BamTools measured a mean coverage of
71x for the tumor (median 65, 16% duplicates) and 29.5x for the normal (median 30, 11% duplicates). That
is about twice T0's tumor depth.

Alignment is bwa-mem3 0.15.0 in its own mode, one task per lane on n2-standard-64:

| Step | Time in the task, peak memory | Cost |
|---|---|---|
| Alignment, per lane | 11-17 min, 61 GB | |
| Alignment, all 12 lanes (five 64-vCPU machines at a time) | 53 min, 3.4 VM hours, 98 CPU hours | $11.09 |
| REDUX, tumor (8 lane BAMs, 1.75 billion reads) | 58 min, 3.3 GB | |
| REDUX, normal (4 lane BAMs, 0.69 billion reads) | 26 min, 2.2 GB | $1.62 (both) |
| SAGE somatic | 11.9 min, 2.2 GB | $0.13 |
| SAGE germline | 1.4 min, 0.8 GB | $0.02 |
| ESVEE | 7.1 min, 9.5 GB | $0.11 |
| Every other step (Java tools; BamTools, AMBER, COBALT and TEAL prep 10-13 min each) | | $1.31 |
| Whole run | 2 h 33 min | $14.28 |

Alignment was 78% of the cost. Everything after it cost $3.19.

The tumor's REDUX was the longest step after alignment. It read 240 GB and wrote 204 GB at 4 of its 16
CPUs, which is in proportion to T0's tumor (193 and 158 GB in 36 minutes) and bound by I/O.

Every step completed. PURPLE QC passed:
- purity 1.00, as expected for a cell line, and ploidy 1.46;
- female by both AMBER and COBALT, and no contamination;
- 423 copy-number segments, MSS, 3.8 mutations per Mb;
- 10,839 PASS somatic small variants and 381 PASS somatic structural variant breakends.

The somatic driver catalog has KRAS, TP53, SMAD4, NF1, SUZ12 and ARID1B mutations, as expected for a
pancreatic adenocarcinoma.

### AMBER, COBALT and BamTools in the pipeline

The HG008 run was resumed twice on 2026-11-27 with AMBER, COBALT and BamTools as Rust ports (image
hmftools 01257d2b90), into new output folders. Alignment, REDUX, SAGE and ESVEE came from the first
run's cache, so these three steps and everything after them ran again. The first rerun used the default
corrected mode, and the second used `--rust_java_compat true`. Both ran on the same machines as the
first run (n2-standard-16):

| Step | Java (first run): time, peak memory | Rust, corrected | Rust, Java mode |
|---|---|---|---|
| BamTools, tumor | 11.5 min, 3.8 GB | 2.7 min, 0.55 GB | 2.9 min, 0.55 GB |
| BamTools, normal | 4.6 min, 2.8 GB | 1.2 min, 0.53 GB | 1.2 min, 0.54 GB |
| COBALT | 10.3 min, 5.7 GB | 3.5 min, 0.85 GB | 4.2 min, 0.86 GB |
| AMBER | 12.2 min, 10.2 GB | 10.7 min, 1.6 GB | 10.1 min, 1.6 GB |
| The rerun, every step | | $0.49 | $0.50 |

**In Java mode**, the 234 output files are identical to the first run's, apart from:
- the `amber.version` and `cobalt.version` files: version `4.4` and `3.1`, against the beta images'
  `4.4-beta.1` and `3.1-beta.1`, and their build dates;
- BamTools' `partition_stats.tsv`: its timing column and row order;
- the dates and random element ids in the LINX, ORANGE and qsee reports;
- the kataegis cluster ids (`KT`) of 104 records in PURPLE's somatic VCF, a one-to-one renaming of 24
  ids. All of PURPLE's inputs were identical, so PURPLE numbers these clusters nondeterministically
  (85 and 75 records changed in the other reruns).

**In corrected mode**, PURPLE's purity, ploidy, QC, copy-number file and driver catalog are the same as
Java's. 140 lines of its segment file change in the observed-ratio columns (COBALT corrections), and
qsee's mapped proportion drops from 0.9967 to 0.98, still a PASS. Java BamTools inverts the QC flag
of the unplaced unmapped reads, so it counts them as QC-failed and leaves them out of the share
(finding BT3).

AMBER used 13-17% of its 16 CPUs. Java AMBER, as the first version of the port, reads each of about
23,000 regions with its own index query. Through the bucket mount each query becomes a new ranged
request. From hmftools 5edcca6f64 the port reads runs of regions with one query each, with every output
unchanged (amber-rs REPORT.md, section 4). A third rerun with that image, in Java mode, took 2.8 minutes
at 1035% CPU and 1.6 GB, against 10.1 minutes before and Java's 12.2. Its outputs, and everything after
AMBER, were again identical to the Java run's, apart from PURPLE's kataegis ids (75 records this time).
That rerun cost $0.25.

### TEAL, LILAC and CIDER in the pipeline

A fourth rerun of HG008 on 2026-11-27 (image hmftools 6a81b6d5f2, `--rust_java_compat true`) ran TEAL,
LILAC and CIDER as ports, with everything before them from the cache. It cost $0.35, on the first run's
machines:

| Step | Java (first run): time, peak memory | Rust, Java mode |
|---|---|---|
| TEAL prep | 13.8 min at 1112% CPU, 8.6 GB | 15.5 min at 165% CPU, 1.1 GB |
| TEAL pipeline | 35 s, 2.3 GB | 7.4 s, 0.4 GB |
| LILAC | 72 s, 13.4 GB | 9 s, 0.5 GB |
| CIDER | 37 s, 6.2 GB | 20 s, 0.52 GB |

Against the first run's Java outputs:
- **TEAL:** tellength and breakend files identical. The telbams have the same records in the same
  order, and their headers differed only in `@HD VN` (below).
- **LILAC:** the same rows in all three files, in another gene group order (finding LI1: Java orders
  the groups by identity hash).
- **CIDER:**
  - the two extra columns of the hmftools 1.2 sources (Versions);
  - one more read in cider.bam (31,408 against 31,407) and one more in the IGH counts. The hmftools
    1.2 jar includes that read too, so this is the image's other build.
  - the `@HD VN` line.
- **Downstream:** ORANGE's JSON differs only in the order of LILAC's alleles.

The `@HD VN` line: Java REDUX, TEAL and CIDER write `VN:1.6` whatever the input BAM says, because
htsjdk's header decoder does not keep the input's version. The three ports kept the input's, which is
`VN:1.3` from bwa-mem3 (so the REDUX BAMs of the Rust runs above say 1.3, the Java run's 1.6) and 1.5
in the HG008 BAMs. Fixed in hmftools 5e11a56926; on a `VN:1.5` slice of the HG008 tumor BAM the three
ports' BAMs are then identical to Java's.

TEAL's prep was I/O bound on Batch at 165% CPU. It read the tumor and then the normal, each as one
stream, at about 160 MB/s through gcsfuse, while Java's 15 reader threads read in parallel.

From hmftools 3bb279d27e, several file handles read each BAM ahead in 32 MB chunks, which the parser
still takes in file order. A fifth rerun ran TEAL and CIDER from that image, for $0.18:
- **TEAL prep:** 7.5 minutes at 420% CPU and 1.3 GB, against Java's 13.8 minutes.
- **TEAL outputs:** the telbams match Java's header and records as SAM text. 11,792 unplaced records still
  differed in their BAM bin: htsjdk writes 0 for a record without a reference, and the inputs had
  htslib's 4680. teal-rs and redux-rs write 0 from hmftools 44720f20ef, and on a public test BAM with
  unplaced reads both are then byte-identical to Java.
- **CIDER:** cider.bam's header now matches Java's; the one extra read and the two columns remain (the
  image's other 1.2 build).

### LINX in the pipeline

A seventh rerun (2026-11-27, hmftools 044933dff1) ran LINX_SOMATIC and LINX_GERMLINE as linx-rs and
TEAL from the same image, in Java mode, for $0.25. Everything before them came from the cache, PAVE
and PURPLE as Java.
- **LINX:** every file is identical to linx-rs run locally on the same inputs, which equals Java with
  `-XX:hashCode=2`.
  - Against the first run's Java LINX: 13 of the 20 files are identical. The rest differ only in
    cluster, chain and breakend ids and in germline row order, which Java takes from identity hashes
    (LX1); sorted and without those ids, they are identical.
  - 2.8 and 3.1 s at 164 and 89 MB, against Java's 6 and 7.5 s at 1 GB.
- **Downstream:**
  - ORANGE's JSON is unchanged from the previous rerun; it differs from the Java run only in LILAC's
    allele order.
  - The LINX visualiser's chr8 plot, the LINXREPORT HTML and ORANGE's PDF show the cluster ids, so
    they differ.
- **TEAL:** both telbams are now byte-identical to the Java run's once decompressed. The prep took 6.9
  minutes at 1.5 GB, against Java's 13.8 minutes at 8.6 GB.

Reruns 3 to 7 had also each copied the REDUX BAMs (151 GB) into their output folders, though they took
REDUX from the cache: the pipeline publishes through its `output` block, so a process `publishDir`
setting does not stop it. Since then, `--publish_alignments false` leaves them out.

### PAVE and PURPLE in the pipeline

An eighth rerun (2026-11-27, hmftools af4949c5f9, `-profile rust_ports,rust_purple`) ran PAVE and PURPLE as
ports in Java mode, with LINX, for $0.16. Everything before them came from the cache.
- **PAVE:** all 17,254 somatic and 846 germline records are identical to the Java run's. Only the
  `##PaveVersion` line differs (2.0 against 2.0-beta.2).
  - Somatic: 6 s and 211 MB, against Java's 31 s and 6.3 GB.
  - Germline: 4.8 s and 189 MB, against 37 s and 9.8 GB.
- **PURPLE:**
  - Every text file is identical.
  - The VCF records are identical, apart from the kataegis ids of 101 somatic records (24 clusters
    renamed one to one; Java numbers them differently from run to run). The headers differ in the
    version lines.
  - 11.9 s and 1.0 GB, against Java's 44 s and 5.3 GB with its plots.
- **Plots and ORANGE:** no plots or Circos files are written. oncoanalyser runs ORANGE only when PURPLE's
  circos plot exists, so it skipped ORANGE, and the run still completed.
- **Later steps:** LINX, CHORD, SIGS, LILAC, NEO, CUPPA, PEACH and Virus Interpreter gave the same
  outputs as with Java PURPLE.

This rerun was the first with `--publish_alignments false`. Its first pass published nothing: a `>>` inside
an `if` in the `output` block is not a publish statement. The filter now sits on the published channel
(fork d248e3d4), and a second pass published the results without the BAMs.

### SIGS, Virus Interpreter, PEACH and CHORD in the pipeline

A ninth rerun (2026-11-27, hmftools a9ca2fafda) ran SIGS, Virus Interpreter, PEACH and CHORD as ports, in
Java mode, for $0.02 in tasks; everything else came from the cache. Every output of the four is identical to
the Java run's, and CUPPA and ORANGE, which ran again on them, gave the same outputs as in the previous
reruns (ORANGE's JSON differs from the Java run only in LILAC's allele order).

| Step | Port | Java |
|---|---|---|
| SIGS | 0.6 s, 2.8 MB | 4 s, 266 MB |
| VIRUSINTERPRETER | 1.6 s, 9.2 MB | 5 s, 598 MB |
| PEACH | 0.6 s, 2.8 MB | 3 s, 80 MB |
| CHORD | 7.4 s, 51 MB (mostly staging) | 9 s, 439 MB (with R) |

CHORD's module still asks R and Rscript for their versions in `versions.yml`; in the image those lines stay
empty and the task succeeds.

### NEO in the pipeline

A tenth rerun (2026-11-27, hmftools 1111f6bc2e) ran NEO_FINDER and NEO_SCORER as ports, in Java mode, for $0.01
in tasks. All three NEO files are identical to the Java run's. The finder took 4.0 s and 105 MB (Java: 6 s,
1.0 GB), the scorer 3.1 s and 134 MB (Java: 16 s, 1.3 GB).

### CUPPA and qsee in the pipeline

An eleventh rerun (2026-11-27, hmftools 3f31dc583a) ran CUPPA and qsee as ports, in Java mode, for $0.02 in tasks.
- **qsee:** the three data files are identical to the Java run's. Without R in the image the report is not drawn.
  2.6 s and 21 MB, against 38 s and 1.1 GB.
- **CUPPA:** the features file is identical. The predictions differ from the pipeline's pycuppa (Linux x86-64)
  only in the last digits (at most 2.4e-13 relative), as pycuppa's own do between machines. No plot is drawn.
  1.7 s and 21 MB, against 16 s and 494 MB.
- **ORANGE** ran on them: its JSON carries those last digits of CUPPA's and has no CUPPA or qsee plot; the rest of
  it is unchanged.

Two image builds before this one failed their checks:
- with rustc 1.100.0 a release build of cuppa-rs mistyped its test VCF's SNVs, so the image is now pinned to
  rustc 1.97.1, the version every port was checked with;
- the qsee wrapper added `-skip_plot` to `qsee -version` too.

### ORANGE in the pipeline

A twelfth rerun (2026-11-27, hmftools 01f01a5f35) ran ORANGE as the port, in Java mode, on rerun 11's inputs (PURPLE
in Java, for its plots); only ORANGE and MultiQC ran. The JSON is identical to rerun 11's from Hartwig's ORANGE image
(5.1-beta.1) but for one key: the image names the virus section `viridian`, the port `virusInterpreter`, as the 5.1
jar of hmftools ca04e6c39d does. Upstream renamed it on 2026-11-23 (f5f660cc9f, which also lets ORANGE read Viridian's
files); the port is to follow. The LINX plot order matched too (both read Linux's directory order in Java mode). No
PDF is written. 5.3 s and 42 MB, against 18 s and 419 MB.

### PURPLE's plots and ORANGE in the pipeline

A fourteenth rerun (2026-11-27, hmftools 617f926b53, `-profile rust_ports,rust_purple`) ran PURPLE as the port
with its 9 plots drawn natively (hmf-plot-rs on BUI's software raster: the 7 R plots and the 2 Circos images, no R,
Perl or Circos in the image), and ORANGE as the port on them, in Java mode; every step after PURPLE ran again on
the images of reruns 7-11. $0.14 in tasks.
- **PURPLE:** the text files are identical to Java's; the VCFs differ only in `##purpleVersion` (4.5 against the
  image's 4.5-beta.5) and in KT on 94 of 17,254 somatic records (PU1); baf.circos is in genome order (PU2). The
  plots differ from R's and Circos's in rendering only (see `purple-rs` REPORT.md). 14.7 s and 868 MB, against
  Java's 43.9 s and 3.7 GB with R and Circos.
- **ORANGE** ran on the native plots: its JSON is byte-identical to Java ORANGE's (rerun 11, the same inputs);
  no PDF is written.
- **Every later step** (LINX, LINX's visualiser, LILAC, TEAL, NEO, CUPPA, CHORD, SIGS, PEACH, Virus Interpreter,
  qsee, the SAGE visualiser) wrote the same files as with Java PURPLE; only gpgr's LINX HTML report differs.

With that, PURPLE joined `rust_ports`. A thirteenth rerun had checked ORANGE after it followed upstream's
f5f660cc9f (the virus section `viridian`): its JSON equalled Java's from the same inputs, byte for byte.

### Alignment: bwa-mem3 mode and the file cache

One normal lane of the pair (25 GB of FASTQ, about 160 million read pairs) was aligned three more times on
2026-11-25, with the rerun's command on the same machine type (n2-standard-64), the three at the same
time so that they shared the bucket's conditions:

| bwa-mem3 0.15.0 | gcsfuse file cache | Task time | FASTQ reading (bwa-mem3's profile) | Output |
|---|---|---|---|---|
| `--compat=bwa-mem2` | off | 20.2 min | 362 s | the rerun's lane BAM, byte for byte |
| `--compat=bwa-mem2` | on (FASTQs cached) | 21.2 min | 577 s | the same |
| own mode | off | 19.8 min | 358 s | every record the same as compat mode's apart from the added `MQ` and `HN` tags; reads that start at one position in another order |

The same lane took 19.9 minutes in the first run (no cache) and 27.8 in the second (cache, five lanes
at a time, 797 s reading). A test job showed that the cache copies a FASTQ in full as soon as it is
read (3.9 GB within 30 seconds of a 200 MB read), and that with `--file-cache-exclude-regex=fastq` it
copies none while still caching the reference files. The test jobs cost $3.87.

The same pair was also run with the Java tools (`-profile` without `rust_ports`), on the same bwa-mem3
alignments taken from the first run's cache, without the gcsfuse cache; against the first Rust run:

| Step | Java: time, peak memory | Rust: time, peak memory |
|---|---|---|
| Alignment, bwa-mem3 0.15.0, 7 lanes (the same tasks) | 53 min, $11.85 | 53 min, $11.85 |
| REDUX, per sample | 1.2-1.4 h, 30 GB | 39.5 min, 2.4 GB |
| SAGE somatic | 3.7 h, 61 GB | 16.5 min, 1.9 GB (14.5 min with the gcsfuse cache) |
| ESVEE | 1.2 h, 108 GB (above its 96 GB request; the machine had 128 GB) | 63.5 min, 12.9 GB (21 min with the cache) |
| Every step after alignment, cost | $11.77 | $5.36 ($3.75 with the cache) |
| Whole run, cost | $23.62 | $17.21 ($16.78 with the cache) |

Results agree: purity, ploidy and contamination are the same; REDUX marks the same duplicates to
within 16 reads of 2.5 billion and BamTools' coverage metrics are identical; of the PASS somatic
small variants 3,755 are shared (3,813 Java, 3,824 Rust); ESVEE has 283 and 298 PASS somatic breakends.
The two runs are not expected to be identical: the Rust run corrects the faults found in the Java
code (the default), and this branch's Java images are pre-release versions older than the code the
ports follow (see Versions). Each was checked against Java of its own version on this pair:

- **ESVEE.** The Java ESVEE task's command rerun in the Rust image with `ESVEE_JAVA_COMPAT=1`, on the
  same inputs: 283 PASS somatic breakends against Java's 283 (281 at the same position, all within 10
  bases) and 27,504 germline against 27,505; the rest differ in ids and in a few support counts, as two
  multi-threaded Java runs do. So the 298 of the default run come from the corrections.
- **SAGE.** Java SAGE built from the commit the port follows, and the port with `SAGE_JAVA_COMPAT=1`,
  on three chromosomes with the task's options: the same 45,774 records, differing only in the local
  phase set ids and their read counts, which Java numbers in thread order. Java took 23 minutes, the
  port 1.6. The SAGE 5.1-beta.3 image differs from that commit by about 30,000 records genome-wide;
  among them, it writes no record for an in-frame deletion in a driver gene that SAGE 5.1, SAGE 4.2 and
  the port all call PASS.

## Checking the wiring

A stub run (no real tools; checks that every step is connected) from FASTQ input:

```bash
mkdir -p temp
for s in tumor.dna normal.dna tumor.rna; do for r in R1 R2; do echo | gzip > temp/subject_a.$s.$r.fastq.gz; done; done
nextflow run <path to this repo> -profile test_stub,rust_ports -stub \
    --input <path to this repo>/tests/rust_ports/fastq.stub.csv --outdir out
```

It runs bwa-mem3 for the two DNA samples, STAR, samtools sort and MarkDuplicates for the RNA sample (with
`--rna_aligner bwa-mem3`: bwa-mem3, TARS and REDUX), no fastp, and every later step.

To check that `--rust_ports_bin` takes these steps out of their containers, run REDUX, SAGE and ESVEE
alone from BAM input under the `docker` profile; with the Docker daemon stopped they still complete,
and without `--rust_ports_bin` they fail to reach Docker:

```bash
mkdir -p temp
for s in tumor normal; do : > temp/subject_a.$s.bam; : > temp/subject_a.$s.bam.bai; done
nextflow run <path to this repo> -profile test_stub,docker,rust_ports -stub \
    --input <path to this repo>/tests/rust_ports/bam.stub.csv --outdir out \
    --processes_manual redux,sage,esvee --rust_ports_bin ~/rust-ports
```

A stub run checks the wiring only; it does not run the tools' commands.
