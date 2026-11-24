# oncoanalyser with bwa-mem3 and the Rust ports of hmftools

This fork of nf-core/oncoanalyser (branch `rust-ports`) runs the pipeline with:

- **bwa-mem3** for read alignment. The branch starts from upstream's `bwa-mem3` branch (Hartwig's
  development line for oncoanalyser 3.1), which already replaces BWA-MEM2 with bwa-mem3. Upstream runs
  it with `--compat=bwa-mem2`, which reads the BWA-MEM2 index and writes the same records as BWA-MEM2.
- **No fastp.** fastp only split FASTQ files into chunks for parallel alignment (or extracted UMIs);
  it changed no reads. `max_fastq_records` defaults to `-1`, so FASTQs are aligned whole and fastp
  does not run unless UMIs are enabled.
- **Rust ports of SAGE, REDUX and ESVEE** in place of the Java tools, chosen with a profile. The ports
  live in the hmftools fork (https://github.com/slowkow/hmftools, branch `rust-ports`): `sage-rs/`,
  `redux-rs/`, `esvee-rs/`. Each takes the Java tool's options, writes the same files, and on every
  test input gives Java's output exactly when asked to (see below).

## Running

```bash
nextflow run slowkow/oncoanalyser -r rust-ports \
    -profile docker,rust_ports \
    --mode wgts --genome GRCh38_hmf \
    --input samplesheet.csv --outdir output/
```

| Profile | Effect |
|---|---|
| `rust_ports` | SAGE, REDUX and ESVEE all run as Rust ports |
| `rust_sage`, `rust_redux`, `rust_esvee` | one tool as a Rust port, the others in Java; combine as needed |

| Parameter | Default | Meaning |
|---|---|---|
| `--rust_ports_container` | `hmftools-rust:dev` | image with the ports (built locally, see below) |
| `--rust_ports_bin` | not set | folder with the ports on the host; these steps then run without a container (see below) |
| `--rust_java_compat` | `false` | reproduce the Java tools' output exactly, including the faults the ports correct by default |
| `--rust_low_base_qual_threshold` | not set (26) | SAGE and ESVEE: bases at or below this quality are low quality |
| `--rust_high_base_qual_threshold` | not set (30) | SAGE and ESVEE: Illumina bases at or above this quality are high quality |

Give `--rust_ports_bin` and `--rust_java_compat` on the command line or in a `-params-file`: they are
read while the configuration is loaded, so a value set in a `-c` config file is not seen.
`--rust_ports_container` can be set in either place (also in a `-c` file). It is read as each task
starts, and Nextflow then warns once at start-up, before the pipeline's parameters are loaded:
`Access to undefined parameter rust_ports_container`. That warning is harmless.

By default the ports correct the faults found in the Java code (each is listed with an id in the
port's `REPORT.md`). `--rust_java_compat true` sets `SAGE_JAVA_COMPAT`, `REDUX_JAVA_COMPAT` and
`ESVEE_JAVA_COMPAT` in the task environment, which makes each port behave as Java does.

The two base-quality parameters are options of the ports, not of the Java tools, so they apply only to
the steps that run as ports. They are for data whose base qualities were rewritten by recalibration and
stop below Q30: Java SAGE then passes no SNV, since it counts only bases at Q30 or above towards
`minTumorSupport`, and Java ESVEE treats most of the tumor's bases as low quality, so the tumor
supports few SVs. Each port's README says what the bands do. When neither is set the task commands are
unchanged (so is the `-resume` cache); they can be given anywhere, including a `-c` file. Hartwig
tuned its filters on data with the usual qualities, so calls made with moved bands need checking.

The profiles also lower the memory requested for these steps (SAGE and REDUX 8 GB, ESVEE 24 GB, each
doubled on retry), since the ports need a few GB where the Java tools are given 36 to 72 GB. More
samples then fit on one machine at once.

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
(`process_high`); the memory requests are this fork's (upstream asks 72 GB). The one-call ESVEE mode
has not been measured on a whole genome yet: assembly keeps its memory until the run ends, so the peak
is expected near assembly's.

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
a container. With `--rust_ports_bin <folder>`, SAGE, REDUX and ESVEE run on the host from that folder,
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

The image holds `sage-rs`, `redux-rs` and `esvee-rs`, wrapper commands `sage`, `redux` and `esvee`
that accept the Java tools' command lines (JVM options are dropped, Java main classes become
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

## Known gaps

- The modules also record the Java version in `versions.yml` (`java --version`). The image has no
  Java, so that entry is empty; the step is not affected.
- On macOS (`--rust_ports_bin`), the modules read tool versions with GNU `sed` syntax, which the
  system `sed` rejects, so the versions in `versions.yml` are empty unless GNU sed is first on `PATH`
  (for example Homebrew's `gnu-sed`). The steps themselves are not affected.
- redux-rs does not draw the base-quality recalibration plot (`-bqr_write_plot`); no later step
  reads it.
- Run end to end on Google Batch only on a 2M-read-pair test (below); not yet on a whole genome or on
  AWS Batch, so the whole-genome memory of each step inside the pipeline (and of the one-call ESVEE
  mode) is not yet measured.

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

## Checking the wiring

A stub run (no real tools; checks that every step is connected) from FASTQ input:

```bash
mkdir -p temp
for s in tumor.dna normal.dna tumor.rna; do for r in R1 R2; do echo | gzip > temp/subject_a.$s.$r.fastq.gz; done; done
nextflow run <path to this repo> -profile test_stub,rust_ports -stub \
    --input <path to this repo>/tests/rust_ports/fastq.stub.csv --outdir out
```

It runs bwa-mem3 for the three samples, no fastp, and every later step.

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
