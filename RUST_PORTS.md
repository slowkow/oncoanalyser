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

By default the ports correct the faults found in the Java code (each is listed with an id in the
port's `REPORT.md`). `--rust_java_compat true` sets `SAGE_JAVA_COMPAT`, `REDUX_JAVA_COMPAT` and
`ESVEE_JAVA_COMPAT` in the task environment, which makes each port behave as Java does.

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
subcommands), and samtools and sambamba, which REDUX and ESVEE are given as `-bamtool`. Build it from
the hmftools fork:

```bash
git clone -b rust-ports https://github.com/slowkow/hmftools.git
cd hmftools
docker build -f rust-docker/Dockerfile -t hmftools-rust:dev .
```

The image is not published to a registry yet; build it on the machine that runs the pipeline, or
push it to your own registry and pass `--rust_ports_container`.

## Versions

The Java images in this branch (`conf/hmf_containers.config`) are Hartwig's pre-release builds:
SAGE 5.1-beta.3 and ESVEE 2.1-beta.3; REDUX is the bioconda 2.0.7. The ports reproduce SAGE 5.1 and
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
- The image has been built for the build machine's architecture only.

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
