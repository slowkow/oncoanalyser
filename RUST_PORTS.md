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
| `--rust_java_compat` | `false` | reproduce the Java tools' output exactly, including the faults the ports correct by default |

By default the ports correct the faults found in the Java code (each is listed with an id in the
port's `REPORT.md`). `--rust_java_compat true` sets `SAGE_JAVA_COMPAT`, `REDUX_JAVA_COMPAT` and
`ESVEE_JAVA_COMPAT` in the task environment, which makes each port behave as Java does.

The profiles also lower the memory requested for these steps (SAGE and REDUX 8 GB, ESVEE 24 GB, each
doubled on retry), since the ports need a few GB where the Java tools are given 36 to 72 GB. More
samples then fit on one machine at once.

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
