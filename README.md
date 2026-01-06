# `ObCheck`
quick and crude reads to pHMM mapping tool to probe for potential Oblin-1-coding reads (suggestive of Obelisk presence)

`ObCheck` is a script-afied (thanks `Claude`) version of the pHMM mapping that was done in the original [Obelisk discovery paper](https://www.cell.com/cell/abstract/S0092-8674(24)01091-2)

This tool iterates over a set of directories, each containing paired-end reads with the naming convention: `$DIRNAME_[12].fastq.gz`

Sequentially, the tool naively translates all reads in all 6 reading frames, and concatenates this to a `$DIRNAME.aa` file

The tool then compares these translated reads against the original `Oblin-1` pHMM (specified by `-h PATH`) using `hmmer` and then picks all reads that map at the default `E-val` threashold of `0.01`

This is likely very permissive and might be susceptible to false-positives.

The resulting counts (and fractional counts) are written to a single `pHMM.counts` tab-seperated file.

One can specify threads with `-t`

THIS TOOL DOES NO READ QC, HAS A HARD CODED E-VAL THRESHOLD, USES THE OUT-OF-DATE/ORIGINAL pHMM, AND LIKELY HAS MANY OTHER ERRORS (AND SPELLING MISTOOKS) - USE WITH CAUTION AND WITH NO GUARANTEES OF FUNCTIONALITY

# What is `ObCheck` for?

`ObCheck` is for a crude first-pass to ask: _"might there be Oblin-1 coding reads (likely derived from Obelisks) in my sequencing data?"_

If you see ample counts, it's likely worthwhile to follow up on these samples much more carefully

## if you have `ObCheck`-positive data, I would love to chat! Please reach out to me!

# What is `ObCheck` NOT for?

`ObCheck` is not suitable for any other downstream analysis, this is really for a first pass check only.

For example, if you want to do counts-based statistics including these putative Obelisks, you will need to _de novo_ (cross-)assemble new Obelisk genomes from your data (appropriately cluster, annotate, and rotate them); classify them against the already known Obelisks (so as to not accidentally invoke a new naming convention); append them to your pre-existing analytical pipeline (e.g. `Kraken2`/`MetaPhlAn`); and re-run your counting/classification prior to any kind of statistical analysis.

This is because 1) this kind of pHMM mapping is substantially different from the counting/classification methods conventionally used and so the counts are not comparable on a technical level, and 2) this kind of mapping attempts to map _accross_ the Oblin-1 space (as described in the intial paper, newer pHMMs are now avaliable) and so you lose any kind "genus"/"species" level information, meaning your counts are likely not comparable on a biological level.

Lastly, `ObCheck` does not prove or disprove the presence of Oblin-1 coding reads in your data - false positives are very possible, as are false negatives.

# Installation / Usage

create `mamba` environment

`mamba create -fy -n ObCheck -c bioconda -c conda-forge eza tree seqkit hmmer pigz bc`

`mamba activate ObCheck`

clone repo

`git clone https://github.com/Zheludev/ObCheck.git`

`chmod +x ObCheck.sh`

run

`$REPOPATH/ObCheck/ObCheck.sh -t 32 -h $REPOPATH/ObCheck`

example file structure

```
├── SRR29668074
│   ├── SRR29668074_1.fastq.gz
│   └── SRR29668074_2.fastq.gz
├── SRR29668075
│   ├── SRR29668075_1.fastq.gz
│   └── SRR29668075_2.fastq.gz
├── SRR29668076
│   ├── SRR29668076_1.fastq.gz
│   └── SRR29668076_2.fastq.gz
├── SRR29668077
│   ├── SRR29668077_1.fastq.gz
│   └── SRR29668077_2.fastq.gz
├── SRR29668078
│   ├── SRR29668078_1.fastq.gz
│   └── SRR29668078_2.fastq.gz
└── SRR29668079
    ├── SRR29668079_1.fastq.gz
    └── SRR29668079_2.fastq.gz
```

example output `.tsv`

```
sample	Ob1	reads	Ob1_fxn
SRR29668074	4489	914765	.00490
SRR29668075	3760	732540	.00513
SRR29668076	2726	710719	.00383
SRR29668077	0	877056	0
SRR29668078	0	929365	0
SRR29668079	0	922258	0
```

# version

as of 09/18/25 this is v0.0.1
