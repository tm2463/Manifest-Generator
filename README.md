# Manifest Generator

Wellcome to Manifest Generator, the swiss army knife of manifest generation! Manifest generator is a command-line tool for generating sample manifests from FASTQ files, supporting short, long, and hybrid sequencing read types.

---

## Requirements

- Python 3.10+
- `pandas`
- [`lrge`](https://github.com/mbhall88/lrge) _(optional, required for genome size estimation)_
- [`shelf`](https://gitlab.internal.sanger.ac.uk/pam/shelf) _(optional, required for `--from_lanes`)_

Shelf must be loaded separately if using '--from_lanes' by running `module load shelf`

---

## Usage

```bash
manifest_generator.py [--from_dir | --from_dir_recursive | --from_paths | --from_lanes | --from_paths_id] [flags]
```

At least one input source is required. Multiple input sources can be chained together. Hybrid read manifests can only be built from --from_paths_id

---

## Input Sources

| Source                                   | Description                                                                             |
| ---------------------------------------- | --------------------------------------------------------------------------------------- |
| `--from_dir <dir> [<dir> ...]`           | Collect FASTQs by searching one or more directories                                     |
| `--from_dir_recursive <dir> [<dir> ...]` | Collect FASTQs by recursively searching one or more directories(requires `--max_depth`) |
| `--from_paths <file>`                    | Supply a list of FASTQ file paths from which to build a manifest                        |
| `--from_lanes <file>`                    | Query FASTQ paths from Shelf using a lanes file                                         |
| `--from_paths_id <file>`                 | Supply a csv containing "ID" and "path" headers                                         |

---

## Flags

| Flag                   | Default        | Description                                             |
| ---------------------- | -------------- | ------------------------------------------------------- |
| `-m`, `--mode`         | `short`        | Manifest type: `short`, `long`, or `hybrid`             |
| `-o`, `--outdir`       | `./`           | Output directory                                        |
| `-n`, `--name`         | `manifest.csv` | Now YOU can name your very own manifest!                |
| `-d`, `--max_depth`    | `None`         | Max directory depth for `--from_dir_recursive`          |
| `--genome_size`        | `False`        | Estimate genome size for long reads using `lrge`        |
| `-t`, `--threads`      | `1`            | Threads for genome size estimation                      |
| `--duplicate_handling` | `error`        | How to handle duplicate IDs: `error`, `warn`, or `drop` |
| `-v`, `--verbose`      | `False`        | Enable debug logging                                    |

---

## Modes

### `short`

Pairs R1 and R2 short reads by sample ID.

Output columns: `ID`, `R1`, `R2`

### `long`

Assigns long reads to sample IDs.

Output columns: `ID`, `long_fastq`, `genome_size`

### `hybrid`

Combines short read pairs with a long read per sample.

Output columns: `ID`, `R1`, `R2`, `long_fastq`, `genome_size`

Hybrid manifests can only be built from --from_paths_id, see Read ID Detection below for details

---

## Read ID Detection

Short read sample IDs are inferred from FASTQ filenames by:

1. Stripping all extensions (e.g. `.fastq.gz`)
2. Removing R1/R2 suffixes matched by the patterns `_R?1` and `_R?2`

Long read sample IDs are inferred from FASTQ filenames by stripping all file extensions only.

Read type (R1, R2, or long) is assigned based on whether the filename matches the R1 or R2 pattern.

Hybrid manifests require the user to provide a sample ID for each read in order to group short and long reads together. Example CSV input:

```csv
ID,read
sample_1,<path/to/sample_1_R1.fq>
sample_1,<path/to/sample_1_R2.fq>
sample_1,<path/to/sample_1_long.fq>
sample_2,<path/to/sample_2_R1.fq>
sample_2,<path/to/sample_2_R2.fq>
sample_2,<path/to/sample_2_long.fq>
```

---

## Duplicate Handling

Duplicates are detected based on expected ID frequency per mode:

| Mode     | Expected occurrences per ID |
| -------- | --------------------------- |
| `long`   | 1                           |
| `short`  | 2 (R1 + R2)                 |
| `hybrid` | 3 (R1 + R2 + long)          |

The `--duplicate_handling` flag controls behaviour when duplicates are found:

| Method  | Behaviour                                     |
| ------- | --------------------------------------------- |
| `error` | Log error and exit (default)                  |
| `warn`  | Log warning and continue with paired reads    |
| `drop`  | Drop duplicates, keep first occurrence per ID |

Reccommended usage:

| Method  | Behaviour                                                                             |
| ------- | ------------------------------------------------------------------------------------- |
| `error` | When you are unsure if duplicates are present                                         |
| `warn`  | When you are unsure if duplicates are present and wish to inspect the output yourself |
| `drop`  | Drop duplicates, keep first occurrence per ID                                         |

---

## Examples

```bash
# Short read manifest from a directory
manifest_generator.py --from_dir /path/to/fastqs -m short -o ./output

# Hybrid manifest, searching recursively up to depth 3
manifest_generator.py --from_dir_recursive /path/to/project -d 3 -m short -o ./output

# Long read manifest with genome size estimation
manifest_generator.py --from_dir /path/to/fastqs -m long --genome_size -t 8

# From a list of paths, warn on duplicates
manifest_generator.py --from_paths reads.txt --duplicate_handling warn

# From Shelf lanes file
manifest_generator.py --from_lanes lanes.txt -m short -o ./output

# Hybrid read manifest
manifest_generator.py --from_paths_id paths_id.csv -m hybrid -o ./output
```

---

## Output

A `.csv` manifest is written to `--outdir` with the filename specified by `--name`. A timestamped log file is also written to the same directory.

```
output/
├── manifest.csv
└── 2025-01-01_12-00-00.log
```
