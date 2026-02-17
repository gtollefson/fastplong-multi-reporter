# fastplong_multireport

Aggregate fastplong QC reports from many samples into a single interactive HTML report for Oxford Nanopore long-read data.

## Overview

[Fastplong](https://github.com/OpenGene/fastplong) produces per-sample JSON QC reports. **fastplong_multireport** scans a directory of these reports, aggregates metrics across samples, and generates a self-contained HTML report with interactive Plotly charts to help identify batch effects, quality trends, and outliers across a run.

## Features

- **Directory-based**: Point at a fastplong results directory; discovers all `*_fastplong_report.json` files (optionally recursive)
- **Self-contained HTML**: Plotly.js is embedded inline—report works offline after transfer from remote servers (scp, rsync, download)
- **Professional layout**: Sticky sidebar navigation, scrollable sections, clean header
- **Interactive charts**: Hover to identify samples, detect batch effects, spot outliers

### Charts included

| Section | Description |
|---------|-------------|
| Summary table | Sample × metrics (reads, retention %, mean length, Q20/Q30, GC %, filtering stats) |
| Quality by position | Mean base quality per position (all samples overlaid; hover to identify) |
| Read counts | Before vs after filtering (grouped bar) |
| Retention rate | Per-sample retention % |
| Mean read length | After filtering (bar chart) |
| Q20/Q30 rates | Base quality rates per sample |
| Length vs quality | Scatter plot (batch effect / outlier detection) |
| Filtering breakdown | Stacked bar (passed, low quality, too short, too long) |

## Requirements

- Python 3.9+
- [pandas](https://pandas.pydata.org/) ≥ 1.5
- [plotly](https://plotly.com/python/) ≥ 5.0

## Installation

### From source (clone or download)

```bash
cd fastplong_multireport
pip install -r requirements.txt
```

### With conda

```bash
conda install -c conda-forge pandas plotly
```

## Usage

### Command line

```bash
# Basic: scan a directory, write report to same dir (default: fastplong_multireport.html)
python -m fastplong_multireport /path/to/fastplong/results/

# Specify output path
python -m fastplong_multireport /path/to/fastplong/results/ -o report.html

# Only search top-level (no subdirs)
python -m fastplong_multireport /path/to/fastplong/results/ --no-recursive

# Custom report title
python -m fastplong_multireport /path/to/fastplong/results/ -t "My Run QC Report"
```

### Python module

```python
from pathlib import Path
from fastplong_multireport.fastplong_multireport import (
    discover_json_reports,
    load_fastplong_reports,
    generate_report,
)

results_dir = Path("results/qc/raw")
json_files = discover_json_reports(results_dir, recursive=True)
data = load_fastplong_reports(json_files)
generate_report(data, Path("report.html"), title="QC Report")
```

## Input

- **Directory** containing `*_fastplong_report.json` files (from fastplong preprocessing)
- Sample names are inferred from filenames (e.g. `sample1_fastplong_report.json` → `sample1`)

## Output

- Single HTML file with all charts embedded
- Fully self-contained (no CDN); works offline after transfer
- Typical size: ~5–10 MB (including embedded Plotly.js)

## Example: Kelch13 haplotype pipeline

If your pipeline writes fastplong JSONs to `results_260205/qc/raw/`:

```bash
python -m fastplong_multireport results_260205/qc/raw/ \
  -o results_260205/qc/fastplong_aggregate_report.html
```

## License

MIT

## Author

George A. Tollefson, 2026. Developed for the Kelch13 haplotype analysis pipeline (Oxford Nanopore K13 flanking amplicon sequencing).
