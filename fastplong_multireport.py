#!/usr/bin/env python3
"""
fastplong_multireport – Aggregate fastplong QC reports into a MultiQC-style interactive HTML report.

Scans a fastplong results directory for *_fastplong_report.json files, aggregates metrics,
and produces a single HTML report with interactive Plotly charts for:
- Summary table (sample × metrics)
- Mean base quality by position (all samples, hover to identify)
- Read count and retention bar charts
- Mean read length, Q20/Q30 rates per sample
- Scatter plots for batch effect detection
- Filtering breakdown per sample

Usage:
    fastplong_multireport /path/to/fastplong/results/
    fastplong_multireport /path/to/fastplong/results/ -o report.html
    fastplong_multireport /path/to/fastplong/results/ --recursive
"""

import argparse
import json
import sys
from pathlib import Path
from datetime import datetime

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


def _get_plotlyjs() -> str:
    """Get Plotly.js as string for embedding. Works across Plotly versions."""
    try:
        from plotly.offline import get_plotlyjs
        return get_plotlyjs()
    except (ImportError, AttributeError):
        pass
    try:
        from plotly.io import get_plotlyjs
        return get_plotlyjs()
    except (ImportError, AttributeError):
        pass
    # Fallback: read from package
    import plotly
    js_path = Path(plotly.__file__).parent / "package_data" / "plotly.min.js"
    if js_path.exists():
        return js_path.read_text(encoding="utf-8", errors="replace")
    raise RuntimeError("Could not load Plotly.js. Ensure plotly is installed.")


def discover_json_reports(results_dir: Path, recursive: bool = True) -> list[Path]:
    """Discover all *_fastplong_report.json files in the results directory."""
    pattern = "**/*_fastplong_report.json" if recursive else "*_fastplong_report.json"
    return sorted(results_dir.glob(pattern), key=lambda p: p.name)


def load_fastplong_reports(qc_files: list[Path]) -> dict:
    """Load all fastplong JSON reports into a structured dict. Sample name from filename."""
    data = {}
    for qc_file in qc_files:
        name = qc_file.stem.replace("_fastplong_report", "")
        try:
            with open(qc_file) as f:
                data[name] = json.load(f)
        except Exception as e:
            print(f"Warning: Could not read {qc_file}: {e}", file=sys.stderr)
    return data


def build_summary_df(data: dict) -> pd.DataFrame:
    """Build sample × metrics DataFrame."""
    rows = []
    for sample, d in data.items():
        s = d.get("summary", {})
        bf = s.get("before_filtering", {})
        af = s.get("after_filtering", {})
        fr = d.get("filtering_result", {})
        passed = fr.get("passed_filter_reads", 0)
        total = bf.get("total_reads", 0)
        retention = (passed / total * 100) if total else 0
        rows.append({
            "Sample": sample,
            "Total reads (before)": bf.get("total_reads", 0),
            "Total reads (after)": af.get("total_reads", 0),
            "Retention %": round(retention, 1),
            "Mean length (bp)": af.get("read_mean_length", 0),
            "Q20 rate (%)": round(af.get("q20_rate", 0) * 100, 1),
            "Q30 rate (%)": round(af.get("q30_rate", 0) * 100, 1),
            "GC %": round(af.get("gc_content", 0) * 100, 2),
            "Low quality": fr.get("low_quality_reads", 0),
            "Too short": fr.get("too_short_reads", 0),
            "Too long": fr.get("too_long_reads", 0),
        })
    return pd.DataFrame(rows)


def downsample_quality_curves(curves: list, max_points: int = 2000) -> tuple[list, list]:
    """Downsample long quality curves for plotting."""
    n = len(curves)
    if n <= max_points:
        return list(range(n)), curves
    step = max(1, n // max_points)
    x = list(range(0, n, step))
    y = [curves[i] for i in x]
    return x, y


def generate_report(data: dict, output_path: Path, title: str = "fastplong QC Report") -> None:
    """Generate MultiQC-style HTML report with Plotly charts."""
    df = build_summary_df(data)
    if df.empty:
        print("No data to plot", file=sys.stderr)
        return

    samples = df["Sample"].tolist()
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    chart_divs = []

    # 1. Summary table
    table = go.Figure(data=[go.Table(
        header=dict(values=list(df.columns), fill_color="paleturquoise", align="left"),
        cells=dict(
            values=[df[c].tolist() for c in df.columns],
            fill_color="lavender",
            align="left",
        ),
    )])
    table.update_layout(
        title="Summary statistics per sample (sortable in full report)",
        margin=dict(l=20, r=20, t=50, b=20),
        height=min(400, 80 + len(df) * 28),
    )
    chart_divs.append(("summary_table", "Summary", table.to_html(full_html=False, include_plotlyjs=False)))

    # 2. Quality by position – all samples overlaid
    fig_qual = go.Figure()
    for sample in samples:
        d = data.get(sample, {})
        rf = d.get("read_after_filtering", d.get("read_before_filtering", {}))
        qc = rf.get("quality_curves", {})
        mean_q = qc.get("mean", [])
        if mean_q:
            x, y = downsample_quality_curves(mean_q)
            fig_qual.add_trace(
                go.Scatter(
                    x=x, y=y, mode="lines", name=sample,
                    line=dict(width=1.5),
                    hovertemplate="<b>%{fullData.name}</b><br>Position: %{x}<br>Mean Q: %{y:.1f}<extra></extra>",
                )
            )
    fig_qual.update_layout(
        title="Mean base quality by position (hover to identify sample)",
        xaxis_title="Position in read",
        yaxis_title="Mean quality (Phred)",
        margin=dict(l=60, r=20, t=50, b=50),
        height=450,
        hovermode="x unified",
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=1.02),
        showlegend=True,
    )
    chart_divs.append(("quality_curves", "Quality by position", fig_qual.to_html(full_html=False, include_plotlyjs=False)))

    # 3. Read count: before vs after
    fig_reads = go.Figure()
    fig_reads.add_trace(go.Bar(name="Before filtering", x=df["Sample"], y=df["Total reads (before)"], marker_color="lightblue"))
    fig_reads.add_trace(go.Bar(name="After filtering", x=df["Sample"], y=df["Total reads (after)"], marker_color="steelblue"))
    fig_reads.update_layout(
        title="Read count before vs after filtering",
        barmode="group",
        xaxis_tickangle=-45,
        margin=dict(l=60, r=20, t=50, b=100),
        height=450,
        hovermode="x unified",
    )
    chart_divs.append(("read_counts", "Read counts", fig_reads.to_html(full_html=False, include_plotlyjs=False)))

    # 4. Retention rate
    fig_ret = px.bar(df, x="Sample", y="Retention %", color="Retention %", color_continuous_scale="Blues")
    fig_ret.update_layout(
        title="Retention rate per sample (%)",
        xaxis_tickangle=-45,
        margin=dict(l=60, r=20, t=50, b=100),
        height=400,
        showlegend=False,
    )
    chart_divs.append(("retention", "Retention rate", fig_ret.to_html(full_html=False, include_plotlyjs=False)))

    # 5. Mean read length
    fig_len = px.bar(df, x="Sample", y="Mean length (bp)", color="Mean length (bp)", color_continuous_scale="Viridis")
    fig_len.update_layout(
        title="Mean read length after filtering",
        xaxis_tickangle=-45,
        margin=dict(l=60, r=20, t=50, b=100),
        height=400,
        showlegend=False,
    )
    chart_divs.append(("mean_length", "Mean read length", fig_len.to_html(full_html=False, include_plotlyjs=False)))

    # 6. Q20 / Q30 rate
    fig_q = go.Figure()
    fig_q.add_trace(go.Bar(name="Q20 rate (%)", x=df["Sample"], y=df["Q20 rate (%)"], marker_color="darkseagreen"))
    fig_q.add_trace(go.Bar(name="Q30 rate (%)", x=df["Sample"], y=df["Q30 rate (%)"], marker_color="seagreen"))
    fig_q.update_layout(
        title="Base quality rates per sample",
        barmode="group",
        xaxis_tickangle=-45,
        margin=dict(l=60, r=20, t=50, b=100),
        height=400,
        hovermode="x unified",
    )
    chart_divs.append(("quality_rates", "Q20/Q30 rates", fig_q.to_html(full_html=False, include_plotlyjs=False)))

    # 7. Scatter: batch effect detection
    fig_scatter = px.scatter(
        df, x="Mean length (bp)", y="Q30 rate (%)", text="Sample", hover_data=["Sample", "Retention %"]
    )
    fig_scatter.update_traces(textposition="top center", textfont_size=10)
    fig_scatter.update_layout(
        title="Mean read length vs Q30 rate (hover to identify outliers/batch effects)",
        margin=dict(l=60, r=20, t=50, b=50),
        height=450,
    )
    chart_divs.append(("scatter_length_q30", "Length vs quality (batch effect)", fig_scatter.to_html(full_html=False, include_plotlyjs=False)))

    # 8. Filtering breakdown
    fig_filt = go.Figure()
    fig_filt.add_trace(go.Bar(name="Passed", x=df["Sample"], y=df["Total reads (after)"], marker_color="forestgreen"))
    fig_filt.add_trace(go.Bar(name="Low quality", x=df["Sample"], y=df["Low quality"], marker_color="coral"))
    fig_filt.add_trace(go.Bar(name="Too short", x=df["Sample"], y=df["Too short"], marker_color="gold"))
    fig_filt.add_trace(go.Bar(name="Too long", x=df["Sample"], y=df["Too long"], marker_color="tomato"))
    fig_filt.update_layout(
        title="Filtering breakdown per sample",
        barmode="stack",
        xaxis_tickangle=-45,
        margin=dict(l=60, r=20, t=50, b=100),
        height=450,
        hovermode="x unified",
    )
    chart_divs.append(("filtering", "Filtering breakdown", fig_filt.to_html(full_html=False, include_plotlyjs=False)))

    # Embed Plotly.js inline so report works offline when transferred from remote servers
    plotly_js = _get_plotlyjs()

    # Build HTML
    nav_items = "".join(f'<li><a href="#{cid}">{label}</a></li>' for cid, label, _ in chart_divs)
    sections = "".join(
        f'<section id="{cid}" class="report-section"><h2>{label}</h2><div class="plot-container">{html}</div></section>'
        for cid, label, html in chart_divs
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>{title}</title>
    <script>{plotly_js}</script>
    <style>
        * {{ box-sizing: border-box; }}
        body {{ font-family: 'Helvetica Neue', Arial, sans-serif; margin: 0; background: #f5f5f5; color: #333; }}
        .header {{ background: linear-gradient(135deg, #1e88e5 0%, #1565c0 100%); color: white; padding: 24px 32px; }}
        .header h1 {{ margin: 0 0 8px 0; font-size: 1.8em; font-weight: 600; }}
        .header p {{ margin: 0; opacity: 0.9; font-size: 0.95em; }}
        .container {{ display: flex; max-width: 1600px; margin: 0 auto; }}
        .nav {{ width: 220px; flex-shrink: 0; background: white; padding: 20px; border-right: 1px solid #ddd; position: sticky; top: 0; height: 100vh; overflow-y: auto; }}
        .nav h3 {{ margin: 0 0 12px 0; font-size: 0.85em; color: #666; text-transform: uppercase; }}
        .nav ul {{ list-style: none; padding: 0; margin: 0; }}
        .nav li {{ margin: 4px 0; }}
        .nav a {{ color: #1565c0; text-decoration: none; font-size: 0.9em; display: block; padding: 6px 10px; border-radius: 4px; }}
        .nav a:hover {{ background: #e3f2fd; }}
        .content {{ flex: 1; padding: 24px 32px; }}
        .report-section {{ background: white; margin-bottom: 24px; padding: 24px; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }}
        .report-section h2 {{ margin: 0 0 20px 0; font-size: 1.2em; color: #1565c0; border-bottom: 2px solid #e3f2fd; padding-bottom: 10px; }}
        .plot-container {{ width: 100%; min-height: 200px; }}
        .plot-container .plotly {{ width: 100% !important; }}
        .footer {{ text-align: center; padding: 20px; color: #666; font-size: 0.85em; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>{title}</h1>
        <p>Aggregated fastplong QC across {len(samples)} samples · Generated {timestamp}</p>
    </div>
    <div class="container">
        <nav class="nav">
            <h3>Report sections</h3>
            <ul>
                {nav_items}
            </ul>
        </nav>
        <main class="content">
            {sections}
        </main>
    </div>
    <div class="footer">
        fastplong_multireport · Aggregated fastplong QC
    </div>
</body>
</html>
"""

    with open(output_path, "w") as f:
        f.write(html)
    print(f"Report written: {output_path}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate fastplong QC reports from a results directory into a MultiQC-style HTML report.",
        epilog="Example: fastplong_multireport results/qc/raw/ -o fastplong_report.html",
    )
    parser.add_argument(
        "results_dir",
        type=Path,
        help="Directory containing fastplong JSON reports (*_fastplong_report.json)",
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=None,
        help="Output HTML path (default: fastplong_multireport.html in results_dir)",
    )
    parser.add_argument(
        "--no-recursive",
        action="store_true",
        help="Only search the top-level of results_dir (default: search recursively)",
    )
    parser.add_argument(
        "-t", "--title",
        default="fastplong QC Report",
        help="Report title (default: fastplong QC Report)",
    )
    args = parser.parse_args()

    results_dir = args.results_dir.resolve()
    if not results_dir.is_dir():
        print(f"Error: {results_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    output_path = args.output
    if output_path is None:
        output_path = results_dir / "fastplong_multireport.html"
    else:
        output_path = output_path.resolve()

    json_files = discover_json_reports(results_dir, recursive=not args.no_recursive)
    if not json_files:
        print(f"Error: No *_fastplong_report.json files found in {results_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(json_files)} fastplong report(s)", file=sys.stderr)

    data = load_fastplong_reports(json_files)
    if not data:
        print("Error: Could not load any fastplong reports", file=sys.stderr)
        sys.exit(1)

    generate_report(data, output_path, title=args.title)


def cli() -> None:
    """CLI entry point."""
    main()


if __name__ == "__main__":
    cli()
