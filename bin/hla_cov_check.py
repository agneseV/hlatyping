#!/usr/bin/env python3
"""Summarize exon depths from SAMTOOLS_BEDCOV into one HLA coverage note."""
import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

CLASS_I = ['HLA-A', 'HLA-B', 'HLA-C']
CLASS_II = [
    'HLA-DRA', 'HLA-DRB5', 'HLA-DRB1', 'HLA-DQA1', 'HLA-DQB1',
    'HLA-DQA2', 'HLA-DQB2', 'HLA-DPA1', 'HLA-DPB1',
]
GENES = CLASS_I + CLASS_II
GENE_RE = re.compile(r'(?:^|;)gene:([^;]+)(?:;|$)')


def check_coverage(input_path, threshold):
    depths = defaultdict(list)
    seen = set()
    with open(input_path, newline='') as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter='\t'), 1):
            if not row or row[0].startswith('#'):
                continue
            if row[:3] == ['chrom', 'start', 'end']:
                continue
            if len(row) != 6:
                raise ValueError(f'line {line_number}: expected six columns, got {len(row)}')
            match = GENE_RE.search(row[3])
            if not match:
                raise ValueError(f'line {line_number}: no gene:<name> in region label')
            gene = match.group(1)
            if gene not in GENES:
                raise ValueError(f'line {line_number}: unexpected gene {gene}')
            try:
                start, end = int(row[1]), int(row[2])
                depth = float(row[5])
            except ValueError as exc:
                raise ValueError(f'line {line_number}: invalid coordinate or mean depth') from exc
            if start < 0 or end <= start or not math.isfinite(depth) or depth < 0:
                raise ValueError(f'line {line_number}: invalid interval or mean depth')
            key = (gene, row[0], start, end)
            if key in seen:
                raise ValueError(f'line {line_number}: duplicate exon interval for {gene}: {key[1:]}')
            seen.add(key)
            depths[gene].append(depth)

    missing = [gene for gene in GENES if not depths[gene]]
    if missing:
        raise ValueError('missing exon coverage for: ' + ', '.join(missing))

    means = {gene: sum(depths[gene]) / len(depths[gene]) for gene in GENES}
    low = {gene for gene in GENES if means[gene] < threshold}
    if len(low) == len(GENES):
        return 'LowHLA', means, depths
    labels = []
    for group, group_label in ((CLASS_I, 'LowHLA-ClassI'), (CLASS_II, 'LowHLA-ClassII')):
        if all(gene in low for gene in group):
            labels.append(group_label)
        else:
            labels.extend(f'Low{gene}' for gene in group if gene in low)
    return ';'.join(labels) if labels else 'OK', means, depths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('coverage', type=Path, help='SAMTOOLS_BEDCOV *.exon_coverage.tsv')
    parser.add_argument('--threshold', type=float, required=True, help='Minimum acceptable mean exon depth (x)')
    parser.add_argument('--output', type=Path, required=True, help='Single-line coverage note')
    parser.add_argument('--summary', type=Path, required=True, help='Per-gene mean-depth table')
    args = parser.parse_args()
    if not math.isfinite(args.threshold) or args.threshold < 0:
        parser.error('--threshold must be a finite, nonnegative number')
    try:
        note, means, depths = check_coverage(args.coverage, args.threshold)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    args.output.write_text(note + '\n')
    with args.summary.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['gene', 'class', 'exons', 'mean_exon_depth', 'threshold', 'status'])
        for gene in GENES:
            writer.writerow([
                gene, 'I' if gene in CLASS_I else 'II', len(depths[gene]),
                f'{means[gene]:.6f}', args.threshold,
                'Low' if means[gene] < args.threshold else 'OK',
            ])


if __name__ == '__main__':
    main()
