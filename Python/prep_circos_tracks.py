#!/usr/bin/env python3

import os
import re
import gzip
from collections import defaultdict


# ============================================================
# settings
# ============================================================

FASTA = (
    '/home/yshin/mendel-nas1/snake_genome_ass/'
    'G_ussuriensis_Chromo/annotation/soft_masked/'
    'Gloydius_ussuriensis_EarlGrey/'
    'Gloydius_ussuriensis_summaryFiles/'
    'Gloydius_ussuriensis.softmasked.fasta'
)

GFF = (
    '/home/yshin/mendel-nas1/snake_genome_ass/'
    'G_ussuriensis_Chromo/annotation/funannotate/'
    'annotate_results/'
    'Gloydius_ussuriensis_AMNH_21010.gff3'
)

OUTDIR = (
    '/home/yshin/mendel-nas1/snake_genome_ass/'
    'G_ussuriensis_Chromo/circos/02_tracks'
)

# 100 kb non-overlapping windows
WINDOW = 100_000

os.makedirs(
    OUTDIR,
    exist_ok=True
)


# ============================================================
# helper functions
# ============================================================

def open_text(path):

    if path.endswith('.gz'):

        return gzip.open(
            path,
            'rt'
        )

    return open(
        path,
        'r'
    )


def fasta_reader(path):

    with open_text(path) as handle:

        name = None
        sequence = []

        for line in handle:

            line = line.rstrip()

            if not line:

                continue

            if line.startswith('>'):

                if name is not None:

                    yield name, ''.join(sequence)

                name = line[1:].split()[0]

                sequence = []

            else:

                sequence.append(line)

        if name is not None:

            yield name, ''.join(sequence)


def chromosome_label(seqid):

    '''
    Recognizes both:

        G_ussuri_chr1
        chr1

        G_ussuri_chrZ
        chrZ

        G_ussuri_chrW
        chrW
    '''

    match = re.search(
        r'(?:^|_)chr([0-9]+|Z|W)$',
        seqid
    )

    if match is None:

        return None

    return match.group(1)


def chromosome_sort_key(seqid):

    label = chromosome_label(seqid)

    if label is None:

        return (
            9999,
            seqid
        )

    if label.isdigit():

        return (
            int(label),
            ''
        )

    if label == 'Z':

        return (
            1000,
            'Z'
        )

    if label == 'W':

        return (
            1001,
            'W'
        )

    return (
        9999,
        label
    )


def chromosome_class(seqid):

    label = chromosome_label(seqid)

    if label is None:

        return 'Unplaced'

    if label in (
        'Z',
        'W'
    ):

        return 'Sex chromosome'

    if label.isdigit():

        number = int(label)

        if number <= 7:

            return 'Macrochromosome'

        return 'Microchromosome'

    return 'Unplaced'


# ============================================================
# output files
# ============================================================

LENGTH_OUT = os.path.join(
    OUTDIR,
    'G_ussuriensis.chromosome_lengths.tsv'
)

CLASS_OUT = os.path.join(
    OUTDIR,
    'G_ussuriensis.chromosome_classes.tsv'
)

TRACK_OUT = os.path.join(
    OUTDIR,
    'G_ussuriensis.100kb.circos_tracks.tsv'
)


# ============================================================
# calculate GC % and Repeat %
# ============================================================

chromosome_lengths = {}

track_rows = []


print()

print(
    'Reading genome and calculating '
    '100 kb GC/repeat tracks...'
)

print()


for seqid, seq in fasta_reader(FASTA):

    label = chromosome_label(seqid)

    # exclude unplaced scaffolds
    if label is None:

        continue

    length = len(seq)

    chromosome_lengths[seqid] = length

    print(
        f'{seqid}: '
        f'{length:,} bp'
    )


    # --------------------------------------------------------
    # divide chromosome into 100 kb windows
    # --------------------------------------------------------

    for start in range(
        0,
        length,
        WINDOW
    ):

        end = min(
            start + WINDOW,
            length
        )

        fragment = seq[
            start:end
        ]


        # ----------------------------------------------------
        # GC %
        #
        # denominator includes only A/C/G/T bases
        # Ns are therefore ignored
        # ----------------------------------------------------

        upper = fragment.upper()

        A = upper.count('A')

        C = upper.count('C')

        G = upper.count('G')

        T = upper.count('T')

        acgt = (
            A +
            C +
            G +
            T
        )

        if acgt > 0:

            gc_pct = (
                (G + C) /
                acgt *
                100.0
            )

        else:

            gc_pct = float(
                'nan'
            )


        # ----------------------------------------------------
        # Repeat %
        #
        # EarlGrey soft-masked FASTA stores repetitive bases
        # as lowercase sequence
        #
        # denominator includes A/C/G/T only
        # ----------------------------------------------------

        repeat_bases = (
            fragment.count('a') +
            fragment.count('c') +
            fragment.count('g') +
            fragment.count('t')
        )

        if acgt > 0:

            repeat_pct = (
                repeat_bases /
                acgt *
                100.0
            )

        else:

            repeat_pct = float(
                'nan'
            )


        # ----------------------------------------------------
        # initialize track row
        #
        # gene_count will be filled later
        # ----------------------------------------------------

        track_rows.append(
            {
                'chr': seqid,
                'start': start,
                'end': end,
                'GC_pct': gc_pct,
                'Repeat_pct': repeat_pct,
                'gene_count': 0
            }
        )


# ============================================================
# sort chromosome order
# ============================================================

chromosome_order = sorted(
    chromosome_lengths,
    key=chromosome_sort_key
)


chromosome_rank = {
    chrom: i
    for i, chrom in enumerate(
        chromosome_order
    )
}


# ============================================================
# count genes per 100 kb window
#
# each gene is counted exactly once using its midpoint
#
# this avoids counting long genes multiple times when a gene
# spans more than one 100 kb window
# ============================================================

print()

print(
    'Reading Funannotate GFF3 and calculating '
    '100 kb gene density...'
)

print()


gene_counts = defaultdict(
    int
)

total_genes = 0

chromosomal_genes = 0

excluded_genes = 0


with open_text(GFF) as handle:

    for line in handle:

        if line.startswith('#'):

            continue

        fields = line.rstrip(
            '\n'
        ).split(
            '\t'
        )

        if len(fields) != 9:

            continue

        seqid = fields[0]

        feature = fields[2]


        # ----------------------------------------------------
        # count gene features only
        # ----------------------------------------------------

        if feature != 'gene':

            continue


        total_genes += 1


        # ----------------------------------------------------
        # exclude genes not on retained chromosomes
        # ----------------------------------------------------

        if seqid not in chromosome_lengths:

            excluded_genes += 1

            continue


        # ----------------------------------------------------
        # GFF3 coordinates are 1-based inclusive
        #
        # convert start to 0-based before calculating midpoint
        # ----------------------------------------------------

        start = (
            int(fields[3]) -
            1
        )

        end = int(
            fields[4]
        )


        midpoint = (
            start +
            end
        ) // 2


        # ----------------------------------------------------
        # identify 100 kb window containing gene midpoint
        # ----------------------------------------------------

        window_start = (
            midpoint //
            WINDOW
        ) * WINDOW


        gene_counts[
            (
                seqid,
                window_start
            )
        ] += 1


        chromosomal_genes += 1


# ============================================================
# add gene counts to track table
# ============================================================

for row in track_rows:

    row['gene_count'] = gene_counts[
        (
            row['chr'],
            row['start']
        )
    ]


# ============================================================
# sort track rows
# ============================================================

track_rows.sort(
    key=lambda x: (
        chromosome_rank[
            x['chr']
        ],
        x['start']
    )
)


# ============================================================
# write chromosome lengths
# ============================================================

with open(
    LENGTH_OUT,
    'w'
) as out:

    out.write(
        'chr\tlength\n'
    )

    for chrom in chromosome_order:

        out.write(
            f'{chrom}\t'
            f'{chromosome_lengths[chrom]}\n'
        )


# ============================================================
# write chromosome classes
# ============================================================

with open(
    CLASS_OUT,
    'w'
) as out:

    out.write(
        'chr\tdisplay\tclass\n'
    )

    for chrom in chromosome_order:

        label = chromosome_label(
            chrom
        )

        cls = chromosome_class(
            chrom
        )

        out.write(
            f'{chrom}\t'
            f'Chr {label}\t'
            f'{cls}\n'
        )


# ============================================================
# write combined 100 kb track file
# ============================================================

with open(
    TRACK_OUT,
    'w'
) as out:

    out.write(
        'chr\t'
        'start\t'
        'end\t'
        'GC_pct\t'
        'Repeat_pct\t'
        'gene_count\n'
    )

    for row in track_rows:

        out.write(
            f'{row["chr"]}\t'
            f'{row["start"]}\t'
            f'{row["end"]}\t'
            f'{row["GC_pct"]:.6f}\t'
            f'{row["Repeat_pct"]:.6f}\t'
            f'{row["gene_count"]}\n'
        )


# ============================================================
# calculate summary statistics
# ============================================================

total_windows = len(
    track_rows
)

total_length = sum(
    chromosome_lengths.values()
)

mean_gc = sum(
    row['GC_pct']
    for row in track_rows
    if row['GC_pct'] == row['GC_pct']
) / sum(
    1
    for row in track_rows
    if row['GC_pct'] == row['GC_pct']
)

mean_repeat = sum(
    row['Repeat_pct']
    for row in track_rows
    if row['Repeat_pct'] == row['Repeat_pct']
) / sum(
    1
    for row in track_rows
    if row['Repeat_pct'] == row['Repeat_pct']
)


# ============================================================
# report
# ============================================================

print()

print(
    '============================================================'
)

print(
    'Finished'
)

print(
    '============================================================'
)

print()


print(
    f'Window size: '
    f'{WINDOW:,} bp'
)

print(
    f'Total windows: '
    f'{total_windows:,}'
)

print(
    f'Chromosomes retained: '
    f'{len(chromosome_order)}'
)

print(
    f'Total plotted sequence: '
    f'{total_length:,} bp'
)

print()


print(
    f'Genes in GFF3: '
    f'{total_genes:,}'
)

print(
    f'Genes on plotted chromosomes: '
    f'{chromosomal_genes:,}'
)

print(
    f'Genes excluded with unplaced scaffolds: '
    f'{excluded_genes:,}'
)

print()


print(
    f'Mean window GC %: '
    f'{mean_gc:.3f}'
)

print(
    f'Mean window Repeat %: '
    f'{mean_repeat:.3f}'
)

print()


print(
    'Chromosomes:'
)

for chrom in chromosome_order:

    print(
        f'  {chrom:25s} '
        f'{chromosome_lengths[chrom]:>12,} bp  '
        f'{chromosome_class(chrom)}'
    )


print()

print(
    'Output files:'
)

print(
    LENGTH_OUT
)

print(
    CLASS_OUT
)

print(
    TRACK_OUT
)

print()