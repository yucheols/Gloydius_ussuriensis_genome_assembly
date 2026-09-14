#!/usr/bin/env python3

# ============================================================
# Build chromosome manifest for GENESPACE
#
# Uses:
#   - NCBI sequence_report.jsonl
#   - chromosome-labelled genome FASTA headers
#   - G. ussuriensis final annotation chromosome IDs
#
# Output:
#   GENESPACE/QC/chromosome_manifest.tsv
#
# ============================================================


import csv
import glob
import json
import re
from pathlib import Path


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------

SYN = Path(
    '/home/yshin/mendel-nas1/snake_genome_ass/'
    'G_ussuriensis_Chromo/synteny'
)

BASE = SYN / 'assemblies_synteny'
GS = SYN / 'GENESPACE'
QC = GS / 'QC'


# ------------------------------------------------------------
# make required directories
# ------------------------------------------------------------

for directory in [
    GS,
    QC,
    GS / 'bed',
    GS / 'peptide',
    GS / 'mapping'
]:

    directory.mkdir(
        parents=True,
        exist_ok=True
    )


# ------------------------------------------------------------
# G. ussuriensis final annotation
# ------------------------------------------------------------

GUSS_GFF = Path(
    '/home/yshin/mendel-nas1/snake_genome_ass/'
    'G_ussuriensis_Chromo/annotation/funannotate/'
    'annotate_results/'
    'Gloydius_ussuriensis_AMNH_21010.gff3'
)


# ------------------------------------------------------------
# species in clean GENESPACE analysis
#
# Crotalus viridis remains excluded because its source
# annotation requires separate duplicated-model cleanup.
# ------------------------------------------------------------

species_order = [
    'Argyrophis_diardii',
    'Xenopeltis_unicolor',
    'Candoia_aspera',
    'Elaphe_schrenckii',
    'Naja_naja',
    'Cerastes_gasperettii',
    'Vipera_berus',
    'Bothrops_insularis',
    'Crotalus_adamanteus',
    'Gloydius_shedaoensis',
    'Gloydius_ussuriensis'
]


# ------------------------------------------------------------
# source annotations
#
# Only root/source annotation files are used.
#
# Files under fasta_gff_ids_check/ are deliberately ignored.
# ------------------------------------------------------------

annotations = {

    'Argyrophis_diardii':
        BASE / 'Argyrophis_diardii' /
        'Argyrophis_diardii.annotation.gff3',

    'Bothrops_insularis':
        BASE / 'Bothrops_insularis' /
        'Bothrops_insularis.annotation.gtf',

    'Candoia_aspera':
        BASE / 'Candoia_aspera' /
        'Candoia_aspera.annotation.gff3',

    'Cerastes_gasperettii':
        BASE / 'Cerastes_gasperettii' /
        'Cerastes_gasperettii.annotation.gff3',

    'Crotalus_adamanteus':
        BASE / 'Crotalus_adamanteus' /
        'Crotalus_adamanteus.annotation.gff3',

    'Elaphe_schrenckii':
        BASE / 'Elaphe_schrenckii' /
        'Elaphe_schrenckii.annotation.gff3',

    'Gloydius_shedaoensis':
        BASE / 'Gloydius_shedaoensis' /
        'Gloydius_shedaoensis.annotation.gff3',

    'Naja_naja':
        BASE / 'Naja_naja' /
        'Naja_naja.annotation.gff3',

    'Vipera_berus':
        BASE / 'Vipera_berus' /
        'Vipera_berus.annotation.gff3',

    'Xenopeltis_unicolor':
        BASE / 'Xenopeltis_unicolor' /
        'Xenopeltis_unicolor.annotation.gff3',

    'Gloydius_ussuriensis':
        GUSS_GFF
}


# ------------------------------------------------------------
# NCBI sequence-report species
# ------------------------------------------------------------

ncbi_species = [
    'Bothrops_insularis',
    'Candoia_aspera',
    'Cerastes_gasperettii',
    'Crotalus_adamanteus',
    'Elaphe_schrenckii',
    'Naja_naja',
    'Vipera_berus'
]


# ------------------------------------------------------------
# species whose FASTA headers explicitly identify chromosomes
# ------------------------------------------------------------

fasta_species = {

    'Argyrophis_diardii':
        BASE / 'Argyrophis_diardii' /
        'Argyrophis_diardii.genome.fa',

    'Gloydius_shedaoensis':
        BASE / 'Gloydius_shedaoensis' /
        'Gloydius_shedaoensis.genome.fa',

    'Xenopeltis_unicolor':
        BASE / 'Xenopeltis_unicolor' /
        'Xenopeltis_unicolor.genome.fa'
}


# ------------------------------------------------------------
# read sequence IDs represented in a source annotation
# ------------------------------------------------------------

def annotation_seqids(path):

    ids = set()

    with path.open() as handle:

        for line in handle:

            if not line.strip():
                continue

            if line.startswith('#'):
                continue

            fields = line.rstrip('\n').split('\t')

            if len(fields) < 9:
                continue

            ids.add(
                fields[0]
            )

    return ids


# ------------------------------------------------------------
# standardize chromosome label formatting
#
# This DOES NOT infer chromosome homology.
#
# Examples:
#
# 1           -> chr1
# 01          -> chr1
# Chromosome1 -> chr1
# Z           -> chrZ
# W           -> chrW
# ------------------------------------------------------------

def standard_chr_name(x):

    x = str(x).strip()

    x = re.sub(
        r'^chromosome[_\s-]*',
        '',
        x,
        flags=re.IGNORECASE
    )

    x = re.sub(
        r'^chr',
        '',
        x,
        flags=re.IGNORECASE
    )

    if re.fullmatch(
        r'0*[0-9]+',
        x
    ):

        return f'chr{int(x)}'

    if x.upper() in {
        'Z',
        'W'
    }:

        return f'chr{x.upper()}'

    return f'chr{x}'


# ------------------------------------------------------------
# verify source annotation files exist
# ------------------------------------------------------------

print()
print('========================================')
print('SOURCE ANNOTATION CHECK')
print('========================================')

for species in species_order:

    path = annotations[species]

    if not path.exists():

        raise FileNotFoundError(
            f'Missing annotation for {species}: {path}'
        )

    print(
        f'{species:25s} PASS  {path}'
    )


# ------------------------------------------------------------
# collect annotation sequence IDs
# ------------------------------------------------------------

annotation_ids = {}

for species in species_order:

    annotation_ids[species] = annotation_seqids(
        annotations[species]
    )


# ------------------------------------------------------------
# manifest rows
# ------------------------------------------------------------

rows = []


# ============================================================
# A. NCBI assemblies
# ============================================================

for species in ncbi_species:

    print()
    print('========================================')
    print(species)
    print('========================================')

    pattern = str(
        BASE /
        species /
        'ncbi_dataset' /
        'data' /
        '*' /
        'sequence_report.jsonl'
    )

    reports = sorted(
        glob.glob(pattern)
    )

    if len(reports) != 1:

        raise RuntimeError(
            f'Expected exactly one sequence_report.jsonl '
            f'for {species}, found {len(reports)}'
        )

    report = Path(
        reports[0]
    )

    print(
        'sequence report:',
        report
    )

    ann_ids = annotation_ids[
        species
    ]

    chromosome_records = 0

    with report.open() as handle:

        for line in handle:

            if not line.strip():
                continue

            rec = json.loads(
                line
            )

            role = str(
                rec.get(
                    'role',
                    ''
                )
            ).strip()

            molecule_type = str(
                rec.get(
                    'assignedMoleculeLocationType',
                    ''
                )
            ).strip()

            chr_name = str(
                rec.get(
                    'chrName',
                    ''
                )
            ).strip()


            # --------------------------------------------
            # retain bona fide assembled chromosome records
            # only
            # --------------------------------------------

            if role != 'assembled-molecule':
                continue

            if molecule_type.lower() != 'chromosome':
                continue

            if chr_name == '':
                continue


            chromosome_records += 1


            # --------------------------------------------
            # accession/name fields
            # --------------------------------------------

            refseq = str(
                rec.get(
                    'refseqAccession',
                    ''
                )
            ).strip()

            genbank = str(
                rec.get(
                    'genbankAccession',
                    ''
                )
            ).strip()

            sequence_name = str(
                rec.get(
                    'sequenceName',
                    ''
                )
            ).strip()


            # --------------------------------------------
            # determine which identifier is actually used
            # by the source annotation
            # --------------------------------------------

            candidates = [
                refseq,
                genbank,
                sequence_name
            ]

            candidates = [
                x
                for x in candidates
                if x != ''
            ]

            matches = [
                x
                for x in candidates
                if x in ann_ids
            ]

            matches = list(
                dict.fromkeys(
                    matches
                )
            )


            if len(matches) == 1:

                annotation_seqid = matches[0]

                status = 'PASS'


            elif len(matches) == 0:

                annotation_seqid = ''

                status = 'NO_ANNOTATION_MATCH'


            else:

                annotation_seqid = '|'.join(
                    matches
                )

                status = 'AMBIGUOUS_MATCH'


            # --------------------------------------------
            # standardized chromosome label
            # --------------------------------------------

            std_chr = standard_chr_name(
                chr_name
            )


            # --------------------------------------------
            # sex chromosome class
            #
            # ONLY explicit source metadata is used.
            #
            # Numerically labelled chromosomes remain
            # unknown/autosome until later evaluation.
            # --------------------------------------------

            if std_chr == 'chrW':

                sex_class = 'W'

            elif std_chr == 'chrZ':

                sex_class = 'Z'

            else:

                sex_class = 'unknown_or_autosome'


            rows.append(
                {
                    'species':
                        species,

                    'annotation_seqid':
                        annotation_seqid,

                    'source_chr':
                        chr_name,

                    'standard_chr':
                        std_chr,

                    'sex_class':
                        sex_class,

                    'length':
                        rec.get(
                            'length',
                            ''
                        ),

                    'role':
                        role,

                    'molecule_type':
                        molecule_type,

                    'refseq_accession':
                        refseq,

                    'genbank_accession':
                        genbank,

                    'sequence_name':
                        sequence_name,

                    'metadata_source':
                        'NCBI_sequence_report',

                    'annotation_match':
                        status
                }
            )


    print(
        'chromosome records:',
        chromosome_records
    )


# ============================================================
# B. FASTA-header assemblies
# ============================================================

for species, fasta in fasta_species.items():

    print()
    print('========================================')
    print(species)
    print('========================================')

    if not fasta.exists():

        raise FileNotFoundError(
            f'Missing genome FASTA for {species}: {fasta}'
        )

    ann_ids = annotation_ids[
        species
    ]

    chromosome_records = 0


    with fasta.open() as handle:

        for line in handle:

            if not line.startswith('>'):
                continue

            header = line.rstrip(
                '\n'
            )


            # --------------------------------------------
            # first token after >
            # --------------------------------------------

            seqid = header[
                1:
            ].split()[0]


            # --------------------------------------------
            # only records explicitly marked Chromosome
            # are retained
            # --------------------------------------------

            m_chr = re.search(
                r'\bChromosome\s+([^\s]+)',
                header,
                flags=re.IGNORECASE
            )

            if m_chr is None:
                continue


            chr_name = m_chr.group(
                1
            )

            chromosome_records += 1


            # --------------------------------------------
            # extract sequence length from header
            # --------------------------------------------

            m_len = re.search(
                r'\bLen=([0-9]+)',
                header
            )

            if m_len:

                length = int(
                    m_len.group(
                        1
                    )
                )

            else:

                length = ''


            # --------------------------------------------
            # annotation ID match
            # --------------------------------------------

            if seqid in ann_ids:

                status = 'PASS'

            else:

                status = 'NO_ANNOTATION_MATCH'


            std_chr = standard_chr_name(
                chr_name
            )


            if std_chr == 'chrW':

                sex_class = 'W'

            elif std_chr == 'chrZ':

                sex_class = 'Z'

            else:

                sex_class = 'unknown_or_autosome'


            rows.append(
                {
                    'species':
                        species,

                    'annotation_seqid':
                        seqid,

                    'source_chr':
                        chr_name,

                    'standard_chr':
                        std_chr,

                    'sex_class':
                        sex_class,

                    'length':
                        length,

                    'role':
                        'assembled-molecule',

                    'molecule_type':
                        'Chromosome',

                    'refseq_accession':
                        '',

                    'genbank_accession':
                        '',

                    'sequence_name':
                        '',

                    'metadata_source':
                        'FASTA_header',

                    'annotation_match':
                        status
                }
            )


    print(
        'chromosome records:',
        chromosome_records
    )


# ============================================================
# C. Gloydius ussuriensis
# ============================================================

species = 'Gloydius_ussuriensis'

print()
print('========================================')
print(species)
print('========================================')


ann_ids = annotation_ids[
    species
]

gussuri_records = []


for seqid in ann_ids:

    m = re.fullmatch(
        r'G_ussuri_chr(.+)',
        seqid
    )

    if m is None:
        continue


    chr_name = m.group(
        1
    )

    std_chr = standard_chr_name(
        chr_name
    )


    if std_chr == 'chrW':

        sex_class = 'W'

    elif std_chr == 'chrZ':

        sex_class = 'Z'

    else:

        sex_class = 'unknown_or_autosome'


    gussuri_records.append(
        {
            'species':
                species,

            'annotation_seqid':
                seqid,

            'source_chr':
                chr_name,

            'standard_chr':
                std_chr,

            'sex_class':
                sex_class,

            'length':
                '',

            'role':
                'assembled-molecule',

            'molecule_type':
                'Chromosome',

            'refseq_accession':
                '',

            'genbank_accession':
                '',

            'sequence_name':
                seqid,

            'metadata_source':
                'G_ussuriensis_annotation',

            'annotation_match':
                'PASS'
        }
    )


rows.extend(
    gussuri_records
)


print(
    'chromosome records:',
    len(
        gussuri_records
    )
)


# ------------------------------------------------------------
# chromosome sorting
# ------------------------------------------------------------

species_rank = {
    species: i
    for i, species in enumerate(
        species_order
    )
}


def chromosome_rank(row):

    chromosome = row[
        'standard_chr'
    ]

    m = re.fullmatch(
        r'chr([0-9]+)',
        chromosome
    )

    if m:

        return (
            species_rank[
                row['species']
            ],
            0,
            int(
                m.group(
                    1
                )
            )
        )


    if chromosome == 'chrZ':

        return (
            species_rank[
                row['species']
            ],
            1,
            0
        )


    if chromosome == 'chrW':

        return (
            species_rank[
                row['species']
            ],
            2,
            0
        )


    return (
        species_rank[
            row['species']
        ],
        3,
        chromosome
    )


rows.sort(
    key=chromosome_rank
)


# ------------------------------------------------------------
# write chromosome manifest
# ------------------------------------------------------------

manifest = (
    QC /
    'chromosome_manifest.tsv'
)


fields = [
    'species',
    'annotation_seqid',
    'source_chr',
    'standard_chr',
    'sex_class',
    'length',
    'role',
    'molecule_type',
    'refseq_accession',
    'genbank_accession',
    'sequence_name',
    'metadata_source',
    'annotation_match'
]


with manifest.open(
    'w',
    newline=''
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter='\t'
    )

    writer.writeheader()

    writer.writerows(
        rows
    )


# ============================================================
# QC summary
# ============================================================

print()
print('========================================')
print('CHROMOSOME MANIFEST SUMMARY')
print('========================================')


problem_count = 0


for species in species_order:

    sp_rows = [
        row
        for row in rows
        if row['species'] == species
    ]

    passes = [
        row
        for row in sp_rows
        if row['annotation_match'] == 'PASS'
    ]

    problems = [
        row
        for row in sp_rows
        if row['annotation_match'] != 'PASS'
    ]

    explicit_z = [
        row
        for row in sp_rows
        if row['sex_class'] == 'Z'
    ]

    explicit_w = [
        row
        for row in sp_rows
        if row['sex_class'] == 'W'
    ]


    problem_count += len(
        problems
    )


    print(
        f'{species:25s} '
        f'chromosomes={len(sp_rows):3d} '
        f'matched={len(passes):3d} '
        f'problems={len(problems):3d} '
        f'explicit_Z={len(explicit_z):2d} '
        f'explicit_W={len(explicit_w):2d}'
    )


# ------------------------------------------------------------
# duplicated standardized chromosome labels
# ------------------------------------------------------------

duplicates = []


for species in species_order:

    seen = {}

    for row in rows:

        if row['species'] != species:
            continue

        chromosome = row[
            'standard_chr'
        ]

        seen.setdefault(
            chromosome,
            []
        ).append(
            row[
                'annotation_seqid'
            ]
        )


    for chromosome, seqids in seen.items():

        if len(seqids) > 1:

            duplicates.append(
                (
                    species,
                    chromosome,
                    seqids
                )
            )


print()
print('========================================')
print('DUPLICATED STANDARDIZED CHROMOSOMES')
print('========================================')


if not duplicates:

    print(
        'none'
    )

else:

    for species, chromosome, seqids in duplicates:

        print(
            species,
            chromosome,
            ','.join(
                seqids
            )
        )


# ------------------------------------------------------------
# show unmatched chromosome records
# ------------------------------------------------------------

problems = [
    row
    for row in rows
    if row['annotation_match'] != 'PASS'
]


print()
print('========================================')
print('UNMATCHED / AMBIGUOUS CHROMOSOMES')
print('========================================')


if not problems:

    print(
        'none'
    )

else:

    for row in problems:

        print(
            row['species'],
            row['source_chr'],
            row['standard_chr'],
            row['annotation_match'],
            row['refseq_accession'],
            row['genbank_accession'],
            row['sequence_name'],
            sep='\t'
        )


# ------------------------------------------------------------
# final report
# ------------------------------------------------------------

print()
print('========================================')
print('FINAL STATUS')
print('========================================')

print()
print(
    'Manifest written to:'
)

print(
    manifest
)

print()


if problem_count == 0:

    print(
        'ALL chromosome records matched source annotations: PASS'
    )

else:

    print(
        f'WARNING: {problem_count} chromosome records '
        'did not match annotation sequence IDs.'
    )


print()

print(
    'Numeric chromosomes have NOT been classified '
    'as autosome, Z, or W.'
)

print(
    'Only explicit chromosome metadata has been used.'
)

print(
    'No previous GENESPACE results were used.'
)