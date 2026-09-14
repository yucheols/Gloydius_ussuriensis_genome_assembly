#!/usr/bin/env python3

# ============================================================
# Build chromosome-only GENESPACE inputs
#
# Uses:
#   - source annotations
#   - source protein FASTAs
#   - GENESPACE/QC/chromosome_manifest.tsv
#
# Rules:
#   - only chromosome-scale pseudomolecules are retained
#   - explicitly annotated W chromosomes are excluded
#   - Z chromosomes are retained
#   - numerically labelled chromosomes are retained for now
#   - longest protein isoform per biological gene/locus
#   - BED chromosome names standardized to chr1, chr2, ..., chrZ
#   - BED column 4 == peptide FASTA ID
#   - "." in protein sequences replaced with "X"
#
# Crotalus viridis is deliberately excluded.
#
# Outputs:
#   GENESPACE/bed/<species>.bed
#   GENESPACE/peptide/<species>.fa
#   GENESPACE/mapping/<species>.longest_isoform.tsv
#   GENESPACE/QC/chromosome_input_summary.tsv
# ============================================================


import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------

ROOT = Path(
    '/home/yshin/mendel-nas1/snake_genome_ass/'
    'G_ussuriensis_Chromo'
)

SYN = ROOT / 'synteny'

BASE = SYN / 'assemblies_synteny'

GS = SYN / 'GENESPACE'

BED_DIR = GS / 'bed'
PEP_DIR = GS / 'peptide'
MAP_DIR = GS / 'mapping'
QC_DIR = GS / 'QC'

MANIFEST = (
    QC_DIR /
    'chromosome_manifest.tsv'
)


for directory in [
    BED_DIR,
    PEP_DIR,
    MAP_DIR,
    QC_DIR
]:

    directory.mkdir(
        parents=True,
        exist_ok=True
    )


# ------------------------------------------------------------
# G. ussuriensis final annotation
# ------------------------------------------------------------

GUSS_DIR = (
    ROOT /
    'annotation' /
    'funannotate' /
    'annotate_results'
)


# ------------------------------------------------------------
# species used in this GENESPACE run
#
# IMPORTANT:
# Crotalus viridis is intentionally NOT included.
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
# source annotation + protein files
# ------------------------------------------------------------

sources = {

    'Argyrophis_diardii': {
        'annotation':
            BASE /
            'Argyrophis_diardii' /
            'Argyrophis_diardii.annotation.gff3',

        'protein':
            BASE /
            'Argyrophis_diardii' /
            'Argyrophis_diardii.protein.faa'
    },

    'Xenopeltis_unicolor': {
        'annotation':
            BASE /
            'Xenopeltis_unicolor' /
            'Xenopeltis_unicolor.annotation.gff3',

        'protein':
            BASE /
            'Xenopeltis_unicolor' /
            'Xenopeltis_unicolor.protein.faa'
    },

    'Candoia_aspera': {
        'annotation':
            BASE /
            'Candoia_aspera' /
            'Candoia_aspera.annotation.gff3',

        'protein':
            BASE /
            'Candoia_aspera' /
            'Candoia_aspera.protein.faa'
    },

    'Elaphe_schrenckii': {
        'annotation':
            BASE /
            'Elaphe_schrenckii' /
            'Elaphe_schrenckii.annotation.gff3',

        'protein':
            BASE /
            'Elaphe_schrenckii' /
            'Elaphe_schrenckii.protein.faa'
    },

    'Naja_naja': {
        'annotation':
            BASE /
            'Naja_naja' /
            'Naja_naja.annotation.gff3',

        'protein':
            BASE /
            'Naja_naja' /
            'Naja_naja.protein.faa'
    },

    'Cerastes_gasperettii': {
        'annotation':
            BASE /
            'Cerastes_gasperettii' /
            'Cerastes_gasperettii.annotation.gff3',

        'protein':
            BASE /
            'Cerastes_gasperettii' /
            'Cerastes_gasperettii.protein.faa'
    },

    'Vipera_berus': {
        'annotation':
            BASE /
            'Vipera_berus' /
            'Vipera_berus.annotation.gff3',

        'protein':
            BASE /
            'Vipera_berus' /
            'Vipera_berus.protein.faa'
    },

    'Bothrops_insularis': {
        'annotation':
            BASE /
            'Bothrops_insularis' /
            'Bothrops_insularis.annotation.gtf',

        'protein':
            BASE /
            'Bothrops_insularis' /
            'Bothrops_insularis.protein.faa'
    },

    'Crotalus_adamanteus': {
        'annotation':
            BASE /
            'Crotalus_adamanteus' /
            'Crotalus_adamanteus.annotation.gff3',

        'protein':
            BASE /
            'Crotalus_adamanteus' /
            'Crotalus_adamanteus.protein.faa'
    },

    'Gloydius_shedaoensis': {
        'annotation':
            BASE /
            'Gloydius_shedaoensis' /
            'Gloydius_shedaoensis.annotation.gff3',

        'protein':
            BASE /
            'Gloydius_shedaoensis' /
            'Gloydius_shedaoensis.protein.faa'
    },

    'Gloydius_ussuriensis': {
        'annotation':
            GUSS_DIR /
            'Gloydius_ussuriensis_AMNH_21010.gff3',

        'protein':
            GUSS_DIR /
            'Gloydius_ussuriensis_AMNH_21010.proteins.fa'
    }
}


# ============================================================
# helper functions
# ============================================================

def parse_gff3_attributes(text):

    out = {}

    for item in text.strip().split(';'):

        if not item:
            continue

        if '=' in item:

            key, value = item.split(
                '=',
                1
            )

            out[key] = value

    return out


def parse_gtf_attributes(text):

    return {
        match.group(1):
            match.group(2)

        for match in re.finditer(
            r'(\S+)\s+"([^"]*)"',
            text
        )
    }


# ------------------------------------------------------------
# chromosome ordering
# ------------------------------------------------------------

def chromosome_key(chrom):

    match = re.fullmatch(
        r'chr([0-9]+)',
        chrom
    )

    if match:

        return (
            0,
            int(
                match.group(1)
            )
        )

    if chrom == 'chrZ':

        return (
            1,
            0
        )

    if chrom == 'chrW':

        return (
            2,
            0
        )

    return (
        3,
        chrom
    )


# ------------------------------------------------------------
# FASTA reader
#
# Also detects duplicate FASTA primary IDs instead of
# silently overwriting them.
# ------------------------------------------------------------

def read_fasta(path):

    records = {}
    headers = {}

    current_id = None
    current_header = None
    sequence = []


    def save_record():

        if current_id is None:
            return

        if current_id in records:

            raise RuntimeError(
                f'Duplicate FASTA ID in {path}: '
                f'{current_id}'
            )

        records[
            current_id
        ] = ''.join(
            sequence
        )

        headers[
            current_id
        ] = current_header


    with path.open() as handle:

        for line in handle:

            line = line.rstrip(
                '\n'
            )

            if line.startswith('>'):

                save_record()

                current_header = line[
                    1:
                ].strip()

                current_id = current_header.split()[0]

                sequence = []

            else:

                sequence.append(
                    line.strip()
                )


    save_record()

    return records, headers


# ------------------------------------------------------------
# sanitize protein sequences
#
# Periods caused DIAMOND makedb failure previously.
#
# Only "." is replaced.
# U is deliberately retained.
# ------------------------------------------------------------

def clean_protein_sequence(sequence):

    sequence = ''.join(
        sequence.split()
    )

    sequence = sequence.upper()

    sequence = sequence.replace(
        '.',
        'X'
    )

    return sequence


# ============================================================
# read chromosome manifest
# ============================================================

if not MANIFEST.exists():

    raise FileNotFoundError(
        f'Missing chromosome manifest: {MANIFEST}'
    )


manifest_rows = []

with MANIFEST.open() as handle:

    reader = csv.DictReader(
        handle,
        delimiter='\t'
    )

    for row in reader:

        manifest_rows.append(
            row
        )


# ------------------------------------------------------------
# build chromosome whitelist
#
# keep_map:
#
#   original annotation sequence ID
#          ->
#   standardized chromosome label
#
# Explicit W chromosomes are excluded here.
#
# Because the manifest contains chromosome pseudomolecules
# only, all unplaced/unlocalized scaffolds are automatically
# absent from keep_map.
# ------------------------------------------------------------

keep_maps = defaultdict(dict)

explicit_w = defaultdict(list)

manifest_kept_chromosomes = defaultdict(set)


for row in manifest_rows:

    species = row[
        'species'
    ]

    if species not in species_order:
        continue


    if row[
        'annotation_match'
    ] != 'PASS':

        raise RuntimeError(
            f'Non-PASS chromosome mapping found for '
            f'{species}: {row}'
        )


    seqid = row[
        'annotation_seqid'
    ]

    standard_chr = row[
        'standard_chr'
    ]

    sex_class = row[
        'sex_class'
    ]


    if not seqid:

        raise RuntimeError(
            f'Blank annotation sequence ID for '
            f'{species}: {row}'
        )


    # --------------------------------------------------------
    # explicit W -> exclude
    # --------------------------------------------------------

    if sex_class == 'W':

        explicit_w[
            species
        ].append(
            seqid
        )

        continue


    # --------------------------------------------------------
    # safety against duplicated mappings
    # --------------------------------------------------------

    if seqid in keep_maps[
        species
    ]:

        raise RuntimeError(
            f'Duplicate manifest sequence ID for '
            f'{species}: {seqid}'
        )


    if standard_chr in manifest_kept_chromosomes[
        species
    ]:

        raise RuntimeError(
            f'Duplicate standardized chromosome for '
            f'{species}: {standard_chr}'
        )


    keep_maps[
        species
    ][
        seqid
    ] = standard_chr


    manifest_kept_chromosomes[
        species
    ].add(
        standard_chr
    )


# ------------------------------------------------------------
# manifest QC
# ------------------------------------------------------------

print()
print('=' * 72)
print('CHROMOSOME MANIFEST FILTER')
print('=' * 72)


for species in species_order:

    retained_labels = sorted(
        manifest_kept_chromosomes[
            species
        ],
        key=chromosome_key
    )


    print()
    print(species)

    print(
        '  chromosomes retained :',
        len(
            retained_labels
        )
    )

    print(
        '  labels               :',
        ', '.join(
            retained_labels
        )
    )

    print(
        '  explicit W removed   :',
        (
            ', '.join(
                explicit_w[
                    species
                ]
            )
            if explicit_w[
                species
            ]
            else 'none'
        )
    )


    if not retained_labels:

        raise RuntimeError(
            f'No chromosomes retained for {species}'
        )


# ============================================================
# parse source annotation
# ============================================================

def parse_annotation(path):

    is_gtf = (
        path.suffix.lower() == '.gtf'
    )


    gene_coords = {}

    transcript_coords = {}

    transcript_to_gene = {}

    # annotation protein_id -> transcript ID(s)

    protein_to_transcript = defaultdict(
        set
    )


    with path.open() as handle:

        for line in handle:

            if (
                not line.strip()
                or
                line.startswith('#')
            ):

                continue


            fields = line.rstrip(
                '\n'
            ).split('\t')


            if len(fields) != 9:
                continue


            chrom = fields[0]

            feature = fields[2]

            start = int(
                fields[3]
            )

            end = int(
                fields[4]
            )

            strand = fields[6]


            if is_gtf:

                attrs = parse_gtf_attributes(
                    fields[8]
                )

            else:

                attrs = parse_gff3_attributes(
                    fields[8]
                )


            # =================================================
            # gene
            # =================================================

            if feature == 'gene':

                if is_gtf:

                    gene_id = attrs.get(
                        'gene_id'
                    )

                else:

                    gene_id = attrs.get(
                        'ID'
                    )


                if gene_id:

                    gene_coords[
                        gene_id
                    ] = (
                        chrom,
                        start,
                        end,
                        strand
                    )


            # =================================================
            # transcript / mRNA
            # =================================================

            if feature in {
                'mRNA',
                'transcript'
            }:

                if is_gtf:

                    transcript_id = attrs.get(
                        'transcript_id'
                    )

                    gene_id = attrs.get(
                        'gene_id'
                    )

                else:

                    transcript_id = attrs.get(
                        'ID'
                    )

                    parent = attrs.get(
                        'Parent'
                    )

                    gene_id = (
                        parent.split(',')[0]
                        if parent
                        else None
                    )


                if transcript_id:

                    transcript_coords[
                        transcript_id
                    ] = (
                        chrom,
                        start,
                        end,
                        strand
                    )


                    if gene_id:

                        transcript_to_gene[
                            transcript_id
                        ] = gene_id


            # =================================================
            # CDS
            #
            # Needed for protein_id -> transcript mappings.
            # =================================================

            if feature == 'CDS':

                if is_gtf:

                    transcript_id = attrs.get(
                        'transcript_id'
                    )

                    gene_id = attrs.get(
                        'gene_id'
                    )

                    protein_id = attrs.get(
                        'protein_id'
                    )


                    if (
                        transcript_id
                        and
                        gene_id
                    ):

                        transcript_to_gene[
                            transcript_id
                        ] = gene_id


                else:

                    parent = attrs.get(
                        'Parent'
                    )

                    transcript_id = (
                        parent.split(',')[0]
                        if parent
                        else None
                    )

                    protein_id = attrs.get(
                        'protein_id'
                    )


                if (
                    protein_id
                    and
                    transcript_id
                ):

                    protein_to_transcript[
                        protein_id
                    ].add(
                        transcript_id
                    )


    return (
        gene_coords,
        transcript_coords,
        transcript_to_gene,
        protein_to_transcript
    )


# ============================================================
# process one species
# ============================================================

def process_species(
    species,
    annotation,
    protein_fasta
):

    print()
    print('=' * 72)
    print(species)
    print('=' * 72)


    keep_map = keep_maps[
        species
    ]

    w_seqids = set(
        explicit_w[
            species
        ]
    )


    # --------------------------------------------------------
    # source files
    # --------------------------------------------------------

    if not annotation.exists():

        raise FileNotFoundError(
            annotation
        )


    if not protein_fasta.exists():

        raise FileNotFoundError(
            protein_fasta
        )


    # --------------------------------------------------------
    # proteins
    # --------------------------------------------------------

    proteins, protein_headers = read_fasta(
        protein_fasta
    )


    # --------------------------------------------------------
    # annotation
    # --------------------------------------------------------

    (
        gene_coords,
        transcript_coords,
        transcript_to_gene,
        protein_to_transcript
    ) = parse_annotation(
        annotation
    )


    # --------------------------------------------------------
    # locus -> protein candidates
    # --------------------------------------------------------

    candidates = defaultdict(
        list
    )


    unmapped = []

    ambiguous = []


    direct_transcript = 0

    via_protein_id = 0

    transcript_as_locus = 0


    # ========================================================
    # map proteins to biological loci
    # ========================================================

    for protein_id, raw_sequence in proteins.items():

        sequence = clean_protein_sequence(
            raw_sequence
        )


        transcript_id = None

        method = None


        # ----------------------------------------------------
        # Pattern 1:
        # FASTA ID == transcript/mRNA ID
        # ----------------------------------------------------

        if (
            protein_id in transcript_coords
            or
            protein_id in transcript_to_gene
        ):

            transcript_id = protein_id

            method = (
                'fasta_id_is_transcript'
            )

            direct_transcript += 1


        # ----------------------------------------------------
        # Pattern 2:
        # FASTA ID == annotation protein_id
        # ----------------------------------------------------

        elif protein_id in protein_to_transcript:

            transcript_ids = sorted(
                protein_to_transcript[
                    protein_id
                ]
            )


            if len(
                transcript_ids
            ) != 1:

                ambiguous.append(
                    (
                        protein_id,
                        ','.join(
                            transcript_ids
                        )
                    )
                )

                continue


            transcript_id = transcript_ids[
                0
            ]

            method = (
                'fasta_id_is_protein_id'
            )

            via_protein_id += 1


        else:

            unmapped.append(
                protein_id
            )

            continue


        # ----------------------------------------------------
        # resolve biological gene/locus
        # ----------------------------------------------------

        gene_id = transcript_to_gene.get(
            transcript_id
        )


        if (
            gene_id
            and
            gene_id in gene_coords
        ):

            locus_id = gene_id

            coords = gene_coords[
                gene_id
            ]


        elif transcript_id in transcript_coords:

            # fallback for annotations without explicit
            # gene features
            #
            # Crotalus viridis would use this extensively,
            # but C. viridis is not part of this run.

            locus_id = transcript_id

            coords = transcript_coords[
                transcript_id
            ]

            transcript_as_locus += 1


        else:

            unmapped.append(
                protein_id
            )

            continue


        candidates[
            locus_id
        ].append(
            (
                protein_id,
                transcript_id,
                sequence,
                method,
                coords
            )
        )


    # ========================================================
    # select longest protein per locus
    # ========================================================

    representatives = {}


    for locus_id, values in candidates.items():

        values = sorted(
            values,
            key=lambda x: (
                -len(
                    x[2]
                ),
                x[0],
                x[1]
            )
        )


        representatives[
            locus_id
        ] = values[
            0
        ]


    # ========================================================
    # chromosome filtering
    #
    # Anything not in keep_map is excluded:
    #
    #   - unplaced scaffolds
    #   - unlocalized scaffolds
    #   - contigs
    #   - explicitly annotated W chromosome
    #
    # Numeric potential W chromosomes remain for now.
    # ========================================================

    retained = []


    removed_explicit_w = 0

    removed_nonchromosome = 0


    for locus_id, value in representatives.items():

        (
            protein_id,
            transcript_id,
            sequence,
            method,
            coords
        ) = value


        source_chrom = coords[
            0
        ]


        if source_chrom in w_seqids:

            removed_explicit_w += 1

            continue


        if source_chrom not in keep_map:

            removed_nonchromosome += 1

            continue


        standard_chrom = keep_map[
            source_chrom
        ]


        retained.append(
            {
                'locus_id':
                    locus_id,

                'protein_id':
                    protein_id,

                'transcript_id':
                    transcript_id,

                'sequence':
                    sequence,

                'method':
                    method,

                'source_chrom':
                    source_chrom,

                'standard_chrom':
                    standard_chrom,

                'start':
                    coords[
                        1
                    ],

                'end':
                    coords[
                        2
                    ],

                'strand':
                    coords[
                        3
                    ],

                'n_isoforms':
                    len(
                        candidates[
                            locus_id
                        ]
                    )
            }
        )


    # --------------------------------------------------------
    # sort retained genes genomically
    # --------------------------------------------------------

    retained.sort(
        key=lambda rec: (
            chromosome_key(
                rec[
                    'standard_chrom'
                ]
            ),
            rec[
                'start'
            ],
            rec[
                'end'
            ],
            rec[
                'locus_id'
            ]
        )
    )


    # ========================================================
    # output paths
    # ========================================================

    bed_path = (
        BED_DIR /
        f'{species}.bed'
    )


    peptide_path = (
        PEP_DIR /
        f'{species}.fa'
    )


    mapping_path = (
        MAP_DIR /
        f'{species}.longest_isoform.tsv'
    )


    # ========================================================
    # BED
    #
    # GFF/GTF:
    #     1-based closed
    #
    # BED:
    #     0-based start
    #
    # Column 1:
    #     standardized chromosome
    #
    # Column 4:
    #     biological locus/gene ID
    # ========================================================

    with bed_path.open(
        'w'
    ) as out:

        for rec in retained:

            out.write(
                f"{rec['standard_chrom']}\t"
                f"{rec['start'] - 1}\t"
                f"{rec['end']}\t"
                f"{rec['locus_id']}\n"
            )


    # ========================================================
    # representative peptide FASTA
    #
    # Protein ID is rewritten to locus ID so:
    #
    # BED column 4 == peptide FASTA ID
    # ========================================================

    with peptide_path.open(
        'w'
    ) as out:

        for rec in retained:

            out.write(
                f">{rec['locus_id']}\n"
            )


            sequence = rec[
                'sequence'
            ]


            for i in range(
                0,
                len(sequence),
                60
            ):

                out.write(
                    sequence[
                        i:
                        i + 60
                    ]
                    +
                    '\n'
                )


    # ========================================================
    # mapping table
    # ========================================================

    with mapping_path.open(
        'w',
        newline=''
    ) as out:

        fields = [
            'locus_id',
            'representative_protein_id',
            'representative_transcript_id',
            'protein_length_aa',
            'source_chromosome',
            'standard_chromosome',
            'start',
            'end',
            'strand',
            'n_protein_isoforms',
            'mapping_method'
        ]


        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter='\t'
        )


        writer.writeheader()


        for rec in retained:

            writer.writerow(
                {
                    'locus_id':
                        rec[
                            'locus_id'
                        ],

                    'representative_protein_id':
                        rec[
                            'protein_id'
                        ],

                    'representative_transcript_id':
                        rec[
                            'transcript_id'
                        ],

                    'protein_length_aa':
                        len(
                            rec[
                                'sequence'
                            ]
                        ),

                    'source_chromosome':
                        rec[
                            'source_chrom'
                        ],

                    'standard_chromosome':
                        rec[
                            'standard_chrom'
                        ],

                    'start':
                        rec[
                            'start'
                        ],

                    'end':
                        rec[
                            'end'
                        ],

                    'strand':
                        rec[
                            'strand'
                        ],

                    'n_protein_isoforms':
                        rec[
                            'n_isoforms'
                        ],

                    'mapping_method':
                        rec[
                            'method'
                        ]
                }
            )


    # ========================================================
    # QC
    # ========================================================

    locus_ids = [
        rec[
            'locus_id'
        ]
        for rec in retained
    ]


    if len(
        locus_ids
    ) != len(
        set(
            locus_ids
        )
    ):

        raise RuntimeError(
            f'Duplicate retained locus IDs: '
            f'{species}'
        )


    # --------------------------------------------------------
    # no explicit chrW may survive
    # --------------------------------------------------------

    if any(
        rec[
            'standard_chrom'
        ] == 'chrW'

        for rec in retained
    ):

        raise RuntimeError(
            f'Explicit W survived filtering: '
            f'{species}'
        )


    # --------------------------------------------------------
    # chromosomes actually represented in outputs
    # --------------------------------------------------------

    retained_chromosomes = sorted(
        {
            rec[
                'standard_chrom'
            ]
            for rec in retained
        },
        key=chromosome_key
    )


    # --------------------------------------------------------
    # BED IDs
    # --------------------------------------------------------

    bed_ids = []


    with bed_path.open() as handle:

        for line in handle:

            fields = line.rstrip(
                '\n'
            ).split('\t')


            if len(fields) != 4:

                raise RuntimeError(
                    f'Invalid BED line for '
                    f'{species}: {line}'
                )


            bed_ids.append(
                fields[
                    3
                ]
            )


    # --------------------------------------------------------
    # peptide IDs
    # --------------------------------------------------------

    output_proteins, _ = read_fasta(
        peptide_path
    )


    peptide_ids = list(
        output_proteins.keys()
    )


    if len(
        bed_ids
    ) != len(
        peptide_ids
    ):

        raise RuntimeError(
            f'BED/peptide count mismatch: '
            f'{species}'
        )


    if set(
        bed_ids
    ) != set(
        peptide_ids
    ):

        raise RuntimeError(
            f'BED/peptide ID mismatch: '
            f'{species}'
        )


    # --------------------------------------------------------
    # invalid peptide characters
    #
    # Letters A-Z and stop "*" are permitted here.
    #
    # "." should already have become X.
    # --------------------------------------------------------

    invalid_characters = Counter()


    allowed = set(
        'ABCDEFGHIJKLMNOPQRSTUVWXYZ*'
    )


    for sequence in output_proteins.values():

        for aa in sequence:

            if aa not in allowed:

                invalid_characters[
                    aa
                ] += 1


    if invalid_characters:

        raise RuntimeError(
            f'Invalid peptide characters for '
            f'{species}: '
            f'{dict(invalid_characters)}'
        )


    # --------------------------------------------------------
    # summary metrics
    # --------------------------------------------------------

    genes_with_protein = len(
        set(
            candidates
        )
        &
        set(
            gene_coords
        )
    )


    explicit_genes_without_protein = len(
        set(
            gene_coords
        )
        -
        set(
            candidates
        )
    )


    print(
        f'Input proteins:                     '
        f'{len(proteins):,}'
    )

    print(
        f'Annotation gene features:           '
        f'{len(gene_coords):,}'
    )

    print(
        f'Annotation transcript features:     '
        f'{len(transcript_coords):,}'
    )

    print(
        f'Mapped directly by transcript ID:   '
        f'{direct_transcript:,}'
    )

    print(
        f'Mapped through CDS protein_id:      '
        f'{via_protein_id:,}'
    )

    print(
        f'Protein-mapped loci:                '
        f'{len(candidates):,}'
    )

    print(
        f'Representative loci before filter:  '
        f'{len(representatives):,}'
    )

    print(
        f'Chromosome loci retained:           '
        f'{len(retained):,}'
    )

    print(
        f'Loci removed on explicit W:         '
        f'{removed_explicit_w:,}'
    )

    print(
        f'Loci removed off chromosomes:       '
        f'{removed_nonchromosome:,}'
    )

    print(
        f'Explicit genes represented:         '
        f'{genes_with_protein:,}'
    )

    print(
        f'Explicit genes without protein:     '
        f'{explicit_genes_without_protein:,}'
    )

    print(
        f'Transcript-as-locus assignments:    '
        f'{transcript_as_locus:,}'
    )

    print(
        f'Unmapped proteins:                  '
        f'{len(unmapped):,}'
    )

    print(
        f'Ambiguous protein mappings:         '
        f'{len(ambiguous):,}'
    )

    print(
        f'Chromosomes represented:            '
        f'{len(retained_chromosomes)}'
    )

    print(
        'Chromosome labels:                 '
        +
        ', '.join(
            retained_chromosomes
        )
    )

    print(
        'Explicit W source sequence(s):     '
        +
        (
            ', '.join(
                sorted(
                    w_seqids
                )
            )
            if w_seqids
            else 'none'
        )
    )

    print(
        'BED/peptide ID QC:                 PASS'
    )

    print(
        'Protein character QC:              PASS'
    )


    if unmapped:

        print()
        print(
            'First unmapped proteins:'
        )

        for protein_id in unmapped[
            :10
        ]:

            print(
                ' ',
                protein_id
            )


    if ambiguous:

        print()
        print(
            'First ambiguous mappings:'
        )

        for protein_id, tids in ambiguous[
            :10
        ]:

            print(
                ' ',
                protein_id,
                '->',
                tids
            )


    return {
        'species':
            species,

        'input_proteins':
            len(
                proteins
            ),

        'gene_features':
            len(
                gene_coords
            ),

        'transcript_features':
            len(
                transcript_coords
            ),

        'protein_mapped_loci':
            len(
                candidates
            ),

        'representative_loci_before_filter':
            len(
                representatives
            ),

        'retained_loci':
            len(
                retained
            ),

        'removed_explicit_W_loci':
            removed_explicit_w,

        'removed_nonchromosome_loci':
            removed_nonchromosome,

        'unmapped_proteins':
            len(
                unmapped
            ),

        'ambiguous_proteins':
            len(
                ambiguous
            ),

        'retained_chromosomes':
            len(
                retained_chromosomes
            ),

        'chromosome_labels':
            ','.join(
                retained_chromosomes
            ),

        'explicit_W_source_sequences':
            ','.join(
                sorted(
                    w_seqids
                )
            )
    }


# ============================================================
# source-file check
# ============================================================

print()
print('=' * 72)
print('SOURCE FILE CHECK')
print('=' * 72)


for species in species_order:

    annotation = sources[
        species
    ][
        'annotation'
    ]

    protein = sources[
        species
    ][
        'protein'
    ]


    if not annotation.exists():

        raise FileNotFoundError(
            f'Missing annotation for '
            f'{species}: {annotation}'
        )


    if not protein.exists():

        raise FileNotFoundError(
            f'Missing protein FASTA for '
            f'{species}: {protein}'
        )


    print(
        f'{species:25s} PASS'
    )


# ============================================================
# run all 11 species
# ============================================================

summary = []


for species in species_order:

    result = process_species(
        species=species,

        annotation=sources[
            species
        ][
            'annotation'
        ],

        protein_fasta=sources[
            species
        ][
            'protein'
        ]
    )


    summary.append(
        result
    )


# ============================================================
# write global summary
# ============================================================

summary_path = (
    QC_DIR /
    'chromosome_input_summary.tsv'
)


summary_fields = [
    'species',
    'input_proteins',
    'gene_features',
    'transcript_features',
    'protein_mapped_loci',
    'representative_loci_before_filter',
    'retained_loci',
    'removed_explicit_W_loci',
    'removed_nonchromosome_loci',
    'unmapped_proteins',
    'ambiguous_proteins',
    'retained_chromosomes',
    'chromosome_labels',
    'explicit_W_source_sequences'
]


with summary_path.open(
    'w',
    newline=''
) as out:

    writer = csv.DictWriter(
        out,
        fieldnames=summary_fields,
        delimiter='\t'
    )


    writer.writeheader()

    writer.writerows(
        summary
    )


# ============================================================
# final summary
# ============================================================

print()
print()
print('=' * 72)
print('FINAL SUMMARY')
print('=' * 72)


print(
    'species\t'
    'retained_loci\t'
    'chromosomes\t'
    'W_loci_removed\t'
    'nonchrom_loci_removed\t'
    'unmapped\t'
    'ambiguous'
)


for row in summary:

    print(
        f"{row['species']}\t"
        f"{row['retained_loci']}\t"
        f"{row['retained_chromosomes']}\t"
        f"{row['removed_explicit_W_loci']}\t"
        f"{row['removed_nonchromosome_loci']}\t"
        f"{row['unmapped_proteins']}\t"
        f"{row['ambiguous_proteins']}"
    )


print()
print(
    'Summary written to:'
)

print(
    summary_path
)

print()
print(
    'GENESPACE chromosome-only inputs generated successfully.'
)

print(
    'Unplaced/unlocalized sequences were excluded.'
)

print(
    'Explicitly annotated W chromosomes were excluded.'
)

print(
    'Z chromosomes were retained.'
)

print(
    'Numerically labelled potential sex chromosomes were retained.'
)