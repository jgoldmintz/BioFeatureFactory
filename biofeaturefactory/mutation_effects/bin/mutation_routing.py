"""Select protein and codon work from validated, per-gene mutation inputs."""

from Bio.Data.CodonTable import TranslationError

from biofeaturefactory.lib.utility import (
    load_validation_failures,
    parse_variant,
    protein_consequence,
    read_fasta,
    should_skip_mutation,
    splice_seq,
    trim_muts,
)


CODON_CLASSES = frozenset({"SYNONYMOUS", "STOP_GAIN", "STOP_LOSS"})


def classify_gene(fasta, mutations, gene, validation_log=None):
    """Return distinct consequence classes, retaining uncertainty as UNKNOWN."""
    tokens = trim_muts(mutations, log=validation_log, gene_name=gene)
    failures = load_validation_failures(validation_log)
    tokens = [token for token in tokens if not should_skip_mutation(gene, token, failures)]
    if not tokens:
        return []

    sequences = read_fasta(fasta)
    orf_sequence = next(
        (sequence for header, sequence in sequences.items() if header.upper() == "ORF"),
        next(iter(sequences.values()), ""),
    ).upper().replace("U", "T")
    classes = set()
    for token in tokens:
        variant = parse_variant(token.upper().replace("U", "T"), is_nt=True)
        if variant is None:
            classes.add("UNKNOWN")
            continue

        first_codon_start = variant.pos0 // 3 * 3
        last_codon_end = (variant.pos0 + len(variant.ref) + 2) // 3 * 3
        if last_codon_end > len(orf_sequence) or any(
            nucleotide not in "ACGT" for nucleotide in orf_sequence[first_codon_start:last_codon_end]
        ):
            classes.add("UNKNOWN")
            continue

        try:
            mutant_sequence = splice_seq(orf_sequence, variant.pos0, variant.ref, variant.alt)
            consequence = protein_consequence(variant, orf_sequence, mut_orf_seq=mutant_sequence)
        except (ValueError, TranslationError):
            consequence = None
        if consequence is None or "X" in consequence["wt_aa"] + consequence["mut_aa"]:
            classes.add("UNKNOWN")
            continue
        label = consequence["aa_consequence"]
        classes.add({"snv": "MISSENSE", "stop_gained": "STOP_GAIN", "stop_lost": "STOP_LOSS"}.get(label, label.upper()))
    return sorted(classes)


def choose_route(classes, protein_explicit=False, codon_explicit=False):
    """Return selection and warnings without printing or launching any work."""
    classes = sorted({label.upper() for label in classes})
    if protein_explicit and codon_explicit:
        mode = "both"
    elif protein_explicit:
        mode = "protein"
    elif codon_explicit:
        mode = "codon"
    else:
        mode = "auto"

    codon_classes = set(classes).intersection(CODON_CLASSES)
    protein_classes = set(classes).difference(CODON_CLASSES)
    if mode == "auto":
        protein = bool(protein_classes)
        codon = bool(codon_classes or "UNKNOWN" in classes)
    else:
        protein = bool(classes) and protein_explicit
        codon = bool(classes) and codon_explicit

    warnings = []
    if mode == "protein" and codon_classes:
        warnings.append(
            "Protein-only routing includes synonymous or stop variants: this will not produce "
            "biologically accurate results because protein-level scores do not capture their codon effects."
        )
    if mode == "codon" and protein_classes:
        warnings.append(
            "Codon scoring of amino-acid-changing or unclassified variants has higher memory and "
            "computational costs and reports per-codon effects, which differ from per-amino-acid effects."
        )
    return {
        "protein": protein,
        "codon": codon,
        "classes": classes,
        "mode": mode,
        "score_missense_codon": bool(mode == "codon" and codon and protein_classes),
        "warnings": warnings,
    }
