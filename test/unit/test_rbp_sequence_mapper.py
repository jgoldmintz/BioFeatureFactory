"""RBP sequence/MSA source precedence regressions."""

from biofeaturefactory.alphafold3.bin.rbp_sequence_mapper import RBPSequenceMapper


def test_msa_is_preferred_while_fasta_supplies_missing_and_corrupt_entries(tmp_path):
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text(
        "Entry\tGene Names\tProtein names\tLength\n"
        "P11111\tWITHMSA\tProtein with MSA\t4\n"
        "P22222\tFASTAONLY\tProtein without MSA\t5\n"
        "P33333\tBROKENMSA\tProtein with corrupt MSA\t6\n"
    )
    sequences = tmp_path / "sequences.fasta"
    sequences.write_text(
        ">sp|P11111|WITHMSA_HUMAN\nFFFF\n"
        ">sp|P22222|FASTAONLY_HUMAN\nMPEPT\n"
        ">sp|P33333|BROKENMSA_HUMAN\nMPEPTI\n"
    )
    msa_dir = tmp_path / "msa"
    msa_dir.mkdir()
    (msa_dir / "AF-P11111-F1-msa_v6.a3m").write_text(
        ">query\nMSSA\n>hit\nMSSV\n"
    )
    (msa_dir / "AF-P33333-F1-msa_v6.a3m").write_text("")

    mapper = RBPSequenceMapper(
        mapping_file=str(mapping),
        sequence_fasta=str(sequences),
        msa_dir=str(msa_dir),
    )

    with_msa = mapper.get_rbp_data("WITHMSA")
    fasta_only = mapper.get_rbp_data("FASTAONLY")
    broken_msa = mapper.get_rbp_data("BROKENMSA")

    assert with_msa.sequence == "MSSA"
    assert with_msa.msa_content is not None
    assert fasta_only.sequence == "MPEPT"
    assert fasta_only.msa_content is None
    assert broken_msa.sequence == "MPEPTI"
    assert broken_msa.msa_content is None
