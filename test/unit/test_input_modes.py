import argparse
import csv
from pathlib import Path
import sys

import pytest

from biofeaturefactory.lib.utility import (
    InputPathAction,
    derive_mutations_root,
    discover_mutation_files,
    mint_pkey,
    validate_input_mode,
)


@pytest.fixture
def inputs(tmp_path):
    root = tmp_path / 'out'
    gene = root / 'GENE'
    for folder in ('fastas', 'mappings/mutations', 'MSA', 'CodonMSA'):
        (gene / folder).mkdir(parents=True, exist_ok=True)
    fasta = gene / 'fastas/GENE.fasta'
    fasta.write_text('>ORF\nATGGCTTAA\n')
    mutations = gene / 'mappings/mutations/GENE_mutations.csv'
    mutations.write_text('mutant\nG4A\n')
    return root, fasta, mutations


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', action=InputPathAction,
                        extensions=('.fa', '.fasta'), type=Path)
    parser.add_argument('-m', '--mutations', action=InputPathAction,
                        extensions=('.csv', '.tsv', '.txt'))
    parser.add_argument('--mapping', action=InputPathAction, extensions=('.csv', '.tsv'))
    parser.add_argument('--reference')
    parser.add_argument('-o', '--output')
    return parser


def test_explicit_files_and_reference_do_not_change_mode(inputs, tmp_path):
    _, fasta, mutations = inputs
    parser = make_parser()
    args = parser.parse_args(['-o', str(tmp_path), '--reference', str(tmp_path),
                              '-m', str(mutations), '-f', str(fasta)])
    assert validate_input_mode(parser, args, ('mutations',)) == 'file'
    assert args.fasta == fasta
    assert args.mutations == str(mutations)
    assert not hasattr(args, '_bff_input_paths')


def test_parent_directory_keeps_automatic_mutations(inputs):
    root, _, mutations = inputs
    parser = make_parser()
    args = parser.parse_args(['-f', str(root)])
    assert validate_input_mode(parser, args, ('mutations',)) == 'directory'
    assert args.mutations is None
    derived = derive_mutations_root(args.mutations, args.fasta)
    assert derived == root
    assert discover_mutation_files(derived) == {'GENE': str(mutations)}


@pytest.mark.parametrize('file_first', [True, False])
def test_first_supplied_flag_controls_mixed_mode_error(inputs, file_first, capsys):
    root, fasta, _ = inputs
    parser = make_parser()
    file_option = ['-f', str(fasta)]
    directory_option = ['-m', str(root)]
    args = parser.parse_args(file_option + directory_option if file_first
                             else directory_option + file_option)
    with pytest.raises(SystemExit) as error:
        validate_input_mode(parser, args)
    assert error.value.code == 2
    expected = '-f selected file mode' if file_first else '-m selected directory mode'
    assert expected in capsys.readouterr().err


@pytest.mark.parametrize('relative', ['GENE', 'GENE/fastas', 'GENE/mappings',
                                      'GENE/mappings/mutations', 'GENE/MSA', 'GENE/CodonMSA'])
def test_nested_directories_are_rejected_with_parent_hint(inputs, relative, capsys):
    root, _, _ = inputs
    parser = make_parser()
    args = parser.parse_args(['-f', str(root / relative)])
    with pytest.raises(SystemExit):
        validate_input_mode(parser, args)
    message = capsys.readouterr().err
    assert 'parent <dir> containing <gene>/...' in message
    assert f'Provide {root} instead.' in message


def test_relative_subdirectory_hint_uses_real_parent(inputs, monkeypatch, capsys):
    root, _, _ = inputs
    monkeypatch.chdir(root / 'GENE')
    parser = make_parser()
    args = parser.parse_args(['-f', 'fastas'])
    with pytest.raises(SystemExit):
        validate_input_mode(parser, args)
    assert f'Provide {root} instead.' in capsys.readouterr().err


def test_file_mode_does_not_derive_mutations_from_fasta(inputs, capsys):
    _, fasta, _ = inputs
    parser = make_parser()
    args = parser.parse_args(['-f', str(fasta)])
    with pytest.raises(SystemExit) as error:
        validate_input_mode(parser, args, ('mutations',))
    assert error.value.code == 2
    assert '--mutations is required in file mode' in capsys.readouterr().err
    assert args.mutations is None


def test_alternative_file_companion(inputs):
    _, fasta, mutations = inputs
    parser = make_parser()
    args = parser.parse_args(['-f', str(fasta), '--mapping', str(mutations)])
    assert validate_input_mode(parser, args, (('mutations', 'mapping'),)) == 'file'


@pytest.mark.parametrize('suffix', ['.csv', '.tsv', '.txt', '.CSV'])
def test_accepted_companion_extensions(inputs, tmp_path, suffix):
    _, fasta, _ = inputs
    mutations = tmp_path / ('GENE' + suffix)
    mutations.write_text('mutant\nG4A\n')
    parser = make_parser()
    args = parser.parse_args(['-m', str(mutations), '-f', str(fasta)])
    assert validate_input_mode(parser, args, ('mutations',)) == 'file'


def test_unsupported_existing_file_is_rejected(inputs, tmp_path, capsys):
    unsupported = tmp_path / 'GENE.bin'
    unsupported.write_bytes(b'not a FASTA')
    parser = make_parser()
    args = parser.parse_args(['-f', str(unsupported)])
    with pytest.raises(SystemExit):
        validate_input_mode(parser, args)
    assert 'unsupported input file type' in capsys.readouterr().err


def test_missing_file_does_not_fall_back_to_directory(tmp_path, capsys):
    parser = make_parser()
    args = parser.parse_args(['-f', str(tmp_path / 'absent.fasta')])
    with pytest.raises(SystemExit):
        validate_input_mode(parser, args)
    assert 'input path does not exist' in capsys.readouterr().err


def test_directory_with_file_extension_is_not_a_file(tmp_path, capsys):
    path = tmp_path / 'GENE.fasta'
    path.mkdir()
    parser = make_parser()
    args = parser.parse_args(['-f', str(path)])
    with pytest.raises(SystemExit):
        validate_input_mode(parser, args)
    assert 'expected an input file, but found a directory' in capsys.readouterr().err


def test_empty_directory_is_not_a_gene_tree(tmp_path, capsys):
    parser = make_parser()
    args = parser.parse_args(['-f', str(tmp_path)])
    with pytest.raises(SystemExit):
        validate_input_mode(parser, args)
    assert 'parent <dir>' in capsys.readouterr().err


def test_same_mode_explicit_roots_are_accepted(inputs):
    root, _, _ = inputs
    parser = make_parser()
    args = parser.parse_args(['-m', str(root), '-f', str(root)])
    assert validate_input_mode(parser, args, ('mutations',)) == 'directory'


def test_positional_fallback_value_is_validated(inputs):
    _, fasta, mutations = inputs
    parser = make_parser()
    args = parser.parse_args(['-m', str(mutations)])
    args.fasta = fasta
    assert validate_input_mode(parser, args, ('mutations',)) == 'file'


def test_home_expansion_preserves_path_types(inputs, monkeypatch):
    root, fasta, mutations = inputs
    monkeypatch.setenv('HOME', str(root))
    parser = make_parser()
    args = parser.parse_args(['-f', '~/GENE/fastas/GENE.fasta',
                              '-m', '~/GENE/mappings/mutations/GENE_mutations.csv'])
    assert validate_input_mode(parser, args, ('mutations',)) == 'file'
    assert args.fasta == fasta
    assert args.mutations == str(mutations)


def test_compressed_and_opaque_file_types(tmp_path):
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', action=InputPathAction, extensions=('.vcf', '.vcf.gz'))
    parser.add_argument('--params', action=InputPathAction)
    vcf = tmp_path / 'GENE.vcf.gz'
    parameters = tmp_path / 'GENE_model_parameters'
    vcf.write_bytes(b'fixture')
    parameters.write_bytes(b'fixture')
    args = parser.parse_args(['--vcf=' + str(vcf), '--params', str(parameters)])
    assert validate_input_mode(parser, args) == 'file'


def test_help_does_not_validate_paths(capsys):
    parser = make_parser()
    with pytest.raises(SystemExit) as error:
        parser.parse_args(['-f', 'missing.fasta', '--help'])
    assert error.value.code == 0
    assert '--fasta' in capsys.readouterr().out


@pytest.mark.parametrize('folder', ['vcf', 'EVmutation', 'adabmDCA'])
@pytest.mark.parametrize('nested', ['', 'GENE', 'GENE/tool'])
def test_partial_gene_trees_require_parent_root(tmp_path, folder, nested, capsys):
    root = tmp_path / 'out'
    (root / 'GENE' / folder).mkdir(parents=True)
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', action=InputPathAction)
    relative = nested.replace('tool', folder)
    args = parser.parse_args(['--source', str(root / relative)])
    if not nested:
        assert validate_input_mode(parser, args) == 'directory'
    else:
        with pytest.raises(SystemExit) as error:
            validate_input_mode(parser, args)
        assert error.value.code == 2
        assert f'Provide {root} instead.' in capsys.readouterr().err


def test_codon_usage_file_and_directory_exports_match(inputs, tmp_path, monkeypatch):
    from biofeaturefactory.codon_usage import codon_usage_pipeline

    root, fasta, mutations = inputs
    exported = []
    for mode, arguments in (
        ('directory', ['-f', str(root)]),
        ('file', ['-f', str(fasta), '-m', str(mutations)]),
    ):
        output = tmp_path / mode
        monkeypatch.setattr(sys, 'argv', ['codon_usage_pipeline.py', *arguments, '-o', str(output)])
        codon_usage_pipeline.main()
        with (output / 'GENE/CodonUsage/GENE.codon_usage.tsv').open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        assert len(rows) == 1
        assert rows[0]['pkey'] == mint_pkey('GENE', 'G4A')
        exported.append(rows)
    assert exported[0] == exported[1]
