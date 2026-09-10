"""RNA and AF3 CLI input-mode validation before prediction or output writes."""

from pathlib import Path

import pytest

from biofeaturefactory.RNAfold import run_viennaRNA_pipeline as rnafold
from biofeaturefactory.miranda import miranda_ensemble as miranda
from biofeaturefactory.genesplicer import genesplicer_ensemble as genesplicer
from biofeaturefactory.alphafold3 import alphafold3_pipeline as af3, burst
from biofeaturefactory.lib.utility import validate_input_mode


TOOLS = (
    ('rnafold', rnafold, '-i', '-tm', 'transcript'),
    ('miranda', miranda, '-i', '-M', 'transcript'),
    ('genesplicer', genesplicer, '-i', '-m', 'gDNA'),
    ('af3', af3, '-f', '-cm', 'chromosome'),
    ('burst', burst, '-f', '-cm', 'chromosome'),
)


@pytest.fixture
def inputs(tmp_path):
    root = tmp_path / 'out'
    gene = root / 'GENE'
    fasta = gene / 'fastas' / 'GENE.fasta'
    fasta.parent.mkdir(parents=True)
    fasta.write_text('>ORF\nAAAAAAA\n>transcript\nAAAAAAA\n>genomic\nAAAAAAA\n')
    mappings = {}
    for kind, filename, value in (
        ('transcript', 'transcript_mapping_GENE.csv', 'transcript'),
        ('gDNA', 'genomic_mapping_GENE.csv', 'genomic'),
        ('chromosome', 'chr_mapping_GENE.csv', 'chromosome'),
    ):
        mapping = gene / 'mappings' / kind / filename
        mapping.parent.mkdir(parents=True)
        mapping.write_text(f'mutant,{value}\nA4G,A4G\n')
        mappings[kind] = mapping
    mutations = gene / 'mappings' / 'mutations' / 'GENE_mutations.csv'
    mutations.parent.mkdir(parents=True)
    mutations.write_text('mutant\nA4G\n')
    vcf = gene / 'vcf' / 'GENE.vcf'
    vcf.parent.mkdir()
    vcf.write_text('#CHROM\tPOS\tID\tREF\tALT\n5\t4\t.\tA\tG\n')
    (gene / 'CodonMSA').mkdir()
    (gene / 'CodonMSA' / 'GENE.codon.msa.fasta').write_text('>focus\nAAAAAAA\n')
    models = tmp_path / 'models'
    models.mkdir()
    reference = tmp_path / 'reference.fasta'
    reference.write_text('>reference\nMKK\n')
    return {
        'root': root, 'gene': gene, 'fasta': fasta, 'mappings': mappings,
        'mutations': mutations, 'vcf': vcf, 'models': models, 'reference': reference,
        'output': tmp_path / 'unwritten',
    }


class ValidationBoundary(BaseException):
    pass


def invoke_validation(tool, inputs, arguments, monkeypatch):
    name, module, _, _, _ = tool
    observed = {}

    def boundary(parser, namespace, required_file_inputs=()):
        observed['mode'] = validate_input_mode(parser, namespace, required_file_inputs)
        observed['args'] = namespace
        raise ValidationBoundary()

    command = [module.__file__]
    if module is burst:
        command.append('ingest')
    command.extend(['-o', str(inputs['output'])])
    if module is miranda:
        command.extend(['-d', str(inputs['reference'])])
    if name in ('af3', 'burst'):
        command.extend(['--postar-db', str(inputs['reference']),
                        '--rbp-mapping', str(inputs['reference']),
                        '-rs', str(inputs['reference'])])
    if module is af3:
        command.extend(['-mdi', str(inputs['models'])])
    command.extend(map(str, arguments))
    monkeypatch.setattr('sys.argv', command)
    monkeypatch.setattr(module, 'validate_input_mode', boundary)
    try:
        module.main()
    except ValidationBoundary:
        pass
    assert not inputs['output'].exists()
    return observed


@pytest.mark.parametrize('tool', TOOLS, ids=[tool[0] for tool in TOOLS])
@pytest.mark.parametrize('location', ['gene', 'fastas', 'CodonMSA'])
def test_nested_directories_are_rejected_with_parent_root_hint(
        tool, inputs, monkeypatch, capsys, location):
    source = inputs['gene'] if location == 'gene' else inputs['gene'] / location
    with pytest.raises(SystemExit) as error:
        invoke_validation(tool, inputs, [tool[2], source], monkeypatch)
    assert error.value.code == 2
    assert str(inputs['root']) in capsys.readouterr().err
    assert not inputs['output'].exists()


@pytest.mark.parametrize('tool', TOOLS, ids=[tool[0] for tool in TOOLS])
def test_file_mode_requires_explicit_companion(tool, inputs, monkeypatch):
    with pytest.raises(SystemExit) as error:
        invoke_validation(tool, inputs, [tool[2], inputs['fasta']], monkeypatch)
    assert error.value.code == 2
    assert not inputs['output'].exists()


@pytest.mark.parametrize('tool', TOOLS, ids=[tool[0] for tool in TOOLS])
@pytest.mark.parametrize('mapping_first', [False, True])
def test_mixed_modes_are_rejected_in_both_argument_orders(
        tool, inputs, monkeypatch, mapping_first):
    fasta_pair = [tool[2], inputs['fasta']]
    mapping_pair = [tool[3], inputs['root']]
    arguments = mapping_pair + fasta_pair if mapping_first else fasta_pair + mapping_pair
    with pytest.raises(SystemExit) as error:
        invoke_validation(tool, inputs, arguments, monkeypatch)
    assert error.value.code == 2
    assert not inputs['output'].exists()


@pytest.mark.parametrize('tool', TOOLS, ids=[tool[0] for tool in TOOLS])
@pytest.mark.parametrize('mapping_first', [False, True])
def test_matching_file_inputs_reach_validation_boundary(
        tool, inputs, monkeypatch, mapping_first):
    fasta_pair = [tool[2], inputs['fasta']]
    mapping_pair = [tool[3], inputs['mappings'][tool[4]]]
    arguments = mapping_pair + fasta_pair if mapping_first else fasta_pair + mapping_pair
    observed = invoke_validation(tool, inputs, arguments, monkeypatch)
    assert observed['mode'] == 'file'


@pytest.mark.parametrize('tool', TOOLS, ids=[tool[0] for tool in TOOLS])
def test_parent_directory_inputs_reach_validation_boundary(tool, inputs, monkeypatch):
    observed = invoke_validation(
        tool, inputs, [tool[3], inputs['root'], tool[2], inputs['root']], monkeypatch)
    assert observed['mode'] == 'directory'


@pytest.mark.parametrize('tool', TOOLS, ids=[tool[0] for tool in TOOLS])
def test_wrong_companion_file_type_is_rejected(tool, inputs, monkeypatch):
    with pytest.raises(SystemExit) as error:
        invoke_validation(
            tool, inputs, [tool[2], inputs['fasta'], tool[3], inputs['fasta']], monkeypatch)
    assert error.value.code == 2


@pytest.mark.parametrize('module', [af3, burst])
def test_af3_mutations_file_can_supply_required_companion(module, inputs, monkeypatch):
    tool = next(candidate for candidate in TOOLS if candidate[1] is module)
    observed = invoke_validation(
        tool, inputs, ['-mu', inputs['mutations'], '-f', inputs['fasta']], monkeypatch)
    assert observed['mode'] == 'file'


@pytest.mark.parametrize('mode', ('directory', 'file'))
def test_rnafold_valid_inputs_still_prepare_jobs(inputs, monkeypatch, mode):
    observed = {}

    def pool_boundary(*args, **kwargs):
        import sys
        frame = sys._getframe(1).f_locals
        observed['jobs'] = frame['work_items']
        observed['intron_mapping'] = frame['args'].intron_premrna_mapping
        raise ValidationBoundary()

    arguments = [rnafold.__file__, '-i', str(inputs['root' if mode == 'directory' else 'fasta']),
                 '-o', str(inputs['output'])]
    if mode == 'file':
        arguments.extend(['-tm', str(inputs['mappings']['transcript'])])
    monkeypatch.setattr('sys.argv', arguments)
    monkeypatch.setattr(rnafold.concurrent.futures, 'ProcessPoolExecutor', pool_boundary)
    with pytest.raises(ValidationBoundary):
        rnafold.main()
    assert len(observed['jobs']) == 1
    assert observed['jobs'][0][7] == 'GENE'
    assert observed['intron_mapping'] is None
    assert not inputs['output'].exists()


@pytest.mark.parametrize('source', ('mutations', 'chromosome'))
def test_af3_file_mode_does_not_derive_missing_companions(inputs, monkeypatch, source):
    observed = {}

    def pipeline_boundary(**kwargs):
        import sys
        observed['args'] = sys._getframe(1).f_locals['args']
        raise ValidationBoundary()

    arguments = [af3.__file__, '-f', str(inputs['fasta']), '-o', str(inputs['output']),
                 '-pd', str(inputs['reference']), '-rm', str(inputs['reference']),
                 '-rs', str(inputs['reference']), '-mdi', str(inputs['models']), '-ch', '5']
    if source == 'mutations':
        arguments.extend(['-mu', str(inputs['mutations']), '-ts', '1'])
    else:
        arguments.extend(['-cm', str(inputs['mappings']['chromosome'])])
    monkeypatch.setattr('sys.argv', arguments)
    monkeypatch.setattr(af3, 'AlphaFold3Pipeline', pipeline_boundary)
    with pytest.raises(ValidationBoundary):
        af3.main()
    namespace = observed['args']
    assert namespace.mutations == (str(inputs['mutations']) if source == 'mutations' else None)
    assert namespace.chromosome_mapping == (str(inputs['mappings']['chromosome'])
                                             if source == 'chromosome' else None)
    assert namespace.premrna_mapping is None
    assert namespace.vcf is None
    assert not inputs['output'].exists()


def test_miranda_file_mode_does_not_reuse_transcript_mapping_as_intron_mapping(inputs, monkeypatch):
    observed = {}

    def output_boundary(path, **kwargs):
        import sys
        assert Path(path) == inputs['output']
        observed['args'] = sys._getframe(1).f_locals['args']
        raise ValidationBoundary()

    monkeypatch.setattr('sys.argv', [miranda.__file__, '-i', str(inputs['fasta']),
                                     '-M', str(inputs['mappings']['transcript']),
                                     '-o', str(inputs['output']), '-d', str(inputs['reference'])])
    monkeypatch.setattr('shutil.which', lambda executable: '/bin/true')
    monkeypatch.setattr(miranda.os, 'makedirs', output_boundary)
    with pytest.raises(ValidationBoundary):
        miranda.main()
    assert observed['args'].intron_premrna_mapping is None
    assert not inputs['output'].exists()
