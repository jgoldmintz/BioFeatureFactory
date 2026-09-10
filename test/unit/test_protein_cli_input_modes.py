import argparse
import ast
from pathlib import Path
import sys

import pytest

from biofeaturefactory.lib.utility import InputPathAction, validate_input_mode


REPO = Path(__file__).resolve().parents[2]
PIPELINES = (
    ('NetSurfP3/netsurfp3_pipeline.py', '-m', ['-M', 'checkpoint.pth', '-c', 'config.yml']),
    ('netMHC/netmhc_pipeline.py', '-m', []),
    ('netNglyc/netnglyc_pipeline.py', '-md', []),
    ('netphos/netphos_pipeline.py', '-md', []),
)


class ValidatedInput(Exception):
    pass


@pytest.fixture
def inputs(tmp_path):
    root = tmp_path / 'out'
    (root / 'GENE/fastas').mkdir(parents=True)
    (root / 'GENE/mappings/mutations').mkdir(parents=True)
    fasta = root / 'GENE/fastas/GENE.fasta'
    mutation_file = root / 'GENE/mappings/mutations/GENE_mutations.csv'
    fasta.write_text('>ORF\nATGGCTTAA\n')
    mutation_file.write_text('mutant\nG4A\n')
    return root, fasta, mutation_file


def invoke_main(relative, arguments, monkeypatch):
    path = REPO / 'biofeaturefactory' / relative
    tree = ast.parse(path.read_text())
    main = next(node for node in tree.body if isinstance(node, ast.FunctionDef)
                and node.name == 'main')

    def stop_after_validation(parser, namespace, required_file_inputs=()):
        mode = validate_input_mode(parser, namespace, required_file_inputs)
        raise ValidatedInput(mode)

    namespace = {
        'argparse': argparse,
        'InputPathAction': InputPathAction,
        'validate_input_mode': stop_after_validation,
    }
    exec(compile(ast.Module(body=[main], type_ignores=[]), str(path), 'exec'), namespace)
    monkeypatch.setattr(sys, 'argv', [str(path), *arguments])
    namespace['main']()


@pytest.mark.parametrize('relative,companion,extra', PIPELINES)
def test_explicit_file_requires_companion(relative, companion, extra, inputs, monkeypatch, capsys):
    root, fasta, _ = inputs
    with pytest.raises(SystemExit) as error:
        invoke_main(relative, ['-i', str(fasta), '-o', str(root), *extra], monkeypatch)
    assert error.value.code == 2
    assert 'required in file mode' in capsys.readouterr().err


@pytest.mark.parametrize('relative,companion,extra', PIPELINES)
@pytest.mark.parametrize('input_mode', ['file', 'directory'])
def test_valid_inputs_pass_without_predictor_dependencies(relative, companion, extra, inputs,
                                                         monkeypatch, input_mode):
    root, fasta, mutation_file = inputs
    arguments = ['-i', str(fasta if input_mode == 'file' else root), '-o', str(root), *extra]
    if input_mode == 'file':
        arguments += [companion, str(mutation_file)]
    with pytest.raises(ValidatedInput) as result:
        invoke_main(relative, arguments, monkeypatch)
    assert result.value.args == (input_mode,)


@pytest.mark.parametrize('relative,companion,extra', PIPELINES)
@pytest.mark.parametrize('companion_first', [True, False])
def test_mixed_inputs_rejected_in_both_orders(relative, companion, extra, inputs,
                                              monkeypatch, capsys, companion_first):
    root, fasta, _ = inputs
    primary = ['-i', str(fasta)]
    secondary = [companion, str(root)]
    arguments = secondary + primary if companion_first else primary + secondary
    with pytest.raises(SystemExit) as error:
        invoke_main(relative, [*arguments, '-o', str(root), *extra], monkeypatch)
    assert error.value.code == 2
    expected = f'{companion} selected directory mode' if companion_first else '-i selected file mode'
    assert expected in capsys.readouterr().err


@pytest.mark.parametrize('relative,companion,extra', PIPELINES)
@pytest.mark.parametrize('nested', ['GENE', 'GENE/fastas'])
def test_nested_directories_rejected_before_backend(relative, companion, extra, inputs,
                                                     monkeypatch, capsys, nested):
    root, _, _ = inputs
    with pytest.raises(SystemExit) as error:
        invoke_main(relative, ['-i', str(root / nested), '-o', str(root), *extra], monkeypatch)
    assert error.value.code == 2
    assert f'Provide {root} instead.' in capsys.readouterr().err


@pytest.mark.parametrize('relative,companion,extra', PIPELINES)
def test_legacy_positional_input_is_validated(relative, companion, extra, inputs, monkeypatch):
    root, fasta, mutation_file = inputs
    with pytest.raises(ValidatedInput) as result:
        invoke_main(relative, [str(fasta), str(root), companion, str(mutation_file), *extra], monkeypatch)
    assert result.value.args == ('file',)
