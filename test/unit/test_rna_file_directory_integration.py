import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import textwrap

import pytest

from biofeaturefactory.lib.utility import mint_pkey


REPO = Path(__file__).resolve().parents[2]
TOKENS = ('A10C', 'A20G', 'A30ACCC')
WT_SEQUENCE = 'A' * 60


def write_file(path, content):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)
    return path


@pytest.fixture
def rna_inputs(tmp_path):
    root = tmp_path / 'input root'
    gene = root / 'TEST'
    fasta = write_file(gene / 'fastas' / 'TEST.fasta',
                       '>ORF\n' + WT_SEQUENCE + '\n>transcript\n' + 'C' * 6 + WT_SEQUENCE
                       + '\n>genomic\n' + 'G' * 9 + WT_SEQUENCE + '\n')
    mappings = {}
    for directory, filename, column, mapped_tokens in (
        ('transcript', 'transcript_mapping_TEST.csv', 'transcript', ('A16C', 'A26G', 'A36ACCC')),
        ('gDNA', 'genomic_mapping_TEST.csv', 'genomic', ('A19C', 'A29G', 'A39ACCC')),
    ):
        mappings[column] = write_file(
            gene / 'mappings' / directory / filename,
            f'pkey,mutant,{column}\n' + ''.join(
                f'{mint_pkey("TEST", token)},{token},{mapped}\n'
                for token, mapped in zip(TOKENS, mapped_tokens)))
    write_file(gene / 'CodonMSA' / 'TEST.codon.msa.fasta', '>focus\nCCCCCC\n>decoy\nCCCCCC\n')
    return {'root': root, 'fasta': fasta, **mappings}


def run_cli(module, arguments, env=None):
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    environment['PYTHONPATH'] = str(REPO) + os.pathsep + environment.get('PYTHONPATH', '')
    environment.update(env or {})
    result = subprocess.run([sys.executable, '-m', module, *map(str, arguments)],
                            cwd=REPO, env=environment, text=True,
                            capture_output=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    return result


def read_tables(directory):
    tables = {}
    for path in sorted(directory.glob('*.tsv')):
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        tables[path.name] = sorted(rows, key=lambda row: json.dumps(row, sort_keys=True))
    return tables


def assert_complete_tables(file_tables, directory_tables, count):
    assert file_tables == directory_tables
    assert len(file_tables) == count
    expected_keys = {mint_pkey('TEST', token) for token in TOKENS}
    for name, rows in file_tables.items():
        assert rows, name
        assert {row['pkey'] for row in rows} == expected_keys, name


def mode_sources(tmp_path, inputs, mode):
    if mode != 'detached':
        return inputs
    detached = tmp_path / 'detached files'
    detached.mkdir(exist_ok=True)
    sources = dict(inputs)
    for key in ('fasta', 'transcript', 'genomic', 'chromosome', 'mutations', 'vcf'):
        if key in inputs:
            path = detached / inputs[key].name
            shutil.copyfile(inputs[key], path)
            sources[key] = path
    return sources


def install_predictor(tmp_path, name):
    executable = tmp_path / 'binary with spaces' / name
    code = '''\
        import json
        import os
        from pathlib import Path
        import sys
        import uuid
        name = Path(sys.argv[0]).name
        fasta = Path(sys.argv[2] if name == 'miranda' else sys.argv[1])
        lines = fasta.read_text().splitlines()
        sequence = ''.join(lines[1:])
        record = {'id': lines[0][1:], 'sequence': sequence}
        Path(os.environ['BFF_TEST_CALLS'], str(uuid.uuid4()) + '.json').write_text(json.dumps(record))
        if name == 'miranda':
            score = 140 + 20 * sequence.count('C') + 10 * sequence.count('G')
            print(f'>miR-test target {score} -20.4 _ 1 20 10 30 20 100 100')
        else:
            assert Path(sys.argv[2], 'config_file').is_file()
            score = 8.25 + sequence.count('C') + 2 * sequence.count('G')
            print(f'10 11 {score} High donor')
    '''
    write_file(executable, f'#!{sys.executable}\n' + textwrap.dedent(code))
    executable.chmod(0o755)
    return executable


def expected_sequences(predictor):
    prefix = 'C' * 6 if predictor == 'miranda' else 'G' * 9
    sequences = [WT_SEQUENCE, 'A' * 9 + 'C' + 'A' * 50,
                 'A' * 19 + 'G' + 'A' * 40,
                 'A' * 30 + 'CCC' + 'A' * 30]
    return sorted(prefix + sequence for sequence in sequences)


@pytest.mark.parametrize('predictor', ['miranda', 'genesplicer'])
def test_native_rna_pipeline_file_directory_equivalence(tmp_path, rna_inputs, predictor):
    executable = install_predictor(tmp_path, predictor)
    reference = write_file(tmp_path / 'reference.fa', '>miR-test\nUUUUUUUUUUUUUUUUUUUU\n')
    model = write_file(tmp_path / 'model with spaces' / 'config_file', 'fixture model\n').parent
    tables = []
    invocations = []
    for mode in ('file', 'directory', 'detached'):
        sources = mode_sources(tmp_path, rna_inputs, mode)
        output = tmp_path / mode
        calls = tmp_path / f'{mode} calls'
        calls.mkdir()
        arguments = ['-i', sources['root' if mode == 'directory' else 'fasta'], '-o', output]
        if predictor == 'miranda':
            module = 'biofeaturefactory.miranda.miranda_ensemble'
            arguments += ['-m', executable.parent, '-d', reference, '-np']
            if mode != 'directory':
                arguments += ['-M', sources['transcript']]
            tool_directory = 'Miranda'
        else:
            module = 'biofeaturefactory.genesplicer.genesplicer_ensemble'
            arguments += ['-g', executable.parent, '--model-dir', model, '-wo', '1']
            if mode != 'directory':
                arguments += ['-m', sources['genomic']]
            tool_directory = 'GeneSplicer'
        run_cli(module, arguments, {'BFF_TEST_CALLS': str(calls)})
        tables.append(read_tables(output / 'TEST' / tool_directory))
        records = [json.loads(path.read_text()) for path in calls.glob('*.json')]
        assert sorted(record['sequence'] for record in records) == expected_sequences(predictor)
        invocations.append(sorted(records, key=lambda record: record['id']))
    assert invocations[0] == invocations[1] == invocations[2]
    for actual in tables[1:]:
        assert_complete_tables(tables[0], actual, count=3)
    events = {row['pkey']: row for row in tables[0]['TEST.events.tsv']}
    expected_deltas = (20.0, 10.0, 60.0) if predictor == 'miranda' else (1.0, 2.0, 3.0)
    expected_distances = (6, 16, 26) if predictor == 'miranda' else (9, 19, 29)
    score_column = 'delta_tot_score' if predictor == 'miranda' else 'dscore'
    for token, delta, distance in zip(TOKENS, expected_deltas, expected_distances):
        event = events[mint_pkey('TEST', token)]
        assert float(event[score_column]) == delta
        assert float(event['distance_to_snv']) == distance
        assert event['cls'] == 'strengthened'


def test_rnafold_file_directory_equivalence_with_real_vienna(tmp_path, rna_inputs):
    pytest.importorskip('RNA')
    tables = []
    for mode in ('file', 'directory', 'detached'):
        sources = mode_sources(tmp_path, rna_inputs, mode)
        output = tmp_path / mode
        arguments = ['-i', sources['root' if mode == 'directory' else 'fasta'],
                     '-o', output, '-w', '31', '-s', '10', '--workers', '1']
        if mode != 'directory':
            arguments += ['-tm', sources['transcript']]
        run_cli('biofeaturefactory.RNAfold.run_viennaRNA_pipeline', arguments)
        directory = output / 'TEST' / 'RNAfold'
        tables.append(read_tables(directory))
        report, = directory.glob('rnafold.run_summary.*.json')
        summary = json.loads(report.read_text())
        assert summary['mutations_successful'] == 3
        assert summary['mutations_unsuccessful'] == 0
    for actual in tables[1:]:
        assert_complete_tables(tables[0], actual, count=2)
    summaries = {row['pkey']: row for row in tables[0]['TEST.rnafold.tsv']}
    for token, coordinate in zip(TOKENS, (16, 26, 36)):
        summary = summaries[mint_pkey('TEST', token)]
        assert int(summary['transcript_pos']) == coordinate
        assert float(summary['ref_mfe_G']) == float(summary['alt_mfe_G']) == 0.0
    positions = tables[0]['TEST.rnafold.positions.tsv']
    assert len(positions) == 96
    inserted = [row for row in positions if row['align_status'] == 'inserted']
    assert len(inserted) == 3
    assert {row['pkey'] for row in inserted} == {mint_pkey('TEST', TOKENS[2])}
    assert all(row['tx_pos'] == row['delta_u'] == '' for row in inserted)


@pytest.fixture
def af3_inputs(tmp_path, rna_inputs):
    write_file(rna_inputs['fasta'], '>ORF\n' + WT_SEQUENCE + '\n>transcript\nCCCCCC\n')
    gene = rna_inputs['root'] / 'TEST'
    chromosome = write_file(
        gene / 'mappings' / 'chromosome' / 'chr_mapping_TEST.csv',
        'pkey,mutant,chromosome\n' + ''.join(
            f'{mint_pkey("TEST", token)},{token},{mapped}\n'
            for token, mapped in zip(TOKENS, ('A110C', 'A120G', 'A130ACCC'))))
    mutations = write_file(gene / 'mappings' / 'mutations' / 'TEST_mutations.csv',
                           'mutant\n' + '\n'.join(TOKENS) + '\n')
    vcf = write_file(gene / 'vcf' / 'TEST.vcf', '#CHROM\tPOS\tID\tREF\tALT\n5\t110\t.\tA\tC\n')
    postar = write_file(tmp_path / 'postar.bed', 'chr5\t100\t170\tfixture\t+\tRBP1\tCLIP\tcell\taccession\t10\n')
    mapping = write_file(tmp_path / 'rbp.tsv', 'Entry\tGene Names\tProtein names\tLength\nP00001\tRBP1\tfixture\t3\n')
    sequences = write_file(tmp_path / 'rbp.fa', '>sp|P00001|RBP1_HUMAN\nMKT\n')
    models = write_file(tmp_path / 'models' / 'af3.bin', 'external predictor fixture\n').parent
    executable = tmp_path / 'af3 binaries' / 'docker'
    code = '''\
        import json
        import os
        from pathlib import Path
        import sys
        import uuid
        arguments = sys.argv[1:]
        if arguments[:2] == ['image', 'inspect']:
            print('sha256:bff-integration-fixture')
            raise SystemExit(0)
        assert arguments[0] == 'run', arguments
        mounts = {}
        for position, argument in enumerate(arguments):
            if argument in ('-v', '--volume'):
                host, container, *options = arguments[position + 1].split(':')
                mounts[container] = Path(host)
        source = mounts['/root/af_input']
        destination = mounts['/root/af_output']
        json_path = next((argument.split('=', 1)[1] for argument in arguments
                          if argument.startswith('--json_path=')), None)
        paths = [source / Path(json_path).name] if json_path else sorted(source.glob('*.json'))
        for path in paths:
            payload = json.loads(path.read_text())
            Path(os.environ['BFF_TEST_CALLS'], str(uuid.uuid4()) + '.json').write_text(json.dumps(payload))
            output_name = payload['name']
            rna = payload['sequences'][0]['rna']['sequence']
            protein = payload['sequences'][1]['protein']['sequence']
            fields = ['group_PDB', 'auth_asym_id', 'auth_seq_id', 'auth_comp_id',
                      'auth_atom_id', 'Cartn_x', 'Cartn_y', 'Cartn_z']
            cif = 'data_fixture\\nloop_\\n' + ''.join('_atom_site.' + field + '\\n' for field in fields)
            atom_name = chr(34) + 'C1' + chr(39) + chr(34)
            for position, residue in enumerate(rna, 1):
                cif += f'ATOM R {position} {residue} {atom_name} {position / 10} 0 0\\n'
            for position in range(1, len(protein) + 1):
                cif += f'ATOM P {position} LYS CA {position / 10} 1 0\\n'
            size = len(rna) + len(protein)
            confidences = {'atom_plddts': [90.0] * size, 'pae': [[2.0] * size for unused in range(size)]}
            pae = 2.0 + rna.count('C') + rna.count('G')
            summary = {'ranking_score': 0.9, 'ptm': 0.8, 'iptm': 0.8,
                       'chain_pair_pae_min': [[0.0, pae], [pae, 0.0]]}
            result_root = destination / output_name
            for directory in [result_root, result_root / 'seed-1_sample-0', result_root / 'seed-1_sample-1']:
                directory.mkdir(parents=True, exist_ok=True)
                (directory / f'{output_name}_model.cif').write_text(cif)
                (directory / f'{output_name}_confidences.json').write_text(json.dumps(confidences))
                (directory / f'{output_name}_summary_confidences.json').write_text(json.dumps(summary))
    '''
    write_file(executable, f'#!{sys.executable}\n' + textwrap.dedent(code))
    executable.chmod(0o755)
    return {**rna_inputs, 'chromosome': chromosome, 'mutations': mutations, 'vcf': vcf,
            'postar': postar, 'rbp_mapping': mapping, 'rbp_sequences': sequences,
            'models': models, 'executable': executable}


def af3_arguments(inputs, mode, output):
    arguments = ['--postar-db', inputs['postar'], '--rbp-mapping', inputs['rbp_mapping'],
                 '-rs', inputs['rbp_sequences'], '-f', inputs['root' if mode == 'directory' else 'fasta'],
                 '-o', output, '-ws', '31', '-rw', '0']
    for flag, name in [('-cm', 'chromosome'), ('-mu', 'mutations'), ('-v', 'vcf')]:
        arguments += [flag, inputs['root'] if mode == 'directory' else inputs[name]]
    return arguments


def assert_af3_payloads(records):
    expected_rna = {'A' * 25, 'A' * 9 + 'C' + 'A' * 15, 'A' * 31,
                    'A' * 15 + 'G' + 'A' * 15,
                    'A' * 16 + 'CCC' + 'A' * 15}
    assert {record['sequences'][0]['rna']['sequence'] for record in records} == expected_rna
    for record in records:
        assert record['sequences'][0]['rna']['id'] == 'R'
        assert record['sequences'][1]['protein']['id'] == 'P'
        assert record['sequences'][1]['protein']['sequence'] == 'MKT'


@pytest.mark.parametrize('driver', ['local', 'burst'])
def test_af3_file_directory_equivalence_through_prediction_and_exports(tmp_path, af3_inputs, driver):
    tables = []
    manifests = []
    for mode in ('file', 'directory', 'detached'):
        sources = mode_sources(tmp_path, af3_inputs, mode)
        output = tmp_path / mode
        calls = tmp_path / f'{mode} calls'
        calls.mkdir()
        environment = {'PATH': str(af3_inputs['executable'].parent) + os.pathsep + os.environ['PATH'],
                       'BFF_TEST_CALLS': str(calls)}
        arguments = af3_arguments(sources, mode, output)
        if driver == 'local':
            run_cli('biofeaturefactory.alphafold3.alphafold3_pipeline',
                    arguments + ['-mdi', af3_inputs['models'], '-mg', '1'], environment)
        else:
            run_cli('biofeaturefactory.alphafold3.burst',
                    ['submit', *arguments, '--model-dir', af3_inputs['models'], '--no-submit'], environment)
            submission, = (output / '.burst' / 'submissions').iterdir()
            manifest_path = submission / 'manifest.tsv'
            with manifest_path.open() as handle:
                manifest = list(csv.DictReader((line for line in handle if not line.startswith('#')), delimiter='\t'))
            assert len(manifest) == 5
            before = manifest_path.read_bytes()
            for index, entry in enumerate(manifest):
                assert entry['array_idx'] == str(index)
                task_environment = dict(os.environ, **environment,
                                        PYTHONDONTWRITEBYTECODE='1', SLURM_ARRAY_TASK_ID=str(index),
                                        CUDA_VISIBLE_DEVICES='0', PYTHONPATH=str(REPO))
                result = subprocess.run(['bash', str(submission / 'run.slurm')],
                                        cwd=REPO, env=task_environment, capture_output=True,
                                        text=True, timeout=60)
                assert result.returncode == 0, result.stdout + result.stderr
            assert manifest_path.read_bytes() == before
            run_cli('biofeaturefactory.alphafold3.burst', ['ingest', *arguments], environment)
            manifests.append([{key: value for key, value in entry.items()
                               if key not in ('input_json_path', 'cache_dir')} for entry in manifest])
        records = [json.loads(path.read_text()) for path in calls.glob('*.json')]
        assert_af3_payloads(records)
        tables.append(read_tables(output / 'TEST' / 'AlphaFold3'))
    if manifests:
        assert manifests[0] == manifests[1] == manifests[2]
    for actual in tables[1:]:
        assert_complete_tables(tables[0], actual, count=3)
    events = {row['pkey']: row for row in tables[0]['TEST.events.tsv']}
    for token, pae_delta, wt_contacts, mut_contacts in zip(TOKENS, (1, 1, 3), (75, 93, 93), (75, 93, 102)):
        event = events[mint_pkey('TEST', token)]
        assert float(event['delta_chain_pair_pae_min']) == pae_delta
        assert int(event['wt_interface_contacts']) == wt_contacts
        assert int(event['mut_interface_contacts']) == mut_contacts
        assert int(event['n_samples_wt']) == int(event['n_samples_mut']) == 2
    sites = tables[0]['TEST.sites.tsv']
    assert len(sites) == 195
    assert all(float(row['contact_frequency']) == 1.0 for row in sites)
    inserted = [row for row in sites if row['align_status'] == 'inserted']
    assert len(inserted) == 3
    assert {row['pkey'] for row in inserted} == {mint_pkey('TEST', TOKENS[2])}
    assert {int(row['res_id']) for row in inserted} == {17, 18, 19}
