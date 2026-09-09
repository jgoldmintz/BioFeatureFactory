"""Regressions for the nine AF3 findings reconciled in September 2026."""

import csv
import json
import subprocess
import sys
import time
from concurrent.futures import Future
from contextlib import nullcontext
from pathlib import Path
from types import SimpleNamespace

import pytest

from biofeaturefactory.alphafold3 import alphafold3_pipeline as local, burst
from biofeaturefactory.alphafold3.bin import burst_manifest as cache
from biofeaturefactory.alphafold3.bin import af3_runner
from biofeaturefactory.alphafold3.bin.af3_parser import (
    AF3Parser, analyze_binding, ensemble_interface_sites,
)
from biofeaturefactory.alphafold3.bin.binding_metrics import (
    BindingMetrics, RnaEditSpan, ThresholdConfig, compute_delta_metrics,
    compute_window_delta, format_sites_rows, format_events_rows,
    aggregate_mutation_summary, qc_flag_for_deltas,
)
from test_af3_runner import FakeAF3, _input, runner_factory


def write_structure(directory, name, rna_positions=(0, 30, 30), atom='"C1\'"'):
    directory.mkdir(parents=True, exist_ok=True)
    headers = ['group_PDB', 'auth_asym_id', 'auth_seq_id', 'auth_comp_id',
               'auth_atom_id', 'Cartn_x', 'Cartn_y', 'Cartn_z']
    rows = ['data_test', 'loop_'] + [f'_atom_site.{header}' for header in headers]
    for residue, position in enumerate(rna_positions, 1):
        rows.extend([f'ATOM R {residue} A P 90 0 0',
                     f'ATOM R {residue} A {atom} {position} 0 0'])
    rows.extend(['ATOM P 1 ALA CA 4 0 0', '#'])
    (directory / f'{name}_model.cif').write_text('\n'.join(rows))
    (directory / f'{name}_confidences.json').write_text(json.dumps({
        'atom_plddts': [90] * (2 * len(rna_positions) + 1), 'pae': [[0, 2], [2, 0]],
    }))
    (directory / f'{name}_summary_confidences.json').write_text(json.dumps({
        'chain_pair_pae_min': [[0, 2], [2, 0]],
    }))


def make_pipeline(tmp_path):
    pipeline = local.AlphaFold3Pipeline.__new__(local.AlphaFold3Pipeline)
    pipeline.threshold_config = ThresholdConfig()
    pipeline.summary_rows = []
    pipeline.events_rows = []
    pipeline.sites_rows = []
    pipeline.n_failed_mutations = 0
    pipeline.output_dir = tmp_path
    pipeline.validation_log = None
    pipeline.multi_window = False
    pipeline.window_size = 9
    pipeline.rbp_window = 50
    return pipeline


def metric(pae=2):
    return BindingMetrics('RBP', pae, 9, 90, 90, True)


def resolved(value):
    future = Future()
    future.set_result(value)
    return future


@pytest.mark.parametrize('atom', ['"C1\'"', "'C1''", "C1'"])
def test_cif_quotes_preserve_terminal_apostrophe_and_ca(tmp_path, atom):
    write_structure(tmp_path, 'test', atom=atom)
    structure = AF3Parser(str(tmp_path)).parse()
    assert list(structure.get_chain('R')[0].atoms) == ['P', "C1'"]
    assert structure.get_chain('R')[0].ca_coord.x == 0
    assert structure.get_chain('P')[0].ca_coord.x == 4
    assert all(residue.plddt == 90 for residue in structure.residues)
    assert analyze_binding(structure).n_contacts == 1
    assert analyze_binding(structure).chain_pair_pae_min == 2


def test_ranked_union_and_zero_frequencies_match_burst_and_local(tmp_path):
    root = tmp_path / 'hash'
    output = root / 'model'
    write_structure(output, 'model', (30, 0, 14))
    write_structure(output / 'seed-1_sample-0', 'sample', (0, 30, 30))
    write_structure(output / 'seed-1_sample-1', 'sample', (30, 0, 30))
    pipeline = make_pipeline(tmp_path)
    parsed = pipeline._parse_af3_output(output, 'RBP')
    local_sites = ensemble_interface_sites(parsed.structures, parsed.aggregation, parsed.ranked)
    burst_metrics, burst_sites, aggregation = burst._parse_cache_entry(root, 'model', 'RBP', ThresholdConfig())
    assert burst_metrics == parsed.metrics
    assert burst_sites == local_sites
    span = RnaEditSpan(1, 1, 1, 3, 3)
    local_rows = format_sites_rows('GENE-key', 'RBP', 'WT', local_sites,
                                  parsed.aggregation.contact_frequency_rna,
                                  parsed.aggregation.contact_frequency_protein, edit_span=span)
    burst_rows = format_sites_rows('GENE-key', 'RBP', 'WT', burst_sites,
                                  aggregation.contact_frequency_rna,
                                  aggregation.contact_frequency_protein, edit_span=span)
    assert local_rows == burst_rows
    rna_rows = [row for row in local_rows if row['chain'] == 'R']
    assert [row['contact_frequency'] for row in rna_rows] == [0.5, 0.5, 0.0]
    assert [row['is_contact'] for row in rna_rows] == [0, 1, 0]
    unknown = format_sites_rows('GENE-key', 'RBP', 'WT', local_sites, edit_span=span)
    assert all(row['contact_frequency'] == '' for row in unknown)


@pytest.mark.parametrize('wt_success,mut_success,paired,qc', [
    ({0, 1, 2}, {0, 1, 2}, 3, 'PASS'),
    ({0}, {2}, 0, 'PARTIAL'),
    ({0, 1}, {1, 2}, 1, 'PARTIAL'),
    ({0}, set(), 0, 'PARTIAL'),
    (set(), set(), 0, 'ALL_FAILED'),
])
def test_live_collector_only_compares_matching_windows(tmp_path, monkeypatch, wt_success, mut_success, paired, qc):
    pipeline = make_pipeline(tmp_path)
    context = local.MutationContext('GENE-key', 'GENE', 'A4G', 4, 'A', 'G',
                                    'A' * 20, 'A' * 9, 'AAAGAAAAA', 3)
    pending = pipeline._PendingRBPAnalysis('RBP', [], 0, n_windows=3,
                                          windows=[('A' * 9, 'G' * 9, index + 1) for index in range(3)],
                                          protein_msa='provided')
    pending.window_wt_futures = [resolved(SimpleNamespace(status='completed', result_path=('WT', index))) for index in range(3)]
    pending.window_mut_futures = [resolved(SimpleNamespace(status='completed', result_path=('MUT', index))) for index in range(3)]

    def parse_output(path, _rbp):
        allele, index = path
        if index not in (wt_success if allele == 'WT' else mut_success):
            return None
        return local._ParsedResult(metric(9 - 3 * index), [], None)

    monkeypatch.setattr(pipeline, '_parse_af3_output', parse_output)
    delta = pipeline._collect_rbp_results(context, pending)
    assert delta.n_windows_paired == paired
    assert delta.n_windows_success_wt == len(wt_success)
    assert delta.n_windows_success_mut == len(mut_success)
    assert delta.delta_chain_pair_pae_min == 0
    assert delta.event_class.value == ('unchanged' if paired else 'incomplete')
    assert qc_flag_for_deltas([delta]) == qc
    pipeline._finalize_mutation_results(context, [delta])
    assert pipeline.summary_rows[0]['qc_flags'] == qc
    assert pipeline.n_failed_mutations == (qc != 'PASS')
    assert pipeline.events_rows[0]['protein_msa'] == 'provided'


def test_all_local_windows_emit_distinct_site_rows(tmp_path, monkeypatch):
    output = tmp_path / 'raw'
    write_structure(output, 'model')
    pipeline = make_pipeline(tmp_path)
    context = local.MutationContext('GENE-key', 'GENE', 'A4AG', 4, 'A', 'AG', 'A' * 20,
                                    'A' * 9, 'A' * 10, 3)
    pending = pipeline._PendingRBPAnalysis('RBP', [], 0, n_windows=2,
                                          windows=[('A' * 9, 'A' * 10, 1), ('A' * 9, 'A' * 10, 5)])
    pending.window_wt_futures = [resolved(SimpleNamespace(status='completed', result_path=output)) for _ in range(2)]
    pending.window_mut_futures = [resolved(SimpleNamespace(status='completed', result_path=output)) for _ in range(2)]
    pipeline._collect_rbp_results(context, pending)
    assert {(row['allele'], row['window_idx'], row['window_edit_offset']) for row in pipeline.sites_rows} == {
        ('WT', 0, 1), ('MUT', 0, 1), ('WT', 1, 5), ('MUT', 1, 5),
    }
    keys = [(row['allele'], row['window_idx'], row['chain'], row['res_id']) for row in pipeline.sites_rows]
    assert len(keys) == len(set(keys)) == 8


@pytest.mark.parametrize('outcome,expected_exit', [('failed', 1), ('exception', 1), ('partial', 1), ('success', 0), ('skip', 0)])
def test_main_counts_real_finalization_across_gene_flushes(tmp_path, monkeypatch, outcome, expected_exit):
    inputs = tmp_path / 'inputs'
    outputs = tmp_path / 'outputs'
    models = tmp_path / 'models'
    models.mkdir()
    for gene in ('GENEA', 'GENEB'):
        fasta = inputs / gene / 'fastas' / f'{gene}.fasta'
        fasta.parent.mkdir(parents=True)
        fasta.write_text('>ORF\nAAAAAAAAA\n')
        mutations = inputs / gene / 'mappings' / 'mutations' / f'{gene}_mutations.csv'
        mutations.parent.mkdir(parents=True)
        mutations.write_text('mutations\nA4G\n')
    pipeline = make_pipeline(outputs)
    pipeline.rbp_db = SimpleNamespace(group_by_rbp=lambda _sites: {'RBP': []})
    pipeline._get_nearby_rbps = lambda *_args: [] if outcome == 'skip' else [object()]

    def submit(context, rbp_name, sites, **_kwargs):
        if outcome == 'exception':
            raise RuntimeError('mock external failure')
        return pipeline._PendingRBPAnalysis(rbp_name, sites, 0,
            wt_future=resolved(SimpleNamespace(status='completed' if outcome in ('partial', 'success') else 'failed', result_path='WT')),
            mut_future=resolved(SimpleNamespace(status='completed' if outcome == 'success' else 'failed', result_path='MUT')))

    pipeline._submit_rbp_jobs = submit
    pipeline._parse_af3_output = lambda *_args: local._ParsedResult(metric(), [], None)
    pipeline.af3_runner = SimpleNamespace(batch_submissions=nullcontext, shutdown=lambda: None)
    monkeypatch.setattr(local, 'AlphaFold3Pipeline', lambda **_kwargs: pipeline)
    monkeypatch.setattr(sys, 'argv', ['af3', '-f', str(inputs), '-o', str(outputs),
                                    '-pd', 'postar', '-rm', 'rbp', '-rs', 'seq', '-mdi', str(models),
                                    '-mu', str(inputs), '-ch', 'chr1', '-ts', '100'])
    assert local.main() == expected_exit
    assert pipeline.n_failed_mutations == 2 * expected_exit
    assert not pipeline.summary_rows
    rows = []
    for gene in ('GENEA', 'GENEB'):
        with (outputs / gene / 'AlphaFold3' / f'{gene}.tsv').open() as handle:
            rows.extend(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 2
    expected_qc = {'failed': 'ALL_FAILED', 'exception': 'FAILED:', 'partial': 'PARTIAL',
                   'success': 'PASS', 'skip': 'no_rbps_in_region'}[outcome]
    assert all(row['qc_flags'].startswith(expected_qc) for row in rows)
    if expected_exit:
        pipeline._record_failed_mutation(rows[0]['pkey'])
        assert pipeline.n_failed_mutations == 2


def test_msa_provenance_matches_submitted_inputs_and_summary(tmp_path):
    pipeline = make_pipeline(tmp_path)
    context = local.MutationContext('GENE-key', 'GENE', 'A4G', 4, 'A', 'G', 'A' * 20,
                                    'A' * 9, 'AAAGAAAAA', 3)
    submitted = []
    def submit(af3_input, **_kwargs):
        submitted.append(af3_input)
        return resolved(SimpleNamespace(status='completed', result_path='raw'))
    pipeline.af3_runner = SimpleNamespace(submit_job_async=submit)
    pipeline._parse_af3_output = lambda *_args: local._ParsedResult(metric(), [], None)
    deltas = []
    for msa in ('>query\nMPEP\n', None):
        pipeline.seq_mapper = SimpleNamespace(get_rbp_data=lambda _rbp: SimpleNamespace(sequence='MPEP', msa_content=msa))
        pending = pipeline._submit_rbp_jobs(context, 'RBP', [])
        deltas.append(pipeline._collect_rbp_results(context, pending))
        assert pending.protein_msa == ('provided' if submitted[-1].protein_msa else 'none')
    assert [row['protein_msa'] for row in format_events_rows('GENE-key', deltas)] == ['provided', 'none']
    summary = aggregate_mutation_summary(deltas)
    assert (summary['protein_msa'], summary['n_rbps_msa_provided'], summary['n_rbps_msa_free']) == ('mixed', 1, 1)


def test_benign_stdout_does_not_amplify_incomplete_batch(runner_factory, monkeypatch):
    runner = runner_factory(batch_size=3)
    def plan(index, names):
        return (set(names[:1]), 1, 'example: permission denied; Error response from daemon') if index == 0 else (set(names), 0, '')
    fake = FakeAF3(plan)
    monkeypatch.setattr(af3_runner.subprocess, 'run', fake)
    with runner.batch_submissions():
        futures = [runner.submit_job_async(_input(f'job{index}'), job_id=f'job{index}') for index in range(3)]
    assert [future.result(timeout=5).status for future in futures] == ['completed'] * 3
    assert len(fake.calls) == 3
    assert runner._docker_failed is False


def test_docker_daemon_stderr_still_aborts(runner_factory, monkeypatch):
    runner = runner_factory()
    fake = FakeAF3(lambda _index, _names: (set(), 125, ''))
    def run(command, **kwargs):
        result = fake(command, **kwargs)
        if command[:2] == ['docker', 'run']:
            kwargs['stderr'].write('docker: Error response from daemon: unavailable\n')
        return result
    monkeypatch.setattr(af3_runner.subprocess, 'run', run)
    assert runner.submit_job(_input('job'), job_id='job').status == 'failed'
    assert runner._docker_failed


def test_clear_preserves_staging_and_invalidates_old_worker(tmp_path):
    root = tmp_path / '.cache' / 'af3'
    destination = root / 'hash'
    stage = cache.begin_cache_stage(destination)
    write_structure(stage / 'model', 'model')
    lock_path = cache.cache_coordination_dir(root) / '.bff-af3-cache.lock'
    inode = lock_path.stat().st_ino
    cache.clear_cache(root)
    assert stage.is_dir()
    assert lock_path.stat().st_ino == inode
    with pytest.raises(ValueError, match='generation was cleared'):
        cache.publish_cache(stage, destination, 'model')
    with pytest.raises(ValueError, match='generation was cleared'):
        cache.begin_cache_stage(destination, '0')
    assert stage.is_dir() and not destination.exists()
    fresh = cache.begin_cache_stage(destination, cache.cache_generation(root))
    write_structure(fresh / 'model', 'model')
    cache.publish_cache(fresh, destination, 'model')
    assert cache.is_cache_complete(destination, 'model')


@pytest.mark.parametrize('operation', ['stage', 'publish'])
def test_clear_serializes_with_other_processes(tmp_path, operation):
    root = tmp_path / '.cache' / 'af3'
    destination = root / 'hash'
    stage = cache.begin_cache_stage(destination)
    write_structure(stage / 'model', 'model')
    helper = Path(cache.__file__)
    if operation == 'publish':
        command = [sys.executable, str(helper), 'publish', str(stage), str(destination), 'model']
    else:
        command = [sys.executable, str(helper), 'stage', str(destination), '0']
    with cache.cache_lock(root):
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        with pytest.raises(subprocess.TimeoutExpired):
            process.communicate(timeout=0.3)
    stdout, stderr = process.communicate(timeout=10)
    assert process.returncode == 0, stderr
    staged_after = Path(stdout.strip()) if operation == 'stage' else stage
    cache.clear_cache(root)
    assert not destination.exists()
    if operation == 'stage':
        assert staged_after.is_dir()


def test_saved_hnrnpc_contact_control_read_only():
    raw_root = Path(__file__).parents[3] / 'non-repo-files' / 'out' / 'af3_runs'
    pipeline = make_pipeline(Path('.'))
    observed = []
    for allele in ('WT', 'MUT'):
        output = raw_root / f'NPM1-4ce2f2ef6b82_HNRNPC_{allele}' / 'output'
        if not output.is_dir():
            pytest.skip('Retained HNRNPC outputs unavailable')
        parsed = pipeline._parse_af3_output(output, 'HNRNPC')
        contacts = [analyze_binding(structure).n_contacts for structure in parsed.structures]
        observed.append(parsed.metrics)
        assert contacts == ([19, 19, 22, 17, 19] if allele == 'WT' else [21, 19, 19, 19, 21])
    assert [item.interface_contacts for item in observed] == [19, 19]
    assert [round(item.chain_pair_pae_min, 2) for item in observed] == [1.92, 2.24]
    assert compute_delta_metrics('HNRNPC', *observed).event_class.value == 'unchanged'


def input_args(root):
    return SimpleNamespace(fasta=str(root), mutations=str(root), chromosome_mapping=str(root),
                           vcf=str(root), chrom=None, validation_log=None)


def test_nested_real_inputs_match_explicit_sources(tmp_path, monkeypatch):
    root = Path(__file__).parents[3] / 'non-repo-files' / 'out'
    if not (root / 'NPM1').is_dir():
        pytest.skip('Retained per-gene inputs unavailable')
    recorded = []
    def record(**kwargs):
        recorded.append((kwargs['gene_name'], kwargs['mutation'], kwargs['chrom'], kwargs['chrom_mapping']))
        return []
    monkeypatch.setattr(burst, '_iterate_mutation', record)
    skipped = []
    list(burst.iterate_inputs(input_args(root), object(), object(), skipped))
    nested = list(recorded)
    recorded.clear()
    for gene in ('NPM1', 'PAM'):
        args = input_args(root)
        args.fasta = str(root / gene / 'fastas' / f'{gene}.fasta')
        args.mutations = str(root / gene / 'mappings' / 'mutations' / f'{gene}_mutations.csv')
        args.chromosome_mapping = str(root / gene / 'mappings' / 'chromosome' / f'chr_mapping_{gene}.csv')
        args.vcf = str(root / gene / 'vcf' / f'{gene}.vcf')
        list(burst.iterate_inputs(args, object(), object(), []))
    assert nested == recorded
    assert len(nested) == 8 and not skipped
    assert [sum(row[0] == gene for row in nested) for gene in ('NPM1', 'PAM')] == [6, 2]
    assert all(row[2] and row[1] in row[3] for row in nested)


def test_mixed_input_layout_avoids_aa_and_prefix_decoys(tmp_path, monkeypatch):
    gene_dir = tmp_path / 'GENE'
    (gene_dir / 'fastas').mkdir(parents=True)
    (gene_dir / 'fastas' / 'GENE.fasta').write_text('>ORF\nAAAAAAAAA\n')
    (tmp_path / 'GENE2.fasta').write_text('>ORF\nAAAAAAAAA\n')
    mutations = gene_dir / 'mappings' / 'mutations'
    mutations.mkdir(parents=True)
    (mutations / 'GENE_mutations.csv').write_text('mutant\nA4G\n')
    (tmp_path / 'GENE2_mutations.csv').write_text('mutant\nA5T\n')
    (tmp_path / 'GENE_aa_mapping.csv').write_text('mutant,aa\nA4G,K2R\n')
    args = input_args(tmp_path)
    args.chrom = 'chr1'
    args.vcf = None
    args.chromosome_mapping = None
    recorded = []
    monkeypatch.setattr(burst, '_iterate_mutation', lambda **kwargs: recorded.append(
        (kwargs['gene_name'], kwargs['mutation'])) or [])
    skipped = []
    list(burst.iterate_inputs(args, object(), object(), skipped))
    assert recorded == [('GENE', 'A4G'), ('GENE2', 'A5T')]
    (mutations / 'GENE_mutations.csv').unlink()
    recorded.clear()
    list(burst.iterate_inputs(args, object(), object(), skipped))
    assert recorded == [('GENE2', 'A5T')]
    assert [(item.gene, item.pkey, item.qc_flag) for item in skipped] == [('GENE', '', 'FAILED:no_mutations_source')]


def test_mixed_chromosome_layout_uses_exact_gene_mapping(tmp_path, monkeypatch):
    fasta = tmp_path / 'GENE' / 'fastas' / 'GENE.fasta'
    fasta.parent.mkdir(parents=True)
    fasta.write_text('>ORF\nAAAAAAAAA\n')
    (tmp_path / 'GENE2.fasta').write_text('>ORF\nAAAAAAAAA\n')
    nested_map = tmp_path / 'GENE' / 'mappings' / 'chromosome' / 'chr_mapping_GENE.csv'
    nested_map.parent.mkdir(parents=True)
    nested_map.write_text('pkey,mutant,orf,chromosome\nGENE-key,A4G,A4G,A104G\n')
    (tmp_path / 'chr_mapping_GENE2.csv').write_text('pkey,mutant,orf,chromosome\nGENE2-key,A5T,A5T,A205T\n')
    (tmp_path / 'chr_mapping_GENE22.csv').write_text('pkey,mutant,orf,chromosome\nGENE22-key,A6C,A6C,A306C\n')
    (tmp_path / 'GENE2_aa_mapping.csv').write_text('mutant,aa\nA5T,K2I\n')
    args = input_args(tmp_path)
    args.mutations = None
    args.vcf = None
    args.chrom = 'chr1'
    recorded = []
    monkeypatch.setattr(burst, '_iterate_mutation', lambda **kwargs: recorded.append(
        (kwargs['gene_name'], kwargs['mutation'], kwargs['chrom_mapping'])) or [])
    list(burst.iterate_inputs(args, object(), object(), []))
    assert recorded == [('GENE', 'A4G', {'A4G': 'A104G'}),
                        ('GENE2', 'A5T', {'A5T': 'A205T'})]


@pytest.mark.parametrize('wt_success,mut_success,paired,qc', [
    ({0, 1}, {0, 1}, 2, 'PASS'), ({0}, {1}, 0, 'PARTIAL'),
    ({0, 1}, {1}, 1, 'PARTIAL'), (set(), set(), 0, 'ALL_FAILED'),
])
def test_burst_pairing_sites_and_exit_match_local(tmp_path, monkeypatch, wt_success, mut_success, paired, qc):
    output = tmp_path / 'out'
    (output / '.cache' / 'af3').mkdir(parents=True)
    args = SimpleNamespace(output=str(output), postar_db='postar', rbp_mapping='rbp', rbp_sequences=None, msa_dir='msa')
    inputs = []
    parsed = {}
    sites = [SimpleNamespace(chain='R', res_id=1, res_name='A', plddt=90,
                             is_contact=True, min_contact_distance=4)]
    for index in range(2):
        for allele, successes in (('WT', wt_success), ('MUT', mut_success)):
            af3_input = af3_runner.AF3Input(f'{allele}{index}', 'A' * (4 + index) if allele == 'WT' else 'G' * (4 + index),
                                           'MPEP', protein_msa='>query\nMPEP\n')
            item = burst.BurstInput('GENE', 'GENE-key', 'RBP', allele, index, af3_input,
                                    edit_span=RnaEditSpan(index + 1, 1, 1, 4 + index, 4 + index))
            inputs.append(item)
            parsed[burst._cache_output_name(item.input_hash)] = (metric(2 + index), sites, None) if index in successes else (None, [], None)
    monkeypatch.setattr(burst, 'POSTAR3Database', lambda _path: object())
    monkeypatch.setattr(burst, 'RBPSequenceMapper', lambda **_kwargs: object())
    monkeypatch.setattr(burst, '_ingest_warn_in_flight_slurm', lambda _path: None)
    monkeypatch.setattr(burst, 'iterate_inputs', lambda *_args, **_kwargs: inputs)
    monkeypatch.setattr(burst, '_parse_cache_entry', lambda _cache, output_name, **_kwargs: parsed[output_name])
    assert burst.cmd_ingest(args) == (qc != 'PASS')
    with (output / 'GENE' / 'AlphaFold3' / 'GENE.events.tsv').open() as handle:
        row = next(csv.DictReader(handle, delimiter='\t'))
    assert (row['qc_flags'], int(row['n_windows_paired']), row['protein_msa']) == (qc, paired, 'provided')
    assert (float(row['delta_chain_pair_pae_min']) == 0) if paired else (row['delta_chain_pair_pae_min'] == '')
    path = output / 'GENE' / 'AlphaFold3' / 'GENE.sites.tsv'
    if wt_success or mut_success:
        with path.open() as handle:
            site_rows = list(csv.DictReader(handle, delimiter='\t'))
        assert {(row['allele'], int(row['window_idx'])) for row in site_rows} == {
            *[('WT', index) for index in wt_success], *[('MUT', index) for index in mut_success],
        }


@pytest.mark.parametrize('operation', ['stage', 'publish'])
def test_actual_clear_waits_for_inflight_stage_or_publish(tmp_path, operation):
    root = tmp_path / '.cache' / 'af3'
    destination = root / 'hash'
    stage = cache.begin_cache_stage(destination)
    write_structure(stage / 'model', 'model')
    ready = tmp_path / 'ready'
    release = tmp_path / 'release'
    worker_code = """
import sys, time
from pathlib import Path
from biofeaturefactory.alphafold3.bin import burst_manifest as cache
operation, root, source, ready, release = sys.argv[1:]
original = cache.tempfile.mkdtemp if operation == 'stage' else cache.is_cache_complete
def delayed(*args, **kwargs):
    Path(ready).touch()
    deadline = time.monotonic() + 10
    while not Path(release).exists():
        if time.monotonic() > deadline:
            raise RuntimeError('test release timeout')
        time.sleep(0.01)
    return original(*args, **kwargs)
if operation == 'stage':
    cache.tempfile.mkdtemp = delayed
    print(cache.begin_cache_stage(Path(root) / 'hash', '0'))
else:
    cache.is_cache_complete = delayed
    cache.publish_cache(Path(source), Path(root) / 'hash', 'model')
"""
    worker = subprocess.Popen([sys.executable, '-c', worker_code, operation, str(root), str(stage), str(ready), str(release)],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    clearer = None
    try:
        deadline = time.monotonic() + 10
        while not ready.exists() and worker.poll() is None and time.monotonic() < deadline:
            time.sleep(0.01)
        assert ready.exists()
        clearer = subprocess.Popen([sys.executable, '-c',
            'import sys; from pathlib import Path; from biofeaturefactory.alphafold3.bin.burst_manifest import clear_cache; clear_cache(Path(sys.argv[1]))', str(root)],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        with pytest.raises(subprocess.TimeoutExpired):
            clearer.communicate(timeout=0.3)
        release.touch()
        stdout, stderr = worker.communicate(timeout=10)
        assert worker.returncode == 0, stderr
        _, stderr = clearer.communicate(timeout=10)
        assert clearer.returncode == 0, stderr
        assert not destination.exists()
        if operation == 'stage':
            stage = Path(stdout.strip())
            assert stage.is_dir()
            write_structure(stage / 'model', 'model')
            with pytest.raises(ValueError, match='generation was cleared'):
                cache.publish_cache(stage, destination, 'model')
    finally:
        release.touch()
        for process in (worker, clearer):
            if process is not None and process.poll() is None:
                process.kill()
                process.communicate(timeout=10)
