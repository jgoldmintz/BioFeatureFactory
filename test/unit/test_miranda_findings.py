import argparse
import concurrent.futures
import json
from pathlib import Path
import subprocess
import shutil
import sys
import textwrap

import pandas as pd
import pytest

from biofeaturefactory.miranda import miranda_ensemble as miranda


def site_row(allele, position, score=140.0, energy=-20.4, join_pos=None, substrate="transcript"):
    return {
        "pkey": "TEST-key", "allele": allele, "mirna_id": "miR-test",
        "site_pos": position, "tot_score": score, "tot_energy": energy,
        "join_pos": position if join_pos is None else join_pos,
        "align_status": "aligned", "locus_id": "m1", "segment_id": "s1",
        "distance_to_snv": 100, "variant_kind": "insertion", "length_delta": 4,
        "substrate": substrate,
    }


def test_tied_representatives_are_order_independent():
    sites = pd.DataFrame([
        site_row("WT", 632), site_row("WT", 641, energy=-18.53),
        site_row("MUT", 641, energy=-18.53), site_row("MUT", 632),
    ])
    expected = miranda.build_events(sites)
    assert expected.iloc[0][["dpos", "delta_energy", "cls"]].tolist() == [0, 0, "none"]
    for seed in range(5):
        pd.testing.assert_frame_equal(expected, miranda.build_events(sites.sample(frac=1, random_state=seed)))


def test_tied_representatives_use_shared_indel_coordinates():
    sites = pd.DataFrame([
        site_row("WT", 100, join_pos=120), site_row("WT", 109, join_pos=129, energy=-18.53),
        site_row("MUT", 129, energy=-18.53), site_row("MUT", 120),
    ])
    event = miranda.build_events(sites).iloc[0]
    assert event[["wt_pos", "mut_pos", "dpos", "delta_energy", "cls"]].tolist() == [100, 120, 0, 0, "none"]


@pytest.mark.parametrize("mutant_rows,expected", [
    ([site_row("MUT", 641)], "shifted"),
    ([site_row("MUT", 632), site_row("MUT", 641, score=150)], "shifted"),
    ([site_row("MUT", 632, score=150)], "strengthened"),
    ([], "lost"),
])
def test_representative_selection_preserves_real_changes(mutant_rows, expected):
    sites = pd.DataFrame([site_row("WT", 632)] + mutant_rows)
    assert miranda.build_events(sites).iloc[0]["cls"] == expected


def test_representatives_do_not_cross_substrates():
    sites = pd.DataFrame([site_row("WT", 632), site_row("MUT", 632, substrate="intron1")])
    assert set(miranda.build_events(sites)["cls"]) == {"lost", "gained"}


def prediction_text(positions):
    return "".join(f">miR-test target 140 -20.4 _ 1 20 {position} {position + 20} 20 100 100\n" for position in positions)


def write_outputs(directory, gene_key, token, positions):
    (directory / f"{gene_key}-wt-miranda.out").write_text(prediction_text(positions))
    (directory / f"{miranda.mint_pkey(gene_key, token)}-mut-miranda.out").write_text(prediction_text(positions))


@pytest.mark.parametrize("substrate", ["transcript", "intron1", "pre_mRNA"])
@pytest.mark.parametrize("with_transcript", [False, True])
@pytest.mark.parametrize("token,statuses", [
    ("A10ATTT", ["aligned", "inserted", "inserted", "aligned"]),
    ("AA10A", ["aligned"] * 4),
    ("invalid", ["unprojected"] * 4),
])
def test_alignment_state_is_record_local(tmp_path, substrate, with_transcript, token, statuses):
    key = miranda._gene_key("TEST", substrate)
    mappings = {key: pd.DataFrame([("A10C", token)], columns=["mutant", "transcript"])}
    write_outputs(tmp_path, key, "A10C", [10, 11, 13, 14])
    if with_transcript and substrate != "transcript":
        mappings["TEST"] = pd.DataFrame([("A100G", "A100G")], columns=["mutant", "transcript"])
        write_outputs(tmp_path, "TEST", "A100G", [100])
    sites, rejected = miranda.build_sites_table_from_outputs(str(tmp_path), mappings, {})
    mutant = sites[(sites["substrate"] == substrate) & (sites["allele"] == "MUT")]
    assert mutant["align_status"].tolist() == statuses
    assert mutant["join_pos"].tolist() == [10, 11, 13, 14]
    assert bool(rejected) == (token == "invalid")


@pytest.mark.parametrize("substrate", ["intron1", "pre_mRNA"])
def test_whole_piece_deletion_projects_without_variant_object(tmp_path, substrate):
    key = miranda._gene_key("TEST", substrate)
    mappings = {key: pd.DataFrame([("AA10A", "AA10del")], columns=["mutant", "transcript"])}
    write_outputs(tmp_path, key, "AA10A", [10, 11, 13])
    sites, rejected = miranda.build_sites_table_from_outputs(str(tmp_path), mappings, {})
    assert not rejected
    assert sites.loc[sites["allele"] == "WT", "align_status"].tolist() == ["deleted", "deleted", "aligned"]
    assert sites.loc[sites["allele"] == "MUT", "align_status"].tolist() == ["aligned"] * 3


@pytest.fixture
def concurrent_predictor(tmp_path):
    binary_dir = tmp_path / "binary with spaces"
    binary_dir.mkdir()
    executable = binary_dir / "miranda"
    executable.write_text(f"#!{sys.executable}\n" + textwrap.dedent("""\
        import os
        from pathlib import Path
        import sys
        import time
        sequence = Path(sys.argv[2]).read_text().splitlines()[1]
        scratch = Path('predictor.tmp')
        scratch.write_text(sequence)
        barrier = Path(Path(sys.argv[1]).read_text())
        (barrier / sequence).write_text(str(Path.cwd()))
        deadline = time.monotonic() + 10
        while len(list(barrier.iterdir())) < 2 and time.monotonic() < deadline:
            time.sleep(0.01)
        time.sleep(0.05)
        print(sequence, scratch.read_text(), Path(sys.argv[2]).read_text().splitlines()[1])
        scratch.unlink()
    """))
    executable.chmod(0o755)
    barrier = tmp_path / "barrier"
    barrier.mkdir()
    database = tmp_path / "mirna.fa"
    database.write_text(str(barrier))
    return binary_dir, barrier, database


@pytest.mark.parametrize("phase", ["WT", "MUT"])
def test_concurrent_predictors_have_private_working_directories(tmp_path, concurrent_predictor, phase):
    binary_dir, barrier, database = concurrent_predictor
    def invoke(sequence):
        output = tmp_path / sequence
        output.mkdir()
        if phase == "WT":
            miranda.run_wt_phase({"TEST": sequence}, str(output), str(binary_dir), str(database))
            return (output / "TEST-wt-miranda.out").read_text().strip()
        result = miranda._run_single_miranda_task(("TEST-key-mut", sequence, str(binary_dir), str(output), str(database), False))
        assert result == ("TEST-key-mut", "ok")
        return (output / "TEST-key-mut-miranda.out").read_text().strip()
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(invoke, ["AAAA", "CCCC"]))
    assert results == ["AAAA AAAA AAAA", "CCCC CCCC CCCC"]
    assert len({entry.read_text() for entry in barrier.iterdir()}) == 2
    assert {entry.name for entry in binary_dir.iterdir()} == {"miranda"}


def test_concurrent_main_invocations_use_private_intermediates(tmp_path):
    runner = textwrap.dedent("""\
        import argparse
        import json
        from pathlib import Path
        import sys
        import time
        import pandas as pd
        from biofeaturefactory.miranda import miranda_ensemble as pipeline
        output, barrier, fasta, mapping, identity = map(Path, sys.argv[1:])
        args = argparse.Namespace(input=str(fasta), output=str(output), miranda_dir=None,
            mirna_db='unused', log=None, mapping_dir=str(mapping), wt_header='transcript',
            intron_premrna_mapping=None, strict_introns=False, no_parallel=True, max_workers=1)
        argparse.ArgumentParser.parse_args = lambda self: args
        import shutil
        shutil.which = lambda binary: '/stub/miranda'
        pipeline.load_transcript_mappings = lambda *args: {'TEST': pd.DataFrame([('A10C', 'A10C')], columns=['mutant', 'transcript'])}
        pipeline.load_wt_sequences = lambda *args, **kwargs: {'TEST': 'A' * 40}
        pipeline._load_substrate_sequences = lambda *args: {}
        pipeline.discover_mapping_files = lambda *args: {}
        def predictor(args):
            seq_id, sequence, binary, workdir, database, strict = args
            workdir = Path(workdir)
            if seq_id.endswith('-wt'):
                (barrier / identity).write_text(str(workdir))
                deadline = time.monotonic() + 10
                while len(list(barrier.iterdir())) < 2 and time.monotonic() < deadline:
                    time.sleep(0.01)
            (workdir / f'{seq_id}-miranda.out').write_text('>miR-test target 140 -20.4 _ 1 20 20 40 20 100 100\\n')
            return seq_id, 'ok'
        pipeline._run_single_miranda_task = predictor
        def wt_phase(wt_sequences, outdir, miranda_dir, mirna_db, strict_substrates):
            predictor(('TEST-wt', 'A' * 40, None, outdir, None, False))
        pipeline.run_wt_phase = wt_phase
        pipeline.main()
    """)
    output = tmp_path / "output"
    output.mkdir()
    barrier = tmp_path / "barrier"
    barrier.mkdir()
    fasta = tmp_path / "TEST.fasta"
    fasta.write_text(">transcript\n" + "A" * 40 + "\n")
    mapping = tmp_path / "transcript_mapping_TEST.csv"
    mapping.write_text("mutant,transcript\nA10C,A10C\n")
    command = [sys.executable, "-c", runner, str(output), str(barrier), str(fasta), str(mapping)]
    processes = [subprocess.Popen(command + [identity], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) for identity in ("first", "second")]
    for process in processes:
        stdout, stderr = process.communicate(timeout=30)
        assert process.returncode == 0, stdout + stderr
    workdirs = {entry.read_text() for entry in barrier.iterdir()}
    assert len(workdirs) == 2
    assert all(Path(directory).parent == output and Path(directory).is_dir() for directory in workdirs)
    summary = json.loads((output / "miranda.run_summary.json").read_text())
    assert summary["pkeys_with_rows"] == 1
    assert summary["parse_rejected"] == 0


def test_failed_invocation_retains_unparsed_outputs(tmp_path, monkeypatch):
    fasta = tmp_path / "TEST.fasta"
    fasta.write_text(">transcript\n" + "A" * 40 + "\n")
    mapping = tmp_path / "transcript_mapping_TEST.csv"
    mapping.write_text("mutant,transcript\nA10C,A10C\n")
    args = argparse.Namespace(
        input=str(fasta), output=str(tmp_path), miranda_dir=None, mirna_db="unused",
        log=None, mapping_dir=str(mapping), wt_header="transcript",
        intron_premrna_mapping=None, strict_introns=False, no_parallel=True, max_workers=1,
    )
    monkeypatch.setattr(argparse.ArgumentParser, "parse_args", lambda self: args)
    monkeypatch.setattr(shutil, "which", lambda binary: "/stub/miranda")
    monkeypatch.setattr(miranda, "load_transcript_mappings", lambda *args: {
        "TEST": pd.DataFrame([("A10C", "A10C")], columns=["mutant", "transcript"]),
    })
    monkeypatch.setattr(miranda, "load_wt_sequences", lambda *args, **kwargs: {"TEST": "A" * 40})
    monkeypatch.setattr(miranda, "_load_substrate_sequences", lambda *args: {})
    monkeypatch.setattr(miranda, "discover_mapping_files", lambda *args: {})

    def fail_after_output(**kwargs):
        (Path(kwargs["outdir"]) / "TEST-wt-miranda.out").write_text("retained diagnostic output")
        raise RuntimeError("predictor failed")

    monkeypatch.setattr(miranda, "run_wt_phase", fail_after_output)
    with pytest.raises(RuntimeError, match="predictor failed"):
        miranda.main()
    retained = list(tmp_path.glob(".miranda-*/TEST-wt-miranda.out"))
    assert len(retained) == 1
    assert retained[0].read_text() == "retained diagnostic output"
