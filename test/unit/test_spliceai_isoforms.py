import os
from pathlib import Path
import shutil
import subprocess
import sys
import textwrap

import pytest


PIPELINE = Path(__file__).resolve().parents[2] / "biofeaturefactory/spliceai/bin/main.nf"


@pytest.mark.skipif(
    os.environ.get("BFF_TEST_NEXTFLOW") != "1" or shutil.which("nextflow") is None,
    reason="Set BFF_TEST_NEXTFLOW=1 with Nextflow installed to test the live process",
)
@pytest.mark.parametrize("missing_annotation", [False, True])
def test_live_spliceai_process_counts_exact_gene_rows(tmp_path, missing_annotation):
    source = PIPELINE.read_text()
    process = "process run_spliceai {" + source.split("process run_spliceai {", 1)[1].split("\n  process parse_results", 1)[0]
    workflow = tmp_path / "main.nf"
    workflow.write_text(textwrap.dedent("""\
        nextflow.enable.dsl=2
        params.maxforks=2
        params.spliceai_accelerator=1
        params.output_dir='results'
        params.vcf_output_dir=null
        params.spliceai_gpus=1
        params.spliceai_mem_growth=true
        params.spliceai_env=''
        params.annotation='annotation.txt'
    """) + process + textwrap.dedent("""\
        workflow {
            inputs = Channel.of('NONE', 'ONE', 'MANY', 'REG.X').map { gene ->
                tuple(gene, file('input.vcf.gz'), file('input.vcf.gz.tbi'))
            }
            run_spliceai(inputs, 'reference.fa', params.annotation, false, 2)
        }
    """))
    annotation = tmp_path / "annotation.txt"
    if not missing_annotation:
        annotation.write_text("#NAME\tCHROM\nONE\t1\nMANY\t1\nMANY\t2\nMANY\t3\nREG.X\t1\nREGAX\t1\n")
    (tmp_path / "input.vcf.gz").write_bytes(b"stub")
    (tmp_path / "input.vcf.gz.tbi").write_bytes(b"stub")
    binaries = tmp_path / "stub-bin"
    binaries.mkdir()
    python = binaries / "python3"
    python.write_text(f"#!{sys.executable}\n" + textwrap.dedent("""\
        from pathlib import Path
        import sys
        args = sys.argv[1:]
        if args[0] == '-c':
            Path(args[args.index('-O') + 1]).write_text('stub prediction\\n')
            print('PREDICTION_ANNOTATION=' + args[args.index('-A') + 1])
        else:
            Path(args[args.index('--output') + 1]).write_text('stub filtered annotation\\n')
            print('FILTER_CALLED')
    """))
    python.chmod(0o755)
    environment = dict(os.environ, NXF_OFFLINE="true", NXF_ANSI_LOG="false")
    environment["PATH"] = str(binaries) + os.pathsep + environment["PATH"]
    result = subprocess.run(
        ["nextflow", "run", str(workflow), "--annotation", str(annotation)],
        cwd=tmp_path, env=environment, capture_output=True, text=True, timeout=120,
    )
    outputs = [path.read_text() for path in (tmp_path / "work").glob("*/*/.command.out")]
    all_output = "\n".join(outputs)
    if missing_annotation:
        assert result.returncode != 0, result.stdout + result.stderr
        assert "isoforms detected" not in all_output
        assert "PREDICTION_ANNOTATION" not in all_output
        errors = "\n".join(path.read_text() for path in (tmp_path / "work").glob("*/*/.command.err"))
        assert "annotation.txt" in errors
    else:
        assert result.returncode == 0, result.stdout + result.stderr
        assert len(outputs) == 4
        for gene, count in (("NONE", 0), ("ONE", 1), ("MANY", 3), ("REG.X", 1)):
            assert f"[run_spliceai] {gene}: {count} isoforms detected" in all_output
        assert all_output.count("FILTER_CALLED") == 1
        assert "PREDICTION_ANNOTATION=MANY_filtered_annotation.txt" in all_output
