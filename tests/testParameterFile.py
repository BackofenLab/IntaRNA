"""Compare JSON/embedded inputs with actual native CLI/config-file predictions."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

binary = str(Path(sys.argv[1]).resolve())
with tempfile.TemporaryDirectory(prefix="intarna-parameters-test-") as temporary:
    root = Path(temporary)
    (root / "tmp").mkdir()
    env = dict(os.environ, TMPDIR=str(root / "tmp"))

    def run(*args, success=True):
        result = subprocess.run([binary, *args], cwd=root, env=env, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        assert (result.returncode == 0) == success, result.stdout
        assert not list((root / "tmp").iterdir()), "Embedded inputs leaked after invocation"
        return result.stdout

    def write(parameters, name="parameters.json"):
        (root / name).write_text(json.dumps({"schema": "IntaRNA-parameters-v1", "parameters": parameters}))
        return f"--parameterFile={name}"

    common = {"outMode": "C", "outCsvCols": "id1,start1,end1,id2,start2,end2,E,ED1,ED2", "out": "result.csv"}
    q = ">query\nGGGAGGAGGG\n"
    t = ">target\nCCCUCCUCCC\n"
    (root / "query.fa").write_text(q)
    (root / "target.fa").write_text(t)
    direct = {**common, "query": "query.fa", "target": "target.fa"}
    run(*(f"--{k}={v}" for k, v in direct.items()))
    baseline = (root / "result.csv").read_text()
    embedded = {**common, "query": {"text": q}, "target": {"text": t}}
    run(write(embedded))
    assert (root / "result.csv").read_text() == baseline
    (root / "legacy.cfg").write_text("\n".join(f"{k} = {v}" for k, v in direct.items()))
    run("--parameterFile=legacy.cfg")
    assert (root / "result.csv").read_text() == baseline

    # Explicit scalar overrides, signed indices and file format detection/BOM.
    run(write({**embedded, "temperature": 15}, "settings.txt"), "--temperature=37")
    assert (root / "result.csv").read_text() == baseline
    path = root / "settings.txt"
    path.write_bytes(b"\xef\xbb\xbf  \n" + path.read_bytes())
    run("--parameterFile", "settings.txt", "--temperature=37")
    assert (root / "result.csv").read_text() == baseline
    indexed = {**direct, "qIdxPos0": -4, "tIdxPos0": 20}
    run(*(f"--{k}={v}" for k, v in indexed.items()))
    expected = (root / "result.csv").read_text()
    run(write({**embedded, "qIdxPos0": -4, "tIdxPos0": 20}))
    assert (root / "result.csv").read_text() == expected

    # Multi-FASTA and repeated outputs preserve native file naming/results.
    (root / "query.fa").write_text(q + ">second\nGGGAGGAGGG\n")
    run(*(f"--{k}={v}" for k, v in direct.items()), "--out=qPu:pu.txt")
    expected = (root / "result.csv").read_text()
    pu = [(root / f"pu-s{i}.txt").read_text() for i in (1, 2)]
    multi = {**embedded, "query": {"text": (root / "query.fa").read_text()},
             "out": ["result.csv", "qPu:embedded-pu.txt"]}
    run(write(multi))
    assert (root / "result.csv").read_text() == expected
    assert (root / "embedded-pu-s2.txt").read_text() == pu[1]
    run(write({**multi, "qAcc": "P", "qAccFile": {"texts": pu}}))
    imported = (root / "result.csv").read_text()
    run(*(f"--{k}={v}" for k, v in direct.items()), "--qAcc=P", "--qAccFile=pu.txt")
    assert (root / "result.csv").read_text() == imported
    (root / "query.fa").write_text(q)

    # RNAplfold opening-energy and probability text, and SHAPE input.
    run(*(f"--{k}={v}" for k, v in direct.items()), "--out=qAcc:ed.txt", "--out=tPu:pu.txt")
    files = {"qAccFile": "ed.txt", "tAccFile": "pu.txt", "qShape": "shape.txt", "tShape": "shape.txt"}
    (root / "shape.txt").write_text("\n".join(f"{i} 0.2" for i in range(1, 11)) + "\n")
    for settings in [{"qAcc": "E", "qAccFile": "ed.txt"}, {"tAcc": "P", "tAccFile": "pu.txt"}, {"qShape": "shape.txt", "tShape": "shape.txt"}]:
        run(*(f"--{k}={v}" for k, v in {**direct, **settings}.items()))
        expected = (root / "result.csv").read_text()
        encoded = {k: {"text": (root / v).read_text()} if k in files else v for k, v in settings.items()}
        run(write({**embedded, **encoded}))
        assert (root / "result.csv").read_text() == expected

    # File personality must select defaults before constructing the CLI parser.
    run(*(f"--{k}={v}" for k, v in direct.items()), "--personality=IntaRNAduplex")
    expected = (root / "result.csv").read_text()
    run(write({**embedded, "personality": "IntaRNAduplex"}))
    assert (root / "result.csv").read_text() == expected
    run(write({**embedded, "personality": "IntaRNAduplex"}), "--personality=IntaRNA")
    assert (root / "result.csv").read_text() == baseline
    evaluation = {"personality": "IntaRNAeval", "rri": "1||||||||||&1||||||||||"}
    run(*(f"--{k}={v}" for k, v in {**direct, **evaluation}.items()))
    expected = (root / "result.csv").read_text()
    run(write({**embedded, **evaluation}))
    assert (root / "result.csv").read_text() == expected

    # Report replay can clear CSV-only settings from the same JSON document.
    sorted_settings = {**embedded, "outCsvSort": "E"}
    run(write(sorted_settings), "--outMode=D", "--outCsvCols=id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E", "--outCsvSort", "")
    assert "-21.4" in (root / "result.csv").read_text()

    # Reject malformed/schema/type errors and ensure cleanup even after errors.
    for bad in [{"unknown": 1}, {"temperature": None}, {"temperature": []}, {"temperature": {}},
                {"temperature": 999}, {"query": {"text": "bad sequence!"}}, {"parameterFile": "other.json"},
                {"out": []}, {"query": {"text": "a\u0000b"}}, {"help": "true"}, {"personality": "unknown"}]:
        run(write({**embedded, **bad}), success=False)
    for text in ['[]', '{"schema":"unknown","parameters":{}}', '{broken',
                 '{"schema":"IntaRNA-parameters-v1","parameters":{"query":"AAA","query":"CCC"}}']:
        (root / "bad.json").write_text(text)
        run("--parameterFile=bad.json", success=False)
print("JSON/legacy parameters, embedded inputs, personality, precedence and errors passed")
