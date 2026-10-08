#!/usr/bin/env python3
"""Consumer checks for the motif workflow and its VCF delivery functions."""
import argparse
import copy
import json
import random
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))


def fails(call, text=None):
    try:
        call()
    except (ValueError, RuntimeError, OSError) as error:
        if text is not None:
            assert text in str(error), str(error)
    else:
        raise AssertionError("Expected consumer-visible failure")


def fixture(destination, mode, image_cache):
    destination.mkdir(parents=True, exist_ok=True)
    rng = random.Random(187)
    contigs = ["segment_A", "7", "scaffold.with-dots", "unit-Z", "pieceQ", "extra_sequence"]
    patterns = ["ACGTTGCACTGA", "GCTAGACCTTGC"]
    sequences = {name: [rng.choice("ACGT") for _ in range(40000)] for name in contigs}
    positions = list(range(201, 24201, 200))
    for name in contigs[:5]:
        for number, pos in enumerate(positions):
            motif = patterns[number % 2]
            sequences[name][pos - 1:pos - 1 + len(motif)] = motif
    reference = destination / "reference.fa"
    with reference.open("w", encoding="utf-8", newline="\n") as stream:
        for name, sequence in sequences.items():
            stream.write(">" + name + "\n")
            text = "".join(sequence)
            for start in range(0, len(text), 80):
                stream.write(text[start:start + 80] + "\n")
    pb = ["PB-library-A", "PB-library-B"]
    il_a = [f"IL-library-A-{i}" for i in range(1, 7)]
    il_b = [f"IL-library-B-{i}" for i in range(1, 6)]
    external = [f"OTHER-library-{i}" for i in range(1, 7)]

    def vcf(path, samples, external=False):
        with path.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write("##fileformat=VCFv4.2\n")
            for name in contigs:
                stream.write(f"##contig=<ID={name},length=40000>\n")
            stream.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            stream.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            for name in contigs:
                selected = [25, *positions, 39975] if name != contigs[-1] else [25, 1001, 39975]
                for index, pos in enumerate(selected):
                    ref = sequences[name][pos - 1]
                    alt = next(x for x in "ACGT" if x != ref)
                    calls = ["0/1"] * len(samples)
                    if samples[-1].startswith("UNLISTED"):
                        calls[-1] = "0/0"
                    if index == 2:
                        calls[0] = "1/."
                        if len(calls) > 1:
                            calls[1] = "0/0"
                    if index == 3 and len(calls) > 1:
                        calls[0], calls[1] = ".", "0/1"
                    if external and index == 4:
                        calls[0] = "./."
                    if name == contigs[-1] and not external:
                        calls = ["0/0"] * len(samples)
                        calls[-1] = "0/1"
                    row = [name, str(pos), f"record-{index}", ref, alt, ".", "PASS", ".", "GT", *calls]
                    stream.write("\t".join(row) + "\n")
                    if index == 5:
                        stream.write("\t".join(row) + "\n")
                    if index == 6:
                        row[4] = alt + "," + next(x for x in "ACGT" if x not in (ref, alt))
                        row[9] = "1/2"
                        stream.write("\t".join(row) + "\n")
    pb_path, il_path, external_path = [destination / name for name in ("pacbio.vcf", "corresponding.vcf", "external.vcf")]
    vcf(pb_path, pb + ["UNLISTED-PB"])
    vcf(il_path, il_a + il_b + ["UNLISTED-IL"])
    vcf(external_path, external, True)
    import yaml
    config = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    config.update(reference_fasta=str(reference), pacbio_vcf=str(pb_path), illumina_vcf=str(il_path),
                  other_illumina=str(external_path), output_dir=str(destination / "results"),
                  scratch_dir=str(destination / "scratch"), pacbio_samples=pb,
                  illumina_corresponding={pb[0]: il_a, pb[1]: il_b})
    config["folds"]["grouping"] = {"enabled": False, "regex": None, "unmatched_policy": "error"}
    config["vcf_outputs"] = {"consensus_fraction": 0.75, "groups": {
        "corresponding_a": {"role": "illumina", "samples": il_a},
        "corresponding_b": {"role": "illumina", "samples": il_b},
        "external": {"role": "other_illumina", "samples": external}}}
    config["execution"].update(use_containers=mode == "containers", image_cache=str(image_cache or destination / "images"))
    config["sample_metadata"] = {
        "pacbio": {p: {"clone_id": f"group-{i}"} for i, p in enumerate(pb)},
        "illumina": {name: {"clone_id": f"group-{i}"} for i, names in enumerate([il_a, il_b]) for name in names},
        "other_illumina": {name: {"clone_id": f"external-{i}"} for i, name in enumerate(external)},
    }
    path = destination / "config.yaml"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return config, path


def check_configuration():
    import yaml
    import prepare
    config = yaml.safe_load((ROOT / "config.yaml").read_text(encoding="utf-8"))
    effective = prepare.validate_config(copy.deepcopy(config), ROOT)
    assert effective["vcf_outputs"] == {"consensus_fraction": 0.75, "groups": {}}
    assert "vcf_outputs" not in prepare.selection_config(effective)
    assert effective["analysis"]["rotations"] == [1, 2, 3, 4, 5]
    assert prepare.dataset_id("pacbio", "same") != prepare.dataset_id("illumina", "same")
    for path, value in [
        (("pacbio_samples",), []),
        (("analysis", "rotations"), [1, 2]),
        (("analysis", "flank_bp"), True),
        (("analysis", "center_radius_bp"), 51),
        (("analysis", "fimo_p_threshold"), 0),
        (("selection", "min_supported_rotations"), 6),
        (("selection", "corresponding_dataset_support_fraction"), 0),
        (("folds", "grouping", "enabled"), "yes"),
        (("folds", "grouping", "unmatched_policy"), "ignore"),
        (("fold_beds",), ["supplied.bed"]),
    ]:
        invalid = copy.deepcopy(config)
        node = invalid
        for key in path[:-1]:
            node = node[key]
        node[path[-1]] = value
        fails(lambda invalid=invalid: prepare.validate_config(invalid, ROOT))
    for minimum in (2, 3, 4):
        changed = copy.deepcopy(config)
        changed["selection"]["min_supported_rotations"] = minimum
        assert prepare.validate_config(changed, ROOT)["selection"]["min_supported_rotations"] == minimum
    valid = copy.deepcopy(effective)
    valid["vcf_outputs"]["groups"] = {
        "group_a": {"role": "illumina", "samples": ["IL-library-A-1", "IL-library-A-2"]},
        "group_b": {"role": "other_illumina", "samples": ["header-name-not-checked-yet"]},
    }
    assert prepare.validate_config(valid, ROOT)["vcf_outputs"] == valid["vcf_outputs"]
    for outputs in (
        {"consensus_fraction": True, "groups": {}},
        {"consensus_fraction": 0, "groups": {}},
        {"consensus_fraction": float("nan"), "groups": {}},
        {"consensus_fraction": 0.75, "groups": {"../bad": {"role": "illumina", "samples": ["IL-library-A-1"]}}},
        {"consensus_fraction": 0.75, "groups": {"empty": {"role": "illumina", "samples": []}}},
        {"consensus_fraction": 0.75, "groups": {"unknown": {"role": "illumina", "samples": ["missing"]}}},
        {"consensus_fraction": 0.75, "groups": {"duplicate": {"role": "other_illumina", "samples": ["x", "x"]}}},
    ):
        invalid = copy.deepcopy(effective)
        invalid["vcf_outputs"] = outputs
        fails(lambda invalid=invalid: prepare.validate_config(invalid, ROOT))
    print("PASS: VCF export defaults, role-scoped groups and frozen-selection isolation")
def check_application_binds():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "motif_profile_workflow", ROOT / "profile" / "workflow.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    config = {
        "output_dir": "/tmp/motif-output",
        "scratch_dir": "/tmp/motif-scratch",
        "reference_fasta": "/source/reference.fa",
        "pacbio_vcf": "/source/pacbio.vcf.gz",
        "illumina_vcf": "/source/corresponding.vcf.gz",
        "other_illumina": "/source/external.vcf.gz",
        "execution": {"image_cache": "/tmp/motif-images", "extra_bind_paths": []},
    }
    for keys, included, excluded in (
        ({"other_illumina"}, ("other_illumina",), ("illumina_vcf", "pacbio_vcf")),
        ({"illumina_vcf"}, ("illumina_vcf",), ("other_illumina", "pacbio_vcf")),
        ({"illumina_vcf", "other_illumina"}, ("illumina_vcf", "other_illumina"), ("pacbio_vcf",)),
    ):
        mounts = module.binds(config, "application", source_keys=keys)
        for key in included:
            assert str(Path(config[key]).absolute()) in mounts
        for key in excluded:
            assert str(Path(config[key]).absolute()) not in mounts
    fails(lambda: module.binds(config, "application"), "requires target-selected source keys")
    fails(lambda: module.binds(config, "application", source_keys={"pacbio_vcf"}), "Unsupported")
    print("PASS: application bind sets exclude unrelated and PacBio sources")




def check_preparation(destination):
    import prepare
    grouped = destination / "fold-fixture"
    grouped.mkdir()
    fai = grouped / "reference.fa.fai"
    contigs = [(f"contig{i}", 100 + 10 * i) for i in range(5)]
    fai.write_text("".join(f"{name}\t{length}\t0\t{length}\t{length + 1}\n" for name, length in contigs), encoding="utf-8")
    manifest = prepare.build_folds(fai, "reference-hash", grouped / "output", {"enabled": False, "regex": None, "unmatched_policy": "error"})
    assert len(manifest["folds"]) == 5 and len(manifest["assignments"]) == 5
    assert all((grouped / "output" / f"folds/fold{i}.bed").is_file() for i in range(1, 6))
    assert prepare.read_fai(fai) == dict(contigs)
    for gt in ("1/.", "0/1", "1|2/.", "2"):
        prepare.parse_gt(gt, 2)
    for gt in ("3/0", "-1", "0//1", "1x2", ""):
        fails(lambda gt=gt: prepare.parse_gt(gt, 2))
    print("PASS: deterministic five-fold preparation and GT validation")


def _frozen_fixture(out):
    import prepare
    import annotate
    library = out / "library"
    reports = out / "reports"
    library.mkdir(parents=True)
    reports.mkdir()
    reference = out / "reference.fa"
    reference.write_text(">chr1\nACGT\n", encoding="utf-8")
    Path(str(reference) + ".fai").write_text("chr1\t4\t6\t4\t5\n", encoding="utf-8")
    artifacts = {
        "final_motifs.meme": "",
        "final_motif_manifest.tsv": "final_motif_id\tfamily_id\tglobal_rank\twidth\n",
        "final_motif_members.tsv": "final_motif_id\tfamily_id\tmember_id\n",
        "scoring_background.bfile": "background\n",
    }
    for name, text in artifacts.items():
        (library / name).write_text(text, encoding="utf-8")
    loci = sqlite3.connect(library / "selection_loci.sqlite")
    loci.execute("CREATE TABLE loci (CHROM TEXT, POS INTEGER, PRIMARY KEY(CHROM, POS))")
    loci.commit()
    loci.close()
    validation = reports / "representative_validation.tsv"
    validation.write_text("candidate_pass\nTrue\n", encoding="utf-8")
    paths = [library / name for name in artifacts] + [library / "selection_loci.sqlite", validation]
    freeze = {
        "schema_version": 1,
        "outputs": {str(path.resolve()): prepare.sha256(path) for path in paths},
        "representative_validation": {"sha256": prepare.sha256(validation)},
        "family_stability": {"rotations": [1, 2, 3, 4, 5]},
        "reference_fasta": str(reference.resolve()),
        "reference_sha256": prepare.sha256(reference),
        "flank_bp": 2,
        "fimo_p_threshold": 0.0001,
        "motif_count": 0,
        "background_sha256": prepare.sha256(library / "scoring_background.bfile"),
        "selection_config": {},
    }
    prepare.atomic_json(library / "freeze.json", freeze)
    return annotate


def check_annotation_roles(destination):
    import annotate
    import prepare
    out = destination / "annotation-role-fixture"
    _frozen_fixture(out)
    for role, name, participant, context, identity in (
        ("other_illumina", "external-A", "false", "external_application_identity_unverified", "header_names_only"),
        ("illumina", "selected-A", "true", "selection_reuse", "frozen_selection_manifest"),
    ):
        did = prepare.dataset_id(role, name)
        app, _ = annotate.annotation_paths(out, role)
        (app / "shards").mkdir(parents=True, exist_ok=True)
        prepare.atomic_json(app / "manifest.json", {
            "identity_check_status": identity, "flank_bp": 2,
            "datasets": [{"dataset_id": did, "role": role, "sample_name": name,
                          "selection_participant": participant, "evidence_context": context,
                          "identity_check_status": identity}],
        })
        prepare.atomic_json(app / "shards/manifest.json", {"shards": []})
        db = sqlite3.connect(app / "shards/windows.sqlite")
        db.execute("CREATE TABLE windows (window_id TEXT PRIMARY KEY, sequence TEXT NOT NULL)")
        db.commit()
        db.close()
        directory = out / "prepared" / did
        directory.mkdir(parents=True, exist_ok=True)
        row = {
            "source_vcf": "source.vcf", "source_record_index": 1, "record_id": "record-1",
            "dataset_id": did, "sample_name": name, "CHROM": "chr1", "POS": 1, "ID": "r1",
            "REF": "A", "ALT": "C", "GT": "0/0", "variant_type": "SNP", "partial_missing": 0,
            "ambiguity_count": 0, "window_id": "", "interval_start": "", "interval_end": "",
            "extraction_status": "eligible", "extraction_reason": "",
        }
        prepare.write_tsv(directory / "genotypes.tsv.gz", [row], list(row))
        annotate.annotations(out, [], role=role)
        rows = list(prepare.read_tsv(app / "genotypes.tsv.gz"))
        assert len(rows) == 1
        assert (rows[0]["role"], rows[0]["selection_participant"], rows[0]["evidence_context"]) == (role, participant, context)
        assert annotate.load_json(app / "annotation.complete.json")["role"] == role
    assert (out / "application/genotypes.tsv.gz").is_file()
    assert (out / "application/corresponding/genotypes.tsv.gz").is_file()
    print("PASS: shared frozen-library annotation keeps external and corresponding namespaces separate")


def check_vcf_outputs(destination):
    import annotate
    import prepare
    import vcf_outputs as exports

    out = destination / "vcf-output-fixture"
    _frozen_fixture(out)
    reference = out / "reference.fa"
    reference.write_text(">chr1\nCAAAAG\n", encoding="utf-8")
    Path(str(reference) + ".fai").write_text("chr1\t6\t6\t6\t7\n", encoding="utf-8")
    freeze_path = out / "library/freeze.json"
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    freeze["reference_sha256"] = prepare.sha256(reference)
    selected, external = "IL-test", "OTHER-test"
    freeze["selection_config"] = {"illumina_corresponding": {"PB-test": [selected]}}
    prepare.atomic_json(freeze_path, freeze)
    sources = {
        "illumina": out / "corresponding.vcf",
        "other_illumina": out / "external.vcf",
    }
    rows = [
        ["chr1", "1", "del-left", "CA", "C", "41", "PASS", "AC=1;AN=2", "DP:GT", "8:0/1"],
        ["chr1", "2", "del-right", "AA", "A", "42", "PASS", "AC=2;AN=2", "GT:DP", "1|1:12"],
        ["chr1", "4", "mixed", "A", "G,<DEL>", "43", "PASS", "AC=1,1;AN=2", "GT:DP", "1/2:7"],
        ["chr1", "5", "breakend", "A", "A[chr1:2[", "44", "PASS", "AC=1;AN=2", "GT:DP", "0/1:5"],
        ["chr1", "5", "spanning", "A", "*", "45", "PASS", "AC=1;AN=2", "GT", "0/1"],
        ["chr1", "5", "unspecified", "A", "<NON_REF>", "46", "PASS", "AC=1;AN=2", "GT", "0/1"],
        ["chr1", "5", "ambiguous", "A", "N", "47", "PASS", "AC=1;AN=2", "GT", "0/1"],
        ["chr1", "6", "reference", "G", ".", "48", "PASS", "AC=0;AN=2", "GT", "0/0"],
    ]

    def write_source(path, sample):
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("##contig=<ID=chr1,length=6>\n")
            handle.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n')
            handle.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">\n')
            handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")
            for row in rows:
                handle.write("\t".join(row) + "\n")

    write_source(sources["illumina"], selected)
    write_source(sources["other_illumina"], external)
    freeze_hash = prepare.sha256(freeze_path)
    for role, name in (("illumina", selected), ("other_illumina", external)):
        source = sources[role]
        did = prepare.dataset_id(role, name)
        source_hash = prepare.sha256(source)
        identity = "frozen_selection_manifest" if role == "illumina" else "header_names_only"
        participant = role == "illumina"
        context = "selection_reuse" if participant else "external_application_identity_unverified"
        app, _ = annotate.annotation_paths(out, role)
        app.mkdir(parents=True, exist_ok=True)
        annotation = app / "genotypes.tsv.gz"
        records = []
        for ordinal, source_row in enumerate(rows, 1):
            keys = source_row[8].split(":")
            values = source_row[9].split(":")
            gt = values[keys.index("GT")] if "GT" in keys else "."
            row = {key: "" for key in annotate.FIELDS}
            row.update(source_vcf=str(source), source_record_index=str(ordinal),
                       record_id=f"{role}_{source_hash}_{ordinal}", dataset_id=did,
                       sample_name=name, role=role,
                       selection_participant=str(participant).lower(), evidence_context=context,
                       CHROM=source_row[0], POS=source_row[1], ID=source_row[2], REF=source_row[3],
                       ALT=source_row[4], GT=gt, annotation_status="promising" if participant and ordinal < 8 else "unclassified",
                       identity_check_status=identity, extraction_status="eligible",
                       interval_start="", interval_end="", window_id="")
            records.append(row)
        prepare.write_tsv(annotation, records, annotate.FIELDS)
        manifest_row = {"dataset_id": did, "role": role, "sample_name": name,
                        "source_vcf": str(source), "corresponding_pacbio": "PB-test" if participant else "",
                        "identity_check_status": identity, "selection_participant": str(participant).lower(),
                        "evidence_context": context}
        application = {
            "schema_version": 1, "datasets": [manifest_row], "source_vcf": str(source),
            "source_sha256": source_hash, "identity_check_status": identity,
            "reference": str(reference), "flank_bp": 2,
        }
        prepare.atomic_json(app / "manifest.json", application)
        complete = app / "annotation.complete.json"
        prepare.atomic_json(complete, {"schema_version": 1, "role": role, "freeze_sha256": freeze_hash})
        prepare.publish_provenance(complete, [source, app / "manifest.json", annotation],
                                   {"fixture": True}, [complete, annotation], scope="application")
    config = {
        "output_dir": str(out), "illumina_vcf": str(sources["illumina"]),
        "other_illumina": str(sources["other_illumina"]),
        "illumina_corresponding": {"PB-test": [selected]},
        "vcf_outputs": {"consensus_fraction": 0.75, "groups": {
            "selected": {"role": "illumina", "samples": [selected]},
            "external": {"role": "other_illumina", "samples": [external]},
        }},
    }
    datasets, groups = exports.manifest(config, out)
    assert [row["sample_name"] for row in datasets] == [selected, external]
    assert [group["required_support"] for group in groups] == [1, 1]
    assert json.loads((out / "vcf_outputs/groups.json").read_text(encoding="utf-8"))["consensus_fraction"] == "3/4"
    selected_id = prepare.dataset_id("illumina", selected)
    external_id = prepare.dataset_id("other_illumina", external)
    exports.sample_inputs(config, out, selected_id)
    exports.sample_inputs(config, out, external_id)
    selected_dir = out / "vcf_outputs" / selected_id
    external_dir = out / "vcf_outputs" / external_id
    selected_data = [line.rstrip("\n") for line in (selected_dir / "promising.raw.vcf").read_text().splitlines()
                     if line and not line.startswith("#")]
    expected_data = ["\t".join(source_row[:9] + [source_row[9]]) for source_row in rows[:7]]
    assert selected_data == expected_data
    raw_header = (selected_dir / "promising.raw.vcf").read_text()
    for metadata in ("##motif_label=prioritization_only", "##motif_selection_participant=true",
                     "##motif_evidence_context=selection_reuse", "##motif_freeze_sha256="):
        assert metadata in raw_header
    assert not [line for line in (external_dir / "promising.raw.vcf").read_text().splitlines()
                if line and not line.startswith("#")]
    source_map = list(prepare.read_tsv(selected_dir / "source_alleles.tsv.gz"))
    assert len(source_map) == 8 and all(row["original_alt_index"] for row in source_map)
    assert json.loads(source_map[2]["original_record_fields"])[4] == "G,<DEL>"
    normalized = selected_dir / "normalized_alleles.vcf"
    subprocess.run(["bcftools", "norm", "--no-version", "-f", str(reference), "-c", "e", "-Ov",
                    "-o", str(normalized), str(selected_dir / "source_alleles.vcf")], check=True)
    exports.capture_alleles(config, out, selected_id)
    alleles = list(prepare.read_tsv(selected_dir / "alleles.tsv.gz"))
    deletion_keys = {(row["canonical_CHROM"], row["canonical_POS"], row["canonical_REF"], row["canonical_ALT"])
                     for row in alleles if row["normalization_status"] == "normalized" and row["original_alt"] in ("C", "A")}
    assert len(deletion_keys) == 1
    assert {row["exclusion_reason"] for row in alleles if row["normalization_status"] == "excluded"} == {
        "symbolic_alt", "breakend_alt", "spanning_deletion", "unspecified_alt", "ambiguous_sequence"}
    version = subprocess.run(["bgzip", "--version"], capture_output=True, text=True, check=True)
    version_path = selected_dir / "bgzip.version.txt"
    version_path.write_text(version.stdout + version.stderr, encoding="utf-8")
    vcf_path = selected_dir / "promising.vcf.gz"
    with vcf_path.open("wb") as handle:
        subprocess.run(["bgzip", "-c", str(selected_dir / "promising.raw.vcf")], stdout=handle, check=True)
    index_path = Path(str(vcf_path) + ".csi")
    subprocess.run(["bcftools", "index", "--csi", "--force", "-o", str(index_path), str(vcf_path)], check=True)
    exports.capture_export(config, out, selected_id, version_path)
    empty_version = external_dir / "bgzip.version.txt"
    empty_version.write_text(version.stdout + version.stderr, encoding="utf-8")
    empty_vcf = external_dir / "promising.vcf.gz"
    with empty_vcf.open("wb") as handle:
        subprocess.run(["bgzip", "-c", str(external_dir / "promising.raw.vcf")], stdout=handle, check=True)
    empty_index = Path(str(empty_vcf) + ".csi")
    subprocess.run(["bcftools", "index", "--csi", "--force", "-o", str(empty_index), str(empty_vcf)], check=True)
    exports.capture_export(config, out, external_id, empty_version)
    normalized_external = external_dir / "normalized_alleles.vcf"
    subprocess.run(["bcftools", "norm", "--no-version", "-f", str(reference), "-c", "e", "-Ov",
                    "-o", str(normalized_external), str(external_dir / "source_alleles.vcf")], check=True)
    exports.capture_alleles(config, out, external_id, normalized_external)
    for group_id in ("selected", "external"):
        group_dir = out / "vcf_outputs" / "groups" / group_id
        details = exports.write_consensus(config, out, group_id)
        assert details["member_count"] == 1 and details["required_support"] == 1
        consensus_raw = group_dir / "consensus.raw.vcf"
        consensus_vcf = group_dir / "consensus.vcf.gz"
        with consensus_vcf.open("wb") as handle:
            subprocess.run(["bgzip", "-c", str(consensus_raw)], stdout=handle, check=True)
        consensus_index = Path(str(consensus_vcf) + ".csi")
        subprocess.run(["bcftools", "index", "--csi", "--force", "-o", str(consensus_index),
                        str(consensus_vcf)], check=True)
        exports.capture_group(config, out, group_id, consensus_raw)
    complete = exports.complete_outputs(config, out)
    assert json.loads(complete.read_text(encoding="utf-8"))["group_count"] == 2
    selected_consensus = out / "vcf_outputs/groups/selected/consensus.vcf.gz"
    consensus_rows = subprocess.run(["bcftools", "view", "-H", str(selected_consensus)],
                                    capture_output=True, text=True, check=True).stdout.splitlines()
    assert len(consensus_rows) == 2
    deletion = next(line.split("\t") for line in consensus_rows if line.split("\t")[3:5] == ["CA", "C"])
    assert deletion[1] == "1" and deletion[9] == "./." and "PROMISING_SUPPORT=1" in deletion[7]
    selected_support = list(prepare.read_tsv(out / "vcf_outputs/groups/selected/support.tsv.gz"))
    deletion_support = next(row for row in selected_support if row["REF"] == "CA" and row["ALT"] == "C")
    assert deletion_support["support_count"] == "1" and deletion_support["genotype_disagreement"] == "true"
    evidence = list(prepare.read_tsv(out / "vcf_outputs/groups/selected/member_evidence.tsv.gz"))
    deletion_evidence = next(row for row in evidence if row["REF"] == "CA" and row["ALT"] == "C")
    assert deletion_evidence["presence_state"] == "present" and deletion_evidence["presence_conflict"] == "false"
    assert deletion_evidence["dosage_disagreement"] == "true" and deletion_evidence["projected_GT"] == "./."
    occurrences = json.loads(deletion_evidence["source_occurrences"])
    assert [row["target_presence"] for row in occurrences] == ["present", "present"]
    mixed_evidence = next(row for row in evidence if row["POS"] == "4" and row["ALT"] == "G")
    mixed_source = json.loads(mixed_evidence["source_occurrences"])[0]
    assert mixed_source["original_format"] == "GT:DP" and mixed_source["original_sample"] == "1/2:7"
    assert [(row["alt_index"], row["allele"]) for row in mixed_source["original_alt_alleles"]] == [
        (1, "G"), (2, "<DEL>")]
    excluded = list(prepare.read_tsv(out / "vcf_outputs/groups/selected/excluded_alleles.tsv.gz"))
    assert {row["exclusion_reason"] for row in excluded} == {
        "symbolic_alt", "breakend_alt", "spanning_deletion", "unspecified_alt", "ambiguous_sequence"}
    empty_consensus = out / "vcf_outputs/groups/external/consensus.vcf.gz"
    assert subprocess.run(["bcftools", "view", "-H", str(empty_consensus)],
                          capture_output=True, text=True, check=True).stdout == ""
    assert subprocess.run(["bcftools", "view", "-H", "-r", "chr1:1-2", str(empty_consensus)],
                          capture_output=True, text=True, check=True).stdout == ""
    assert json.loads((out / "vcf_outputs/outputs.complete.json").read_text())["dataset_count"] == 2
    query = subprocess.run(["bcftools", "view", "-H", "-r", "chr1:1-2", str(vcf_path)],
                           capture_output=True, text=True, check=True)
    assert len(query.stdout.splitlines()) == 2
    assert subprocess.run(["bcftools", "view", "-H", str(empty_vcf)],
                          capture_output=True, text=True, check=True).stdout == ""
    bad = copy.deepcopy(config)
    bad["vcf_outputs"]["groups"]["external"]["samples"] = ["unknown"]
    fails(lambda: exports.manifest(bad, out), "Unknown other_illumina members")
    datasets_path = out / "vcf_outputs/datasets.json"
    frozen_before = prepare.sha256(out / "library/freeze.json")
    datasets_before = datasets_path.read_bytes()
    annotation_hashes = [prepare.sha256(annotate.annotation_paths(out, role)[0] / "genotypes.tsv.gz")
                         for role in ("illumina", "other_illumina")]
    changed_selection = copy.deepcopy(config)
    changed_selection["selection"] = {"corresponding_dataset_support_fraction": 0.2}
    _datasets, same_groups = exports.manifest(changed_selection, out)
    assert same_groups[0]["required_support"] == 1
    assert datasets_path.read_bytes() == datasets_before
    assert prepare.sha256(out / "library/freeze.json") == frozen_before
    assert annotation_hashes == [prepare.sha256(annotate.annotation_paths(out, role)[0] / "genotypes.tsv.gz")
                                for role in ("illumina", "other_illumina")]
    no_groups = copy.deepcopy(config)
    no_groups["vcf_outputs"]["groups"] = {}
    _datasets, groups = exports.manifest(no_groups, out)
    assert groups == [] and datasets_path.read_bytes() == datasets_before
    assert json.loads(exports.complete_outputs(no_groups, out).read_text())["group_count"] == 0
    check_six_member_consensus(destination)
    print("PASS: raw VCF projection, annotation identity, indel normalization, exclusions, BGZF/CSI and empty export")



def check_six_member_consensus(destination):
    import annotate
    import prepare
    import vcf_outputs as exports

    out = destination / "six-member-consensus"
    _frozen_fixture(out)
    reference = out / "reference.fa"
    reference.write_text(">chr1\nCAAAAG\n", encoding="utf-8")
    Path(str(reference) + ".fai").write_text("chr1\t6\t6\t6\t7\n", encoding="utf-8")
    selected = "IL-fixed"
    samples = [f"OTHER-{i}" for i in range(1, 6)] + ["OTHER,6"]
    freeze_path = out / "library/freeze.json"
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    freeze["reference_sha256"] = prepare.sha256(reference)
    freeze["selection_config"] = {"illumina_corresponding": {"PB-test": [selected]}}
    prepare.atomic_json(freeze_path, freeze)
    corresponding = out / "corresponding.vcf"
    external = out / "external.vcf"
    corresponding_rows = [
        ["chr1", "3", "selected", "A", "C", "30", "PASS", "AC=1;AN=2", "GT:DP", "0/1:10"]
    ]
    external_rows = [
        ["chr1", "3", "five", "A", "C", "31", "PASS", "AC=5;AN=12", "DP:GT",
         *[f"8:{'0/1' if i < 5 else '0/0'}" for i in range(6)]],
        ["chr1", "4", "four", "A", "G", "32", "PASS", "AC=4;AN=12", "DP:GT",
         *[f"8:{'0/1' if i < 4 else '0/0'}" for i in range(6)]],
        ["chr1", "5", "duplicate-a", "A", "T", "33", "PASS", "AC=6;AN=12", "DP:GT",
         *["8:0/1" for _ in range(6)]],
        ["chr1", "5", "duplicate-b", "A", "T", "34", "PASS", "AC=6;AN=12", "DP:GT",
         *[f"8:{'0/1' if i < 4 else ('1/1' if i == 4 else '0/0')}" for i in range(6)]],
        ["chr1", "6", "missing-and-conflict-a", "G", "A", "35", "PASS", "AC=5;AN=10", "DP:GT",
         *[f"8:{'0/1' if i < 4 or i == 5 else './.'}" for i in range(6)]],
        ["chr1", "6", "missing-and-conflict-b", "G", "A", "36", "PASS", "AC=4;AN=12", "DP:GT",
         *[f"8:{'0/1' if i < 4 else ('./.' if i == 4 else '0/0')}" for i in range(6)]],
    ]

    def write_source(path, names, rows):
        with path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=6>\n")
            handle.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">\n')
            handle.write('##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">\n')
            handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n')
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
                         "\t".join(names) + "\n")
            for row in rows:
                handle.write("\t".join(row) + "\n")

    write_source(corresponding, [selected], corresponding_rows)
    write_source(external, samples, external_rows)
    sources = {"illumina": corresponding, "other_illumina": external}
    source_rows = {"illumina": corresponding_rows, "other_illumina": external_rows}
    identities = {"illumina": "frozen_selection_manifest", "other_illumina": "header_names_only"}
    contexts = {"illumina": "selection_reuse", "other_illumina": "external_application_identity_unverified"}
    for role, names in (("illumina", [selected]), ("other_illumina", samples)):
        source = sources[role]
        source_hash = prepare.sha256(source)
        app, _ = annotate.annotation_paths(out, role)
        app.mkdir(parents=True, exist_ok=True)
        annotation = app / "genotypes.tsv.gz"
        annotation_rows = []
        datasets = []
        for name in names:
            did = prepare.dataset_id(role, name)
            datasets.append(dict(dataset_id=did, role=role, sample_name=name, source_vcf=str(source),
                                 corresponding_pacbio="PB-test" if role == "illumina" else "",
                                 identity_check_status=identities[role],
                                 selection_participant=str(role == "illumina").lower(),
                                 evidence_context=contexts[role]))
            for ordinal, fields in enumerate(source_rows[role], 1):
                keys, values = fields[8].split(":"), fields[9 + names.index(name)].split(":")
                gt = values[keys.index("GT")]
                record = {field: "" for field in annotate.FIELDS}
                record.update(source_vcf=str(source), source_record_index=str(ordinal),
                              record_id=f"{role}_{source_hash}_{ordinal}", dataset_id=did,
                              sample_name=name, role=role,
                              selection_participant=str(role == "illumina").lower(),
                              evidence_context=contexts[role], CHROM=fields[0], POS=fields[1],
                              ID=fields[2], REF=fields[3], ALT=fields[4], GT=gt,
                              annotation_status="promising", identity_check_status=identities[role],
                              extraction_status="eligible")
                annotation_rows.append(record)
        prepare.write_tsv(annotation, annotation_rows, annotate.FIELDS)
        application = dict(schema_version=1, datasets=datasets, source_vcf=str(source),
                           source_sha256=source_hash, identity_check_status=identities[role],
                           reference=str(reference), flank_bp=2)
        prepare.atomic_json(app / "manifest.json", application)
        complete = app / "annotation.complete.json"
        prepare.atomic_json(complete, {"schema_version": 1, "role": role,
                                       "freeze_sha256": prepare.sha256(freeze_path)})
        prepare.publish_provenance(complete, [source, app / "manifest.json", annotation],
                                   {"fixture": True}, [complete, annotation], scope="application")

    config = {
        "output_dir": str(out), "illumina_vcf": str(corresponding), "other_illumina": str(external),
        "illumina_corresponding": {"PB-test": [selected]},
        "vcf_outputs": {"consensus_fraction": 0.75,
                        "groups": {"six_external": {"role": "other_illumina", "samples": samples}}},
    }
    datasets, groups = exports.manifest(config, out)
    assert len(datasets) == 7 and groups[0]["n"] == 6 and groups[0]["required_support"] == 5
    version = subprocess.run(["bgzip", "--version"], capture_output=True, text=True, check=True)
    for dataset in datasets:
        did = dataset["dataset_id"]
        exports.sample_inputs(config, out, did)
        sample_dir = out / "vcf_outputs" / did
        version_path = sample_dir / "bgzip.version.txt"
        version_path.write_text(version.stdout + version.stderr, encoding="utf-8")
        vcf = sample_dir / "promising.vcf.gz"
        with vcf.open("wb") as handle:
            subprocess.run(["bgzip", "-c", str(sample_dir / "promising.raw.vcf")], stdout=handle, check=True)
        index = Path(str(vcf) + ".csi")
        subprocess.run(["bcftools", "index", "--csi", "--force", "-o", str(index), str(vcf)], check=True)
        exports.capture_export(config, out, did, version_path)
        if dataset["role"] == "other_illumina":
            normalized = sample_dir / "normalized_alleles.vcf"
            subprocess.run(["bcftools", "norm", "--no-version", "-f", str(reference), "-c", "e", "-Ov",
                            "-o", str(normalized), str(sample_dir / "source_alleles.vcf")], check=True)
            exports.capture_alleles(config, out, did, normalized)
    details = exports.write_consensus(config, out, "six_external")
    group_dir = out / "vcf_outputs/groups/six_external"
    consensus_raw = group_dir / "consensus.raw.vcf"
    consensus_vcf = group_dir / "consensus.vcf.gz"
    with consensus_vcf.open("wb") as handle:
        subprocess.run(["bgzip", "-c", str(consensus_raw)], stdout=handle, check=True)
    consensus_index = Path(str(consensus_vcf) + ".csi")
    subprocess.run(["bcftools", "index", "--csi", "--force", "-o", str(consensus_index),
                    str(consensus_vcf)], check=True)
    group_complete = exports.capture_group(config, out, "six_external", consensus_raw)
    summary = {row["POS"]: row for row in prepare.read_tsv(group_dir / "support.tsv.gz")}
    assert summary["3"]["support_count"] == "5" and summary["3"]["retained"] == "true"
    assert summary["3"]["support_percent"] == "83.333333"
    assert json.loads(summary["3"]["supporting_samples"]) == samples[:5]
    assert summary["4"]["support_count"] == "4" and summary["4"]["retained"] == "false"
    assert summary["4"]["support_percent"] == "66.666667"
    assert summary["5"]["support_count"] == "5" and summary["5"]["retained"] == "true"
    assert summary["5"]["genotype_disagreement"] == "true"
    assert summary["6"]["group_n"] == "6" and summary["6"]["required_support"] == "5"
    assert json.loads(summary["5"]["conflicting_samples"]) == ["OTHER,6"]
    assert json.loads(summary["6"]["conflicting_samples"]) == ["OTHER,6"]
    assert summary["6"]["support_count"] == "4" and summary["6"]["retained"] == "false"
    consensus_rows = subprocess.run(["bcftools", "view", "-H", str(consensus_vcf)],
                                    capture_output=True, text=True, check=True).stdout.splitlines()
    assert [line.split("\t")[1] for line in consensus_rows] == ["3", "5"]
    assert "PROMISING_SUPPORT=5" in consensus_rows[0] and "GT_DISAGREEMENT" in consensus_rows[0]
    assert "CONFLICTING_SAMPLES=OTHER%2C6" in consensus_rows[1]
    assert consensus_rows[1].split("\t")[14] == "./."
    assert len(consensus_rows[0].split("\t")) == 15
    evidence = list(prepare.read_tsv(group_dir / "member_evidence.tsv.gz"))
    missing = next(row for row in evidence if row["POS"] == "6" and row["sample_name"] == "OTHER-5")
    conflict = next(row for row in evidence if row["POS"] == "6" and row["sample_name"] == "OTHER,6")
    assert missing["presence_state"] == "unknown" and missing["support"] == "false"
    assert conflict["presence_state"] == "conflict" and conflict["presence_conflict"] == "true"
    conflict_five = next(row for row in evidence if row["POS"] == "5" and row["sample_name"] == "OTHER,6")
    assert conflict_five["presence_state"] == "conflict" and conflict_five["support"] == "false"
    dosage = next(row for row in evidence if row["POS"] == "5" and row["sample_name"] == "OTHER-5")
    assert dosage["projected_GT"] == "./."
    assert dosage["support"] == "true" and dosage["dosage_disagreement"] == "true"
    completed_hash = prepare.sha256(group_complete)
    member_marker = out / "vcf_outputs" / prepare.dataset_id("other_illumina", samples[-1]) / "alleles.complete.json"
    saved_marker = member_marker.read_bytes()
    member_marker.unlink()
    fails(lambda: exports.write_consensus(config, out, "six_external"), "Incomplete or stale artifact")
    fails(lambda: exports.complete_outputs(config, out), "Incomplete or stale artifact")
    assert prepare.sha256(group_complete) == completed_hash
    member_marker.write_bytes(saved_marker)
    global_marker = exports.complete_outputs(config, out)
    assert json.loads(global_marker.read_text(encoding="utf-8"))["group_count"] == 1
    print("PASS: six-member ceil threshold, fixed denominator, exact supporters and incomplete-member guard")
def check_vcf_logic():
    import vcf_outputs as exports

    def call(gt, ordinal=1, promised=True, targets=(1,), alt_count=2):
        return dict(source_record_index=str(ordinal), GT=gt, alt_count=alt_count,
                    target_indexes=set(targets), annotation_status="promising" if promised else "unclassified",
                    allele_identities={2: ("canonical", "chr1", "4", "A", "G")})

    for original, indexes, count, expected in (
        ("0/2", {1}, 2, "0/."),
        ("1/2", {1}, 2, "1/."),
        ("2/2", {1}, 2, "./."),
        ("1/.", {1}, 2, "1/."),
        ("0/.", {1}, 2, "0/."),
        ("0|2/1", {1}, 2, "0|./1"),
        ("1/2", {1, 2}, 2, "1/1"),
        ("0/1/2", {1}, 2, "0/1/."),
    ):
        assert exports.project_gt(original, indexes, count) == expected
    fails(lambda: exports.project_gt("0/3", {1}, 2), "Out-of-range GT")

    both_present = exports.resolve_member([call("0/1", 1), call("1/1", 2)])
    assert both_present["support"] and not both_present["presence_conflict"]
    assert both_present["dosage_disagreement"] and both_present["duplicate_genotype_disagreement"]
    assert both_present["projected_GT"] == "./."

    contradiction = exports.resolve_member([call("0/1", 1), call("0/0", 2)])
    assert not contradiction["support"] and contradiction["presence_conflict"]
    assert contradiction["presence_state"] == "conflict" and contradiction["projected_GT"] == "./."

    known_missing = exports.resolve_member([call("0/1", 1), call("./.", 2)])
    assert known_missing["support"] and not known_missing["presence_conflict"]
    assert known_missing["missingness"] and not known_missing["dosage_disagreement"]
    assert known_missing["projected_GT"] == "./."

    identical = exports.resolve_member([call("0/1", 1), call("0/1", 2)])
    assert identical["support"] and identical["projected_GT"] == "0/1"
    assert not identical["duplicate_genotype_disagreement"]

    unknown = exports.resolve_member([call("0/.", 1), call("./.", 2)])
    assert not unknown["support"] and unknown["presence_state"] == "unknown"
    assert unknown["missingness"] and not unknown["presence_conflict"]
    assert unknown["projected_GT"] == "./."

    identical_partial = exports.resolve_member([call("0/.", 1), call("0/.", 2)])
    assert identical_partial["projected_GT"] == "0/." and not identical_partial["support"]

    ploidy = exports.resolve_member([call("0/1", 1), call("1", 2)])
    assert ploidy["support"] and ploidy["ploidy_disagreement"]
    assert ploidy["duplicate_genotype_disagreement"] and ploidy["projected_GT"] == "."

    partial_target = exports.resolve_member([call("1/.", 1), call("./.", 2)])
    assert partial_target["support"] and partial_target["projected_GT"] == "./."
    assert not partial_target["presence_conflict"]

    non_target = exports.resolve_member([call("2/2", 1)])
    assert non_target["presence_state"] == "absent" and non_target["call_state"] == "other_alt_only"
    assert non_target["projected_GT"] == "./." and not non_target["support"]

    unpromising = exports.resolve_member([call("0/1", 1, promised=False)])
    assert not unpromising["support"] and unpromising["presence_state"] == "present"
    no_occurrence = exports.resolve_member([])
    assert no_occurrence["call_state"] == "absent" and no_occurrence["presence_state"] == "unknown"
    mixed_labels = exports.resolve_member([call("0/1", 1, promised=True), call("0/1", 2, promised=False)])
    assert mixed_labels["support"] and mixed_labels["annotation_disagreement"]
    unknown_promising = exports.resolve_member([call("0/.", 1, promised=True)])
    assert not unknown_promising["support"] and unknown_promising["presence_state"] == "unknown"
    collapsed = exports.resolve_member([call("1/2", 1, targets=(1, 2))])
    assert collapsed["support"] and collapsed["projected_GT"] == "1/1" and collapsed["duplicate_count"] == 1
    print("PASS: target projection, duplicate presence policy, unresolved GTs and independent support votes")

def check_science():
    import evidence
    import annotate
    import math
    hit, distant = {"start": 48, "stop": 59}, {"start": 4, "stop": 15}
    positive = evidence.metric([hit, distant], [], 4, 50, 10)
    assert positive["real_coverage"] == 0.5 and math.isinf(positive["coverage_ratio"])
    assert evidence.support(positive, 1) and not evidence.support(positive, 1, central=True)
    high = {"start": 4, "stop": 15, "strand": "-", "score": "20", "p-value": "0.00001"}
    low = {"start": 48, "stop": 59, "strand": "+", "score": "19", "p-value": "0.000001"}
    assert annotate.winner_key(low, 1) < annotate.winner_key(high, 2)
    assert evidence.complete_link_families(["C", "A", "B"], {frozenset(("A", "B")), frozenset(("B", "C"))})
    print("PASS: core motif evidence and deterministic annotation winner")


def _load_config(path):
    import yaml
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def _run_snakemake(config_path, targets, workflow_scope=None, config_overrides=()):
    config = _load_config(config_path)
    command = [
        "snakemake", *targets, "--snakefile", str(ROOT / "Snakefile"),
        "--directory", str(ROOT), "--configfile", str(config_path),
        "--profile", "none", "--executor", "local", "--cores", "2",
    ]
    overrides = ([f"workflow_scope={workflow_scope}"] if workflow_scope else []) + list(config_overrides)
    if overrides:
        command.extend(["--config", *overrides])
    if config["execution"]["use_containers"]:
        command.extend([
            "--software-deployment-method", "apptainer",
            "--apptainer-prefix", config["execution"]["image_cache"],
        ])
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    if result.stdout:
        sys.stdout.write(result.stdout)
    if result.stderr:
        sys.stderr.write(result.stderr)
    return result


def _selection_snapshot(config):
    import prepare
    out = Path(config["output_dir"])
    paths = [out / "library/freeze.json", out / "manifests/selection_fingerprints.json"]
    for sample in config["pacbio_samples"]:
        rotation = out / "evidence" / prepare.dataset_id("pacbio", sample) / "r01"
        paths.extend([rotation / "candidates.json", rotation / "transfer.complete.json"])
    snapshot = {}
    for path in paths:
        assert path.is_file(), f"Missing frozen-selection output: {path}"
        snapshot[path] = (path.stat().st_mtime_ns, prepare.sha256(path))
    return snapshot


def _check_workflow_outputs(config):
    import annotate
    import prepare
    import vcf_outputs as exports
    out = Path(config["output_dir"])
    freeze = out / "library/freeze.json"
    assert freeze.is_file()
    for role, directory in (
        ("other_illumina", out / "application"),
        ("illumina", out / "application/corresponding"),
    ):
        complete = directory / "annotation.complete.json"
        annotations = directory / "genotypes.tsv.gz"
        exports.verify_artifact(complete, [complete, annotations])
        rows = iter(prepare.read_tsv(annotations))
        first = next(rows)
        assert first["role"] == role
        if role == "illumina":
            assert first["selection_participant"] == "true"
            assert first["evidence_context"] == "selection_reuse"
        else:
            assert first["selection_participant"] == "false"

    datasets = json.loads((out / "vcf_outputs/datasets.json").read_text(encoding="utf-8"))["datasets"]
    groups = json.loads((out / "vcf_outputs/groups.json").read_text(encoding="utf-8"))["groups"]
    assert {row["role"] for row in datasets} == {"illumina", "other_illumina"}
    assert [(group["group_id"], group["n"], group["required_support"]) for group in groups] == [
        ("corresponding_a", 6, 5), ("corresponding_b", 5, 4), ("external", 6, 5)]
    for row in datasets:
        annotation_dir, _ = annotate.annotation_paths(out, row["role"])
        assert Path(row["annotation_path"]) == annotation_dir / "genotypes.tsv.gz"
        assert Path(row["annotation_complete"]).is_file()
        dataset_id = row["dataset_id"]
        sample = out / "vcf_outputs" / dataset_id
        marker = sample / "export.complete.json"
        vcf = sample / "promising.vcf.gz"
        index = Path(str(vcf) + ".csi")
        exports.verify_artifact(marker, [marker, vcf, index])
        subprocess.run(
            ["bcftools", "view", "-H", "-r", "segment_A:1-40000", str(vcf)],
            check=True, capture_output=True, text=True,
        )
    for group in groups:
        directory = out / "vcf_outputs/groups" / group["group_id"]
        marker = directory / "complete.json"
        vcf = directory / "consensus.vcf.gz"
        index = Path(str(vcf) + ".csi")
        artifacts = [marker, vcf, index, *(directory / name for name in (
            "support.tsv.gz", "member_evidence.tsv.gz", "excluded_alleles.tsv.gz"))]
        exports.verify_artifact(marker, artifacts)
        subprocess.run(
            ["bcftools", "view", "-H", "-r", "segment_A:1-40000", str(vcf)],
            check=True, capture_output=True, text=True,
        )
    complete = exports.complete_outputs(config, out)
    summary = json.loads(complete.read_text(encoding="utf-8"))
    assert summary["dataset_count"] == len(datasets) and summary["group_count"] == len(groups)


def _run_application_only(config_path):
    import annotate
    config = _load_config(config_path)
    out = Path(config["output_dir"])
    annotate.verify_freeze(out)
    before = _selection_snapshot(config)
    missing_selection = str(out / "unavailable-selection.vcf.gz")
    missing_external = str(out / "unavailable-external.vcf.gz")
    for target, overrides in (
        ("annotate_other", [f"pacbio_vcf={missing_selection}", f"illumina_vcf={missing_selection}"]),
        ("annotate_corresponding", [f"other_illumina={missing_external}"]),
    ):
        result = _run_snakemake(
            config_path, [target], workflow_scope="application", config_overrides=overrides)
        if result.returncode:
            raise RuntimeError(f"Application-only target failed without unrelated source: {target}")
        assert _selection_snapshot(config) == before, "Application-only target changed frozen-selection outputs"
    result = _run_snakemake(
        config_path,
        ["annotate_corresponding", "export_promising", "consensus_groups", "vcf_deliverables"],
        workflow_scope="application",
    )
    if result.returncode:
        raise RuntimeError("Frozen-library application-only run failed")
    assert _selection_snapshot(config) == before, "Application-only run changed frozen-selection outputs"
    _check_workflow_outputs(config)


def _check_full_workflow_guards(config_path):
    import prepare
    import vcf_outputs as exports
    config = _load_config(config_path)
    out = Path(config["output_dir"])
    groups = json.loads((out / "vcf_outputs/groups.json").read_text(encoding="utf-8"))["groups"]
    external_group = next(group for group in groups if group["role"] == "other_illumina")
    directory = out / "vcf_outputs/groups" / external_group["group_id"]
    support = directory / "support.tsv.gz"
    saved_support = support.read_bytes()
    support.write_bytes(saved_support + b"x")
    fails(lambda: exports.complete_outputs(config, out))
    support.write_bytes(saved_support)
    exports.complete_outputs(config, out)

    member_id = external_group["member_dataset_ids"][-1]
    sample = out / "vcf_outputs" / member_id
    sites = sample / "source_alleles.vcf"
    lines = sites.read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(lines):
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        fields[3] = next(base for base in "ACGT" if base != fields[3])
        lines[index] = "\t".join(fields)
        break
    else:
        raise AssertionError("No concrete allele available for the invalid-REF guard")
    sites.write_text("\n".join(lines) + "\n", encoding="utf-8")
    for name in (
        "normalized_alleles.vcf", "alleles.tsv.gz", "alleles.tsv.gz.provenance.json",
        "alleles.complete.json", "alleles.complete.json.provenance.json",
    ):
        (sample / name).unlink(missing_ok=True)
    previous_group_marker = prepare.sha256(directory / "complete.json")
    result = _run_snakemake(config_path, ["consensus_groups"], workflow_scope="application")
    assert result.returncode != 0, "Consensus accepted a member with invalid-reference normalization"
    assert prepare.sha256(directory / "complete.json") == previous_group_marker
    fails(lambda: exports.complete_outputs(config, out))


def check_full_workflow(destination, mode, image_cache):
    config, config_path = fixture(destination, mode, image_cache)
    result = _run_snakemake(config_path, ["all"])
    if result.returncode:
        raise RuntimeError("Full rule-all workflow failed")
    _check_workflow_outputs(config)
    _run_application_only(config_path)
    _check_full_workflow_guards(config_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["containers", "modules"], default="modules")
    parser.add_argument("--image-cache", type=Path)
    parser.add_argument("--keep-fixtures", type=Path)
    parser.add_argument("--logic-only", action="store_true")
    parser.add_argument("--vcf-outputs-only", action="store_true")
    parser.add_argument("--application-only", type=Path)
    args = parser.parse_args()
    if args.application_only:
        if args.keep_fixtures:
            raise ValueError("--application-only and --keep-fixtures cannot be combined")
        config_path = args.application_only / "config.yaml" if args.application_only.is_dir() else args.application_only
        if not config_path.is_file():
            raise FileNotFoundError(config_path)
        _run_application_only(config_path)
        return

    def run_checks(destination):
        if args.vcf_outputs_only:
            check_vcf_logic()
            check_vcf_outputs(destination)
            return
        check_configuration()
        check_application_binds()
        check_preparation(destination)
        check_science()
        check_vcf_logic()
        check_annotation_roles(destination)
        if args.logic_only:
            return
        check_full_workflow(destination, args.mode, args.image_cache)

    if args.keep_fixtures:
        args.keep_fixtures.mkdir(parents=True, exist_ok=True)
        if not args.logic_only and not args.vcf_outputs_only and any(args.keep_fixtures.iterdir()):
            raise ValueError("--keep-fixtures must be an empty directory for a full workflow run")
        run_checks(args.keep_fixtures)
    else:
        with tempfile.TemporaryDirectory(prefix="motif-consumer-") as temp:
            run_checks(Path(temp))


if __name__ == "__main__":
    main()
