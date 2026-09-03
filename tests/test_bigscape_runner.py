import sqlite3
from pathlib import Path

from bagpan.bigscape_runner import _extract_gcf_mapping, stage_bgc_regions
from bagpan.discover import SpeciesFiles


def _make_db(tmp_path: Path, table: str, id_col: str, family_col: str, rows) -> Path:
    db_path = tmp_path / "results.db"
    conn = sqlite3.connect(str(db_path))
    conn.execute(f"CREATE TABLE {table} ({id_col} TEXT, {family_col} TEXT)")
    conn.executemany(f"INSERT INTO {table} VALUES (?, ?)", rows)
    # an unrelated table that shouldn't confuse the introspection
    conn.execute("CREATE TABLE run_parameters (key TEXT, value TEXT)")
    conn.commit()
    conn.close()
    return db_path


def test_extract_gcf_mapping_finds_recognizable_schema(tmp_path):
    db_path = _make_db(
        tmp_path, "gcf_membership", "bgc_id", "family_id",
        [("region001", "GCF_1"), ("region002", "GCF_1"), ("region003", "GCF_2")],
    )
    mapping = _extract_gcf_mapping(db_path)
    assert mapping == {"region001": "GCF_1", "region002": "GCF_1", "region003": "GCF_2"}


def test_extract_gcf_mapping_unrecognized_schema_returns_empty(tmp_path, capsys):
    db_path = tmp_path / "results.db"
    conn = sqlite3.connect(str(db_path))
    conn.execute("CREATE TABLE something_else (foo TEXT, bar TEXT)")
    conn.commit()
    conn.close()

    mapping = _extract_gcf_mapping(db_path)
    assert mapping == {}
    assert "WARNING" in capsys.readouterr().out


def test_extract_gcf_mapping_handles_corrupt_db(tmp_path, capsys):
    db_path = tmp_path / "not_a_db.db"
    db_path.write_text("this is not a sqlite database")
    mapping = _extract_gcf_mapping(db_path)
    assert mapping == {}
    assert "WARNING" in capsys.readouterr().out


def test_stage_bgc_regions_prefixes_by_species_and_skips_missing(tmp_path):
    species_a_regions = tmp_path / "species_a_regions"
    species_a_regions.mkdir()
    (species_a_regions / "ctg1.region001.gbk").write_text("LOCUS fake\n")
    (species_a_regions / "ctg1.region002.gbk").write_text("LOCUS fake\n")

    sf_a = SpeciesFiles(
        name="species_a", root=tmp_path, proteome=tmp_path / "a.fa", locus_prefix="A",
        antismash_regions_dir=species_a_regions,
    )
    sf_b = SpeciesFiles(
        name="species_b", root=tmp_path, proteome=tmp_path / "b.fa", locus_prefix="B",
        antismash_regions_dir=None,
    )

    input_dir = tmp_path / "input"
    usable = stage_bgc_regions([sf_a, sf_b], input_dir)

    assert [sf.name for sf in usable] == ["species_a"]
    staged = sorted(p.name for p in input_dir.glob("*.gbk"))
    assert staged == ["species_a__ctg1.region001.gbk", "species_a__ctg1.region002.gbk"]