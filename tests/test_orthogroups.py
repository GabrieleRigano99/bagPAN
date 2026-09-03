from pathlib import Path

import pytest

from bagpan.orthogroups import OrthogroupsError, classify_orthogroups, parse_orthogroups_txt

FIXTURES = Path(__file__).parent / "fixtures"

LOCUS_PREFIX_TO_SPECIES = {"SPA": "species_a", "SPB": "species_b", "SPC": "species_c"}


def test_parse_orthogroups_txt():
    clusters = parse_orthogroups_txt(FIXTURES / "Orthogroups.txt", LOCUS_PREFIX_TO_SPECIES)
    assert set(clusters) == {"OG0000001", "OG0000002", "OG0000003", "OG0000004"}
    assert clusters["OG0000001"] == ["SPA_000001-T1", "SPB_000001-T1", "SPC_000001-T1"]
    assert clusters["OG0000003"] == ["SPA_000003-T1"]


def test_parse_orthogroups_txt_rejects_unknown_prefix():
    with pytest.raises(OrthogroupsError):
        parse_orthogroups_txt(FIXTURES / "Orthogroups.txt", {"SPA": "species_a"})


def test_classify_orthogroups_core_accessory_singleton():
    clusters = parse_orthogroups_txt(FIXTURES / "Orthogroups.txt", LOCUS_PREFIX_TO_SPECIES)
    result = classify_orthogroups(
        clusters, LOCUS_PREFIX_TO_SPECIES, ["species_a", "species_b", "species_c"]
    )
    assert result.category["OG0000001"] == "Core"
    assert result.category["OG0000002"] == "Accessory"
    assert result.category["OG0000003"] == "Singleton"
    assert result.category["OG0000004"] == "Singleton"
    assert result.species_present["OG0000001"] == {"species_a", "species_b", "species_c"}
    assert result.species_present["OG0000002"] == {"species_a", "species_b"}