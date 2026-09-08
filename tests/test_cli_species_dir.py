from pathlib import Path

from bagpan.cli import _discover_species_dir


def test_discover_species_dir_lists_subdirectories_sorted(tmp_path):
    (tmp_path / "species_b").mkdir()
    (tmp_path / "species_a").mkdir()
    (tmp_path / "species_c").mkdir()
    (tmp_path / "not_a_species.txt").write_text("x")  # files are skipped
    (tmp_path / ".hidden_dir").mkdir()  # hidden dirs are skipped

    pairs = _discover_species_dir(tmp_path)

    assert pairs == [
        ("species_a", str(tmp_path / "species_a")),
        ("species_b", str(tmp_path / "species_b")),
        ("species_c", str(tmp_path / "species_c")),
    ]


def test_discover_species_dir_empty(tmp_path):
    assert _discover_species_dir(tmp_path) == []


def test_discover_species_dir_returns_string_paths(tmp_path):
    (tmp_path / "sx").mkdir()
    name, path = _discover_species_dir(tmp_path)[0]
    assert name == "sx"
    assert isinstance(path, str)
    assert Path(path) == tmp_path / "sx"
