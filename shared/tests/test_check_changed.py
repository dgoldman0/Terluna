"""Selection must preserve affected consumers without replaying unrelated work.

Small temporary repositories exercise Git changes and dependency chains. No
scientific model, external data or complete repository suite runs in these tests.
"""
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from shared import check_changed as checks


@pytest.fixture
def repo(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    put(tmp_path, ".gitignore", ".pytest_cache/\n")
    put(tmp_path, "pyproject.toml", '''[tool.pytest.ini_options]
testpaths = ["shared", "atmosphere", "biosphere", "climate", "geography", "resources", "research/studies"]
norecursedirs = ["tools", "archive", ".*"]
''')
    commit(tmp_path)
    return tmp_path


def put(root, path, content):
    target = root / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(content)


def commit(root):
    subprocess.run(["git", "add", "."], cwd=root, check=True)
    subprocess.run(["git", "-c", "user.name=Check fixture", "-c", "user.email=check@example.invalid",
                    "commit", "-qm", "Fixture state"], cwd=root, check=True)


def select(root, *changed):
    files, _ = checks.inventory(root, "HEAD")
    return checks.plan(root, files | set(changed), set(changed))


def test_prose_links_do_not_select_the_linked_models(repo):
    put(repo, "README.md", "See atmosphere/model.py for the calculation.\n")
    put(repo, "atmosphere/model.py", "VALUE = 1\n")
    put(repo, "atmosphere/tests/test_model.py", "from atmosphere import model\n")
    result = select(repo, "README.md")
    assert not result.python and not result.targets and not result.full


def test_document_read_by_one_test_selects_only_that_test(repo):
    put(repo, "research/decisions.md", "A population range.\n")
    put(repo, "resources/tests/test_population.py", '''from pathlib import Path
def test_population():
    assert 'range' in Path('research/decisions.md').read_text()
def test_expensive_simulation():
    raise AssertionError('must not be selected for a prose edit')
''')
    result = select(repo, "research/decisions.md")
    assert result.python == ["resources/tests/test_population.py::test_population"]
    assert result.targets == []


def test_module_level_document_input_selects_the_whole_file(repo):
    put(repo, "research/decisions.md", "A population range.\n")
    put(repo, "resources/tests/test_population.py", '''from pathlib import Path
CONTENT = Path('research/decisions.md').read_text()
def test_population():
    assert CONTENT
''')
    assert select(repo, "research/decisions.md").python == ["resources/tests/test_population.py"]


def test_planning_json_avoids_owner_fallback_but_still_selects_actual_readers(repo):
    path = "habitation/crater_gardens/buildout.json"
    put(repo, path, "{}\n")
    put(repo, "resources/tests/test_unrelated.py", "def test_other(): pass\n")
    assert not select(repo, path).python
    put(repo, "resources/tests/test_buildout.py", f"INPUT = {path!r}\n")
    assert select(repo, path).python == ["resources/tests/test_buildout.py"]


def test_routing_tables_do_not_make_selector_a_scientific_consumer(repo):
    put(repo, "research/status.json", "{}\n")
    put(repo, "shared/check_changed.py", "RECORDS = {'research/status.json'}\n")
    put(repo, "shared/tests/test_check_changed.py", "from shared import check_changed\n")
    assert not select(repo, "research/status.json").python


def test_transitive_imports_reach_consumers_without_sibling_tests(repo):
    put(repo, "atmosphere/thermo.py", "VALUE = 1\n")
    put(repo, "biosphere/growth.py", "from atmosphere import thermo\n")
    put(repo, "biosphere/tests/test_growth.py", "from biosphere import growth\n")
    put(repo, "biosphere/tests/test_unrelated.py", "def test_other(): pass\n")
    result = select(repo, "atmosphere/thermo.py")
    assert result.python == ["biosphere/tests/test_growth.py"]
    assert result.reasons[result.python[0]] == [
        "atmosphere/thermo.py", "biosphere/growth.py", "biosphere/tests/test_growth.py"]


def test_relative_imports_and_package_initializers(repo):
    put(repo, "atmosphere/air/__init__.py", "VALUE = 1\n")
    put(repo, "atmosphere/air/thermo.py", "VALUE = 1\n")
    put(repo, "atmosphere/air/model.py", "from . import thermo\n")
    put(repo, "atmosphere/tests/test_air.py", "import atmosphere.air.model\n")
    for changed in ("atmosphere/air/thermo.py", "atmosphere/air/__init__.py"):
        assert select(repo, changed).python == ["atmosphere/tests/test_air.py"]


def test_product_metadata_keys_connect_producers_to_downstream_tests(repo):
    put(repo, "atmosphere/thermo.py", "VALUE = 1\n")
    put(repo, "atmosphere/results/air.json", json.dumps({
        "producer": {"files": {"atmosphere/thermo.py": "digest"}}}))
    put(repo, "research/studies/habitat/results/habitat.json", json.dumps({
        "producer": {"inputs": {"atmosphere/results/air.json": "digest"}}}))
    put(repo, "research/studies/habitat/run.py", '''from pathlib import Path
HERE = Path(__file__).resolve().parent
OUT = HERE / 'results' / 'habitat.json'
''')
    put(repo, "research/studies/habitat/test_run.py", "from research.studies.habitat import run\n")
    put(repo, "research/studies/unrelated/test_run.py", "def test_other(): pass\n")
    result = select(repo, "atmosphere/thermo.py")
    assert result.python == ["research/studies/habitat/test_run.py"]
    assert "atmosphere/results/air.json" in result.reasons[result.python[0]]
    assert "research/studies/habitat/results/habitat.json" in result.reasons[result.python[0]]


def test_json_citation_is_not_a_dependency_but_a_pinned_document_is(repo):
    put(repo, "research/decisions.md", "A population range.\n")
    put(repo, "atmosphere/results/air.json", json.dumps({
        "source": "research/decisions.md", "description": "See research/decisions.md for the decision.",
        "producer": {"cited_files": ["research/decisions.md"]}}))
    put(repo, "atmosphere/tests/test_air.py", "INPUT = 'atmosphere/results/air.json'\n")
    assert not select(repo, "research/decisions.md").python
    put(repo, "atmosphere/results/air.json", json.dumps({
        "producer": {"inputs": {"research/decisions.md": "digest"}}}))
    assert select(repo, "research/decisions.md").python == ["atmosphere/tests/test_air.py"]


def test_historical_verification_hashes_are_not_live_product_inputs(repo):
    put(repo, "research/status.json", "{}\n")
    put(repo, "research/studies/habitat/checks.json", json.dumps({
        "verification_sources": {"research/status.json": "a past digest"}}))
    put(repo, "research/studies/habitat/test_checks.py", "INPUT = 'research/studies/habitat/checks.json'\n")
    assert not select(repo, "research/status.json").python


def test_source_citations_are_not_reads(repo):
    put(repo, "research/decisions.md", "A population range.\n")
    put(repo, "biosphere/taxa.py", "CITED = 'research/decisions.md'\nSOURCE = 'See research/decisions.md for the names.'\n")
    put(repo, "biosphere/tests/test_taxa.py", "from biosphere import taxa\n")
    assert not select(repo, "research/decisions.md").python


def test_named_document_path_read_is_a_dependency(repo):
    put(repo, "research/decisions.md", "A population range.\n")
    put(repo, "biosphere/taxa.py", '''from pathlib import Path
SOURCE = Path('research/decisions.md')
def read_names():
    return SOURCE.read_text()
''')
    put(repo, "biosphere/tests/test_taxa.py", "from biosphere import taxa\n")
    assert select(repo, "research/decisions.md").python == ["biosphere/tests/test_taxa.py"]


def test_bare_star_is_not_an_input_glob(repo):
    put(repo, "atmosphere/input.json", "{}\n")
    put(repo, "atmosphere/tests/test_names.py", "NAMES = {'DATA', '*'}\n")
    files, _ = checks.inventory(repo, "HEAD")
    reverse = checks.graph(repo, files, checks.test_files(repo, files))
    assert "atmosphere/tests/test_names.py" not in reverse.get("atmosphere/input.json", set())


@pytest.mark.parametrize("pattern", ["*.json", "*/*.json", "*/snapshot.json"])
def test_runtime_glob_does_not_link_unrelated_products(repo, pattern):
    put(repo, "atmosphere/results/snapshot.json", "{}\n")
    put(repo, "atmosphere/tests/test_glob.py", f'''def test_values():
    return external_runtime_directory.glob({pattern!r})
''')
    files, _ = checks.inventory(repo, "HEAD")
    reverse = checks.graph(repo, files, checks.test_files(repo, files))
    assert not reverse.get("atmosphere/results/snapshot.json")


def test_glob_with_known_directory_selects_its_readers(repo):
    put(repo, "atmosphere/results/case1.json", "{}\n")
    put(repo, "atmosphere/tests/test_glob.py", '''from pathlib import Path
RESULTS = Path(__file__).resolve().parents[1] / 'results'
def test_values():
    return list(RESULTS.glob('case*.json'))
''')
    assert select(repo, "atmosphere/results/case1.json").python == ["atmosphere/tests/test_glob.py::test_values"]


def test_external_makefile_identity_does_not_bind_repository_makefile(repo):
    put(repo, "Makefile", "check:\n\ttrue\n")
    put(repo, "climate/results/build.json", json.dumps({"files": {"Makefile": "external digest"}}))
    put(repo, "climate/build.py", "PATCH = ('Makefile', 'external model')\n")
    files, _ = checks.inventory(repo, "HEAD")
    reverse = checks.graph(repo, files, checks.test_files(repo, files))
    assert not reverse.get("Makefile")


def test_joined_paths_resolve_from_repository_root(repo):
    put(repo, "geography/results/water.json", "{}\n")
    put(repo, "biosphere/lake.py", '''from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / 'geography' / 'results' / 'water.json'
''')
    put(repo, "biosphere/tests/test_lake.py", "from biosphere import lake\n")
    assert select(repo, "geography/results/water.json").python == ["biosphere/tests/test_lake.py"]


def test_ambiguous_basenames_select_all_possible_consumers(repo):
    put(repo, "geography/results/water.json", "{}\n")
    put(repo, "climate/results/water.json", "{}\n")
    put(repo, "biosphere/tests/test_lake.py", "INPUT = 'water.json'\n")
    for path in ("geography/results/water.json", "climate/results/water.json"):
        assert select(repo, path).python == ["biosphere/tests/test_lake.py"]


def test_source_and_its_test_changed_together_do_not_expand_to_owner(repo):
    put(repo, "atmosphere/thermo.py", "VALUE = 1\n")
    put(repo, "atmosphere/tests/test_thermo.py", "from atmosphere import thermo\n")
    put(repo, "atmosphere/tests/test_unrelated.py", "def test_other(): pass\n")
    result = select(repo, "atmosphere/thermo.py", "atmosphere/tests/test_thermo.py")
    assert result.python == ["atmosphere/tests/test_thermo.py"]


def test_unresolved_code_falls_back_to_owner_suite(repo):
    put(repo, "climate/new_solver.f90", "! a compiled dependency\n")
    put(repo, "climate/tests/test_build.py", "def test_build(): pass\n")
    put(repo, "biosphere/tests/test_growth.py", "def test_growth(): pass\n")
    result = select(repo, "climate/new_solver.f90")
    assert result.python == ["climate/tests/test_build.py"]
    assert "owner fallback: climate" in result.reasons[result.python[0]]


@pytest.mark.parametrize("path", ["shared/constants.json", "shared/provenance.py", "pyproject.toml",
                                  "research/requirements.txt", "unknown-global.cfg"])
def test_shared_and_unknown_global_changes_require_full_suite(repo, path):
    assert select(repo, path).full


def test_conftest_affects_only_descendant_tests(repo):
    put(repo, "atmosphere/rain/conftest.py", "VALUE = 1\n")
    put(repo, "atmosphere/rain/test_rain.py", "def test_rain(): pass\n")
    put(repo, "atmosphere/tests/test_other.py", "def test_other(): pass\n")
    assert select(repo, "atmosphere/rain/conftest.py").python == ["atmosphere/rain/test_rain.py"]


def test_javascript_viewer_does_not_bake_unrelated_immersion(repo):
    put(repo, "visualization/atmospheric-columns/viewer.mjs", "export const value = 1;\n")
    result = select(repo, "visualization/atmospheric-columns/viewer.mjs")
    assert result.targets == ["check-layers", "test-column-viewer"]
    assert not result.full


def test_javascript_imports_follow_a_python_produced_input(repo):
    put(repo, "geography/terrain.py", "VALUE = 1\n")
    put(repo, "geography/results/terrain.json", json.dumps({"files": {"geography/terrain.py": "digest"}}))
    put(repo, "immersion/bake/terrain.mjs", "const input = '../../geography/results/terrain.json';\n")
    result = select(repo, "geography/terrain.py")
    assert "test-immersion" in result.targets


def test_fixture_data_and_heavy_excluded_harnesses_are_not_collected(repo):
    put(repo, "atmosphere/tools/test_heavy.py", "def test_heavy(): pass\n")
    put(repo, "atmosphere/tests/test_small.py", "def test_small(): pass\n")
    files, _ = checks.inventory(repo, "HEAD")
    assert checks.test_files(repo, files) == {"atmosphere/tests/test_small.py"}


def test_git_inventory_includes_staged_unstaged_untracked_and_deleted(repo):
    for name in ("staged.md", "unstaged.md", "deleted.md"):
        put(repo, name, "before\n")
    commit(repo)
    put(repo, "staged.md", "after\n")
    subprocess.run(["git", "add", "staged.md"], cwd=repo, check=True)
    put(repo, "unstaged.md", "after\n")
    put(repo, "new file.md", "new\n")
    (repo / "deleted.md").unlink()
    files, changed = checks.inventory(repo, "HEAD")
    assert changed == {"staged.md", "unstaged.md", "deleted.md", "new file.md"}
    assert "deleted.md" in files


def test_deleted_module_still_selects_surviving_importers(repo):
    put(repo, "atmosphere/thermo.py", "VALUE = 1\n")
    put(repo, "atmosphere/tests/test_thermo.py", "from atmosphere import thermo\n")
    commit(repo)
    (repo / "atmosphere/thermo.py").unlink()
    files, changed = checks.inventory(repo, "HEAD")
    assert checks.plan(repo, files, changed).python == ["atmosphere/tests/test_thermo.py"]


def test_base_revision_checks_committed_work_and_both_sides_of_rename(repo):
    put(repo, "before.md", "text\n")
    commit(repo)
    base = checks.git(repo, "rev-parse", "HEAD").strip()
    subprocess.run(["git", "mv", "before.md", "after.md"], cwd=repo, check=True)
    commit(repo)
    assert checks.inventory(repo, "HEAD")[1] == set()
    assert checks.inventory(repo, base)[1] == {"before.md", "after.md"}


def test_explicit_ignored_input_reaches_its_consumer(repo, monkeypatch, capsys):
    put(repo, ".gitignore", "research/runs/\n")
    put(repo, "research/runs/forcing.json", "{}\n")
    put(repo, "biosphere/tests/test_forcing.py", "INPUT = 'research/runs/forcing.json'\n")
    commit(repo)
    monkeypatch.setattr(checks, "ROOT", repo)
    assert checks.main(["--json", "--include", "research/runs/forcing.json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["python"] == ["biosphere/tests/test_forcing.py"]


def test_plan_does_not_execute_checks(repo, monkeypatch):
    monkeypatch.setattr(checks, "ROOT", repo)
    marker = repo / "executed"
    monkeypatch.setattr(checks, "commands", lambda *_: [[sys.executable, "-c",
                       "from pathlib import Path; Path('executed').touch()"]])
    assert checks.main(["--plan"]) == 0
    assert not marker.exists()


def test_failed_check_stops_and_returns_failure(repo, monkeypatch):
    monkeypatch.setattr(checks, "ROOT", repo)
    monkeypatch.setattr(checks, "commands", lambda *_: [
        [sys.executable, "-c", "raise SystemExit(7)"],
        [sys.executable, "-c", "from pathlib import Path; Path('incorrectly-continued').touch()"],
    ])
    assert checks.main([]) == 1
    assert not (repo / "incorrectly-continued").exists()


def test_malformed_json_fails_before_running_selected_checks(repo, monkeypatch):
    put(repo, "research/status.json", "{ invalid json\n")
    monkeypatch.setattr(checks, "ROOT", repo)
    assert checks.main([]) == 1


def test_bad_python_cannot_silently_drop_dependencies(repo, monkeypatch):
    put(repo, "atmosphere/broken.py", "def broken(\n")
    monkeypatch.setattr(checks, "ROOT", repo)
    assert checks.main(["--plan"]) == 1


def test_dependency_cache_invalidates_on_content_and_new_files(repo):
    put(repo, "atmosphere/a.py", "VALUE = 1\n")
    put(repo, "atmosphere/tests/test_a.py", "from atmosphere import a\n")
    assert select(repo, "atmosphere/a.py").python == ["atmosphere/tests/test_a.py"]
    put(repo, "atmosphere/tests/test_b.py", "from atmosphere import a\n")
    assert select(repo, "atmosphere/a.py").python == ["atmosphere/tests/test_a.py", "atmosphere/tests/test_b.py"]
    stamp = (repo / "atmosphere/tests/test_b.py").stat()
    put(repo, "atmosphere/tests/test_b.py", "def test_unrelated(): pass\n")
    os.utime(repo / "atmosphere/tests/test_b.py", ns=(stamp.st_atime_ns, stamp.st_mtime_ns))
    assert select(repo, "atmosphere/a.py").python == ["atmosphere/tests/test_a.py"]


def test_bad_dependency_cache_is_rebuilt(repo):
    put(repo, "atmosphere/a.py", "VALUE = 1\n")
    put(repo, "atmosphere/tests/test_a.py", "from atmosphere import a\n")
    put(repo, ".pytest_cache/check_changed/dependencies.json", "{broken\n")
    assert select(repo, "atmosphere/a.py").python == ["atmosphere/tests/test_a.py"]


def test_full_target_retains_all_original_checks_and_local_target_uses_selector():
    root = Path(checks.__file__).resolve().parents[1]
    full = subprocess.check_output(["make", "-n", "check-full"], cwd=root, text=True)
    for command in ("shared/check_layers.py", "-m pytest", "npm test", "atmosphere/column/tests/",
                    "illumination/references/tests/", "visualization/atmospheric-columns/tests/",
                    "visualization/light-calendar/build.py", "illumination/calendar/tests/",
                    "research/check.py", "ensemble/tools/validate_setup.py"):
        assert command in full
    assert "check_changed.py" not in full
    local = subprocess.check_output(["make", "-n", "check"], cwd=root, text=True)
    assert "shared/check_changed.py" in local
    assert "-m pytest" not in local and "npm test" not in local
