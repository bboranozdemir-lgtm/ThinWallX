"""T31–T36: build and install real artifacts, without modifying the checkout."""
from pathlib import Path
import os
import re
import shutil
import subprocess
import sys
import tarfile
import venv
import zipfile

from thinwallx import __version__

ROOT = Path(__file__).resolve().parents[1]
DOCS = ("THEORY_AND_CONVENTIONS.md", "CLI_REFERENCE.md", "VERIFICATION_BENCHMARKS.md")


def checked(args: list[str], cwd: Path, env: dict[str, str]) -> str:
    proc = subprocess.run(args, cwd=cwd, env=env, text=True, capture_output=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return proc.stdout


def test_t31_t32_version_entrypoint() -> None:
    metadata = (ROOT / "pyproject.toml").read_text()
    assert f'version = "{__version__}"' in metadata
    assert __version__ == "1.0.0"
    assert 'thinwallx = "thinwallx.cli:main"' in metadata


def test_t34_no_unresolved_markers() -> None:
    for source in (ROOT / "src/thinwallx").glob("*.py"):
        assert not re.search(r"\b(TODO|FIXME|NotImplementedError)\b", source.read_text(encoding="utf-8")), source
    # Returning NotImplemented from binary operator dispatch is required Python behavior.


def test_t35_documentation_links() -> None:
    for name in DOCS:
        path = ROOT / "docs" / name
        text = path.read_text(encoding="utf-8")
        assert len(text) > 2500
        for target in re.findall(r"\]\(([^)]+)\)", text):
            if "://" not in target and not target.startswith("#"):
                assert (path.parent / target.split("#")[0]).exists(), target


def test_t33_t36_build_install_wheel_and_sdist(tmp_path: Path) -> None:
    project = tmp_path / "project"
    project.mkdir()
    for name in ("pyproject.toml", "README.md", "MANIFEST.in"):
        shutil.copy2(ROOT / name, project / name)
    for name in ("src", "docs", "schemas", "tests"):
        shutil.copytree(ROOT / name, project / name, ignore=shutil.ignore_patterns("__pycache__", "*.egg-info"))
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["PYTHONIOENCODING"] = "utf-8"
    checked([sys.executable, "-W", "error", "-c",
             "from setuptools.build_meta import build_wheel, build_sdist; build_wheel('dist'); build_sdist('dist')"], project, env)
    wheel = next((project / "dist").glob("*.whl"))
    archive = next((project / "dist").glob("*.tar.gz"))
    with zipfile.ZipFile(wheel) as z:
        assert "thinwallx/cli.py" in z.namelist()
        assert any(n.endswith("THEORY_AND_CONVENTIONS.md") for n in z.namelist())
    with tarfile.open(archive) as t:
        assert any(n.endswith("tests/test_cli.py") for n in t.getnames())
        assert any(n.endswith("schemas/thinwallx-0.8.schema.json") for n in t.getnames())
    target = tmp_path / "installed"
    venv.EnvBuilder(with_pip=True, system_site_packages=True).create(target)
    binary = target / ("Scripts" if os.name == "nt" else "bin")
    python = binary / ("python.exe" if os.name == "nt" else "python")
    checked([str(python), "-m", "pip", "install", "--no-deps", "--no-index", str(wheel)], tmp_path, env)
    console = binary / ("thinwallx.exe" if os.name == "nt" else "thinwallx")
    assert checked([str(console), "--version"], tmp_path, env).strip() == "ThinWallX 1.0.0"
    assert "inspect" in checked([str(python), "-m", "thinwallx", "--help"], tmp_path, env)
    imported = checked([str(python), "-c", "import thinwallx; print(thinwallx.__file__)"], tmp_path, env)
    assert str(target).lower() in imported.lower()
