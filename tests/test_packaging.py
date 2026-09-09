"""Build, validate, and install release artifacts without modifying the checkout."""
from pathlib import Path
import os
import re
import shutil
import subprocess
import sys
import tarfile
import venv
import zipfile

from sectalix import __version__

ROOT = Path(__file__).resolve().parents[1]
DOCS = ("THEORY_AND_CONVENTIONS.md", "CLI_REFERENCE.md", "VERIFICATION_BENCHMARKS.md")


def checked(args: list[str], cwd: Path, env: dict[str, str]) -> str:
    proc = subprocess.run(args, cwd=cwd, env=env, text=True, capture_output=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return proc.stdout


def test_version_and_entrypoint() -> None:
    metadata = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert f'version = "{__version__}"' in metadata
    assert __version__ == "1.0.1"
    assert 'sectalix = "sectalix.cli:main"' in metadata


def test_no_unresolved_source_markers() -> None:
    for source in (ROOT / "src/sectalix").glob("*.py"):
        assert not re.search(r"\b(TODO|FIXME|NotImplementedError)\b", source.read_text(encoding="utf-8")), source
    # Returning NotImplemented from binary operator dispatch is required Python behavior.


def test_documentation_links() -> None:
    for name in DOCS:
        path = ROOT / "docs" / name
        text = path.read_text(encoding="utf-8")
        assert len(text) > 2500
        for target in re.findall(r"\]\(([^)]+)\)", text):
            if "://" not in target and not target.startswith("#"):
                assert (path.parent / target.split("#")[0]).exists(), target


def test_build_install_wheel_and_sdist(tmp_path: Path) -> None:
    project = tmp_path / "project"
    project.mkdir()
    for name in ("pyproject.toml", "README.md", "MANIFEST.in", "LICENSE"):
        if (ROOT / name).exists():
            shutil.copy2(ROOT / name, project / name)
    for name in ("src", "docs", "schemas", "tests"):
        shutil.copytree(ROOT / name, project / name, ignore=shutil.ignore_patterns("__pycache__", "*.egg-info"))
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["PYTHONIOENCODING"] = "utf-8"
    # Only the known CPython 3.10/3.11 Windows distutils import warning is exempt.
    # Other messages, categories, modules and Python versions remain errors.
    build_code = "import sys, warnings\n"
    build_code += (
        "if sys.platform == 'win32' and (3, 10) <= sys.version_info[:2] < (3, 12):\n"
        "    warnings.filterwarnings('ignore', "
        "message=r'^The distutils package is deprecated and slated for removal in Python 3\\.12\\. ', "
        "category=DeprecationWarning, module=r'^distutils$')\n"
        "from setuptools.build_meta import build_wheel, build_sdist\n"
        "build_wheel('dist'); build_sdist('dist')\n"
    )
    checked([sys.executable, "-W", "error", "-c", build_code], project, env)
    wheel = next((project / "dist").glob("*.whl"))
    archive = next((project / "dist").glob("*.tar.gz"))
    with zipfile.ZipFile(wheel) as z:
        assert "sectalix/cli.py" in z.namelist()
        assert not any(name.startswith("thinwallx/") for name in z.namelist())
        assert any(n.endswith("THEORY_AND_CONVENTIONS.md") for n in z.namelist())
    with tarfile.open(archive) as t:
        assert any(n.endswith("tests/test_cli.py") for n in t.getnames())
        assert any(n.endswith("schemas/sectalix-0.8.schema.json") for n in t.getnames())
        assert not any("/src/thinwallx/" in name for name in t.getnames())
    target = tmp_path / "installed"
    venv.EnvBuilder(with_pip=True, symlinks=(os.name != "nt"), system_site_packages=True).create(target)
    binary = target / ("Scripts" if os.name == "nt" else "bin")
    python = binary / ("python.exe" if os.name == "nt" else "python")
    checked([str(python), "-m", "pip", "install", "--no-deps", "--no-index", str(wheel)], tmp_path, env)
    console = binary / ("sectalix.exe" if os.name == "nt" else "sectalix")
    assert checked([str(console), "--version"], tmp_path, env).strip() == "Sectalix 1.0.1"
    assert "inspect" in checked([str(python), "-m", "sectalix", "--help"], tmp_path, env)
    imported = checked([str(python), "-c", "import sectalix; print(sectalix.__file__)"], tmp_path, env)
    assert target.resolve().as_posix().lower() in Path(imported.strip()).resolve().as_posix().lower()
