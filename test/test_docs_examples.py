"""The documentation's examples are executable, and the pages really include them.

Prose cannot be trusted to stay true on its own. Every code block on the site is a real
module under `docs/examples`, run by this test and asserting its own claims, so a renamed
argument or a changed answer fails here rather than misleading a reader.

Two failure modes look different and are covered separately:

* an example stops working -- caught by running it;
* a page stops *including* an example -- `pymdownx.snippets` fails open on an unknown
  region, rendering an empty code block in a build that still succeeds, so the references
  are checked against the files they name.
"""
import importlib.util
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / 'docs'
EXAMPLES = DOCS / 'examples'

# Mirrors pymdownx.snippets.base_path in zensical.toml
SNIPPET_BASE_PATHS = (EXAMPLES, ROOT)

# --8<-- "file.py"  or  --8<-- "file.py:region"
SNIPPET = re.compile(r'^\s*--8<--\s+"(?P<path>[^":]+)(?::(?P<region>[\w-]+))?"', re.M)

# `# --8<-- [start:name]` / `[end:name]`
REGION = re.compile(r'--8<--\s+\[(?P<edge>start|end):(?P<name>[\w-]+)\]')

docs_present = pytest.mark.skipif(
    not EXAMPLES.is_dir(), reason='the documentation sources are not in this checkout')


def example_paths():
    return sorted(EXAMPLES.glob('*.py')) if EXAMPLES.is_dir() else []


def page_paths():
    return sorted(DOCS.rglob('*.md')) if DOCS.is_dir() else []


def resolve(path):
    for base in SNIPPET_BASE_PATHS:
        candidate = base / path
        if candidate.is_file():
            return candidate
    return None


@docs_present
@pytest.mark.parametrize('path', example_paths(), ids=lambda p: p.stem)
def test_example_runs(path):
    """Every example executes and its own assertions hold.

    An example resting on numpy skips where numpy is not installed -- the wheel test job
    installs only pytest -- and only for that reason: the skip is decided by which module
    the ImportError names, so an example broken any other way still fails here.
    """
    spec = importlib.util.spec_from_file_location(f'docs_example_{path.stem}', path)
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except ImportError as error:
        missing = (error.name or '').split('.')[0]
        if missing in {'numpy', 'scipp'}:
            pytest.skip(f'{path.name} needs {missing}, which is not installed')
        raise


@docs_present
@pytest.mark.parametrize('page', page_paths(), ids=lambda p: str(p.relative_to(DOCS)))
def test_every_snippet_reference_resolves(page):
    """A page including a file or a region that does not exist renders an empty block.

    `pymdownx.snippets` does not fail the build for it, so nothing else would notice.
    """
    for match in SNIPPET.finditer(page.read_text()):
        target = resolve(match.group('path'))
        assert target is not None, (
            f'{page.relative_to(DOCS)} includes {match.group("path")!r}, which is not '
            f'under any snippet base path')
        region = match.group('region')
        if region is None:
            continue
        text = target.read_text()
        for edge in ('start', 'end'):
            assert f'--8<-- [{edge}:{region}]' in text, (
                f'{page.relative_to(DOCS)} includes region {region!r} of {target.name}, '
                f'which has no [{edge}:{region}] marker')


@docs_present
@pytest.mark.parametrize('path', example_paths(), ids=lambda p: p.stem)
def test_every_region_is_used_somewhere(path):
    """An example region nothing includes is either dead or a page forgot it."""
    pages = ' '.join(page.read_text() for page in page_paths())
    defined = {m.group('name') for m in REGION.finditer(path.read_text())}
    for name in sorted(defined):
        assert f'{path.name}:{name}"' in pages, (
            f'{path.name} defines region {name!r} and no page includes it')


@docs_present
def test_every_example_is_referenced():
    """A whole example nothing includes would go stale unnoticed."""
    pages = ' '.join(page.read_text() for page in page_paths())
    for path in example_paths():
        assert f'"{path.name}:' in pages or f'"{path.name}"' in pages, (
            f'{path.name} is not included by any page')


@docs_present
def test_regions_are_balanced():
    """Every start has an end, so an include cannot run to the bottom of the file."""
    for path in example_paths():
        edges = [(m.group('edge'), m.group('name')) for m in REGION.finditer(path.read_text())]
        starts = [name for edge, name in edges if edge == 'start']
        ends = [name for edge, name in edges if edge == 'end']
        assert sorted(starts) == sorted(ends), f'{path.name} has unbalanced regions'
        assert len(starts) == len(set(starts)), f'{path.name} repeats a region name'

# --- the pages only use Markdown the configuration actually enables -------------------
#
# A feature whose extension is not loaded does not fail the build: `:material-wrench:`
# renders as the literal text `:material-wrench:`, an admonition renders as a paragraph
# beginning `!!!`, and the site publishes looking wrong. That is how the icons on the home
# page reached production. These check the sources against zensical.toml instead of the
# built HTML, so they need no build.

FEATURES = (
    (re.compile(r':(?:material|fontawesome|octicons|simple)-[a-z0-9-]+:'),
     'pymdownx.emoji', 'icon shortcodes'),
    (re.compile(r'^\s*!!!\s', re.M), 'admonition', 'admonitions'),
    (re.compile(r'^\s*\?\?\?\s', re.M), 'pymdownx.details', 'collapsible admonitions'),
    (re.compile(r'^\s*--8<--\s', re.M), 'pymdownx.snippets', 'snippet includes'),
    (re.compile(r'^\s*\|.+\|\s*$', re.M), 'tables', 'tables'),
    (re.compile(r'<[a-z]+[^>]*\bmarkdown\b'), 'md_in_html', 'Markdown inside HTML'),
    (re.compile(r'^\s*=== "', re.M), 'pymdownx.tabbed', 'content tabs'),
    (re.compile(r'^\s*```mermaid', re.M), 'pymdownx.superfences', 'mermaid diagrams'),
)

ICON = re.compile(r':((?:material|fontawesome|octicons|simple)-[a-z0-9-]+):')


def enabled_extensions():
    """Every extension named in zensical.toml, as dotted names.

    Dotted TOML keys nest, so `pymdownx.emoji.emoji_index` arrives as
    {'pymdownx': {'emoji': {'emoji_index': ...}}}; flatten it back.
    """
    import tomllib
    with open(ROOT / 'zensical.toml', 'rb') as file:
        config = tomllib.load(file)
    names = set()

    def walk(mapping, prefix=()):
        for key, value in mapping.items():
            path = (*prefix, key)
            names.add('.'.join(path))
            if isinstance(value, dict):
                walk(value, path)

    walk(config['project']['markdown_extensions'])
    return names


@docs_present
@pytest.mark.parametrize('page', page_paths(), ids=lambda p: str(p.relative_to(DOCS)))
def test_page_only_uses_enabled_markdown(page):
    enabled = enabled_extensions()
    text = page.read_text()
    for pattern, extension, description in FEATURES:
        match = pattern.search(text)
        if match is None:
            continue
        assert any(name.startswith(extension) for name in enabled), (
            f'{page.relative_to(DOCS)} uses {description} ({match.group(0).strip()!r}) but '
            f'{extension} is not enabled in zensical.toml, so it will render as literal '
            f'text in a build that still succeeds')


@docs_present
def test_every_icon_shortcode_names_a_real_icon():
    """A misspelled icon renders as its own shortcode, and the build says nothing."""
    zensical = pytest.importorskip('zensical', reason='needs the documentation toolchain')
    icons = Path(zensical.__file__).parent / 'templates' / '.icons'
    if not icons.is_dir():
        pytest.skip('this zensical does not ship an icon set where expected')
    for page in page_paths():
        for match in ICON.finditer(page.read_text()):
            family, _, name = match.group(1).partition('-')
            assert (icons / family / f'{name}.svg').is_file(), (
                f'{page.relative_to(DOCS)} uses :{match.group(1)}:, which is not in '
                f'{icons}/{family}')
