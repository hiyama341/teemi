# Releasing

This repository uses `versioneer` and Git tags for versioning.

Important details for this repo:

- tags are plain version numbers like `1.0.2`
- do not use a `v` prefix unless the versioning setup is changed
- PyPI publishing is triggered by pushing a tag

## Pre-release checklist

1. Make sure your working tree is clean enough for release.
2. Run the test suite.
3. Build the package locally.
4. Confirm the computed version looks correct.
5. Create and push a Git tag for the release.
6. Verify the GitHub Actions publish job succeeds.

## Check the current computed version

Use:

```bash
python setup.py version
```

Because this repo uses `versioneer`, the version is derived from Git tags.

If `HEAD` is already on the release tag, the version should resolve to that exact tag.
If `HEAD` is ahead of a tag, the version will include local version metadata.

## Recommended release flow

Run tests:

```bash
pytest
```

Build the distributions:

```bash
python setup.py sdist bdist_wheel
```

Create the release tag:

```bash
git tag 1.0.3
```

Push the branch and the tag:

```bash
git push origin main
git push origin 1.0.3
```

Once the tag is pushed, GitHub Actions should publish the package to PyPI using the configured `PYPI_TOKEN_TEEMI` secret.

## Verify the workflow configuration

Publishing is controlled by:

- [`.github/workflows/main.yml`](./.github/workflows/main.yml)

The publish step runs only when the pushed ref is a tag.

## UV-based developer workflow

The package can still be released with the existing build commands, but `uv` can be used for local environment management.

Example:

```bash
uv venv
source .venv/bin/activate
uv pip install -e .[dev]
uv run pytest
```

If `uv` becomes the standard developer workflow later, the release build steps can also be updated.

## Notes

- Keep release tags consistent with the historical format already used in the repo.
- Avoid tagging from a dirty tree.
- If packaging metadata changes, test both install and build before tagging.
