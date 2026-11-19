# Releasing

This repository uses `versioneer` and Git tags for versioning.

Important details for this repo:

- tags are plain version numbers like `1.0.2`
- do not use a `v` prefix unless the versioning setup is changed
- PyPI publishing is triggered by pushing a tag

## Pre-release checklist

1. Cut the `Unreleased` section of `HISTORY.rst` to the version being released.
2. Make sure your working tree is clean enough for release.
3. Run the test suite.
4. Build the package locally.
5. Confirm the computed version looks correct.
6. Land the work on `main` through a pull request.
7. Create and push a Git tag for the release, from `main`.
8. Verify the GitHub Actions publish job succeeds.

`HISTORY.rst` is rendered as the second half of the PyPI long description and a
published version cannot be re-uploaded, so step 1 has to be right before the
tag exists. Rename `Unreleased` to the version and its date, and open a fresh
`Unreleased` heading above it.

## Check the current computed version

Use:

```bash
uv run python -c 'import teemi; print(teemi.__version__)'
```

Because this repo uses `versioneer`, the version is derived from Git tags.

If `HEAD` is already on the release tag, the version should resolve to that exact tag.
If `HEAD` is ahead of a tag, the version will include local version metadata.

## Recommended release flow

Run tests:

```bash
uv run --locked --extra dev pytest tests
```

Build the distributions:

```bash
uv build
```

Land the work on `main`. Branch protection on `main` requires signed commits
and one approving review, and forbids force pushes, so `git push origin main` is
rejected for commits made without a signing key. Merge the pull request with the
GitHub merge button instead, which signs the merge commit with GitHub's own key:

```bash
git push origin your-branch
```

Then tag the merge commit from an up-to-date `main`:

```bash
git switch main
git pull
git tag 1.2.0
git push origin 1.2.0
```

Once the tag is pushed, GitHub Actions should publish the package to PyPI using the configured `PYPI_TOKEN_TEEMI` secret.

## Verify the workflow configuration

Publishing is controlled by:

- [`.github/workflows/main.yml`](./.github/workflows/main.yml)

The publish step runs only when the pushed ref is a tag, and only from the
Python 3.11 matrix entry. Note that nothing requires the tagged commit to be on
`main`: the workflow triggers on any push, so a tag pushed from a feature branch
publishes that branch to PyPI. Always tag from `main`.

## UV-based developer workflow

Use the committed lockfile for local development and CI. uv manages `.venv` and installs Teemi in editable mode.

Example:

```bash
uv sync --locked --extra dev
uv run --locked --extra dev pytest tests
uv build
```

Dependencies live in `pyproject.toml`. After changing them, run `uv lock` and commit both `pyproject.toml` and `uv.lock`.

## Notes

- Keep release tags consistent with the historical format already used in the repo.
- Avoid tagging from a dirty tree.
- If packaging metadata changes, test both install and build before tagging.
