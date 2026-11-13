# targetlist-generator

- Package in `src/targetlist_generator/`; `pipeline.make_targetlist` is the single entry for `cli.py` and `gui.py` (PySide6, optional extra).
- Test: `pytest -m "not network"` (GUI tests need `QT_QPA_PLATFORM=offscreen` headless); lint: `ruff check src tests`.
- Output file formats (IRTF starcat, ING `.cat`) are consumed by telescope software: keep columns and separators unchanged.
- Fail loudly on bad input; no silent fallbacks.
- `CITATION.cff` is tracked by Zenodo: bump its version and date only in the commit preparing a full release.
- British spelling in docs and messages.
