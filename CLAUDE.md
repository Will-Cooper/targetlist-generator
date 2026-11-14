# targetlist-generator

- Package in `src/targetlist_generator/`; `pipeline.make_targetlist` is the single entry for `cli.py` and `gui.py` (PySide6, optional extra).
- Test: `pytest -m "not network"` (GUI tests need `QT_QPA_PLATFORM=offscreen` headless); lint: `ruff check src tests`.
- Output formats (`formats.py`: IRTF, ING `.cat`, Keck, Lick, ESO p2, Subaru, MMT, Magellan, ...) are read by telescope software: each writer follows its observatory's documented spec (URL in its docstring); keep columns and separators unchanged, and do not add a format whose spec cannot be confirmed.
- The planner never drops a visible target (extras are listed at their best time); standards prefer short slews (`--std-max-sep`) over perfect airmass matches.
- Fail loudly on bad input; no silent fallbacks.
- `CITATION.cff` is tracked by Zenodo: bump its version and date only in the commit preparing a full release.
- British spelling in docs and messages.
