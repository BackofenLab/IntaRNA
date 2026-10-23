# Optional standalone IntaRNApvalue executable

From the repository root, in a Python 3.10+ environment:

```bash
python3 -m pip install ./python pyinstaller
python3 python/bin/build.py
python/bin/build/IntaRNApvalue --help
```

The builder uses the invoking interpreter, resolves source paths independently
of the working directory, and creates `python/bin/build/IntaRNApvalue` (with an
`.exe` suffix on Windows). Build metadata is removed automatically. An existing
executable with that name is replaced by PyInstaller. Build on each target OS;
no precompiled binaries are shipped. The resulting executable still needs the
native IntaRNA program beside it or on PATH.
