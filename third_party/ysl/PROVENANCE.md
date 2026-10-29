# Vendored upstream package: `ysl`

- Upstream: `/Users/colmclaude/claude_experiments/YEAST_SLs/src/ysl` (repository YEAST_SLs, commit `9b5dce913342e8da7f763bed96d568a017f19712`; upstream working tree clean at copy time)
- Copied: 2026-10-29 with `rsync -a --exclude __pycache__` (read-only on the source). Licence: same as upstream (no licence file upstream; internal research use).
- Upstream SHA-256 below = `git show 9b5dce9:src/ysl/<file>` (blob at the pinned commit); vendored SHA-256 = file in this directory.

## Patches (only path roots; no model, label, split or evaluation logic changed)

Why: upstream derives every root from the package location (`Path(__file__).resolve().parents[N]` = the YEAST_SLs
checkout). In a vendored copy that would point at humanSLmapping, and the runner would write results/ledger/predictions
inside YEAST_SLs. The patches make data/config reads go to the read-only upstream root (from `config/paths.yaml`) and
any runner output go to humanSLmapping.

1. **New file `_paths.py`** (not upstream): resolves `UPSTREAM_ROOT` (env `YSL_UPSTREAM_ROOT`, else `hsl.paths` →
   `config/paths.yaml: yeast_sls_root`) and `OUTPUT_ROOT` (env `YSL_OUTPUT_ROOT`, else `hsl.paths`); importing `hsl` also
   installs the YEAST_SLs write guard. Refuses an OUTPUT_ROOT inside UPSTREAM_ROOT.
2. `data/costanzo.py`, `data/features.py`, `data/sgd.py`, `data/expression_meta.py`, `eval/splits.py`, `features/pairs.py`,
   `features/cofit.py`, `models/tier2.py`, `predict.py`: the single line `ROOT = Path(__file__).resolve().parents[N]` →
   `from ysl._paths import UPSTREAM_ROOT as ROOT` (one-line patch each, marked `# humanSLmapping PATCH`).
3. `runner.py`: `ROOT = UPSTREAM_ROOT` (data guard + recipes, read-only); `RESULTS = OUTPUT_ROOT/results/ysl_runner`;
   `PRED_DIR = OUTPUT_ROOT/data/ysl_predictions`; new `UPSTREAM_PRED_DIR` kept in the method data guard's protected set
   (`_protected`), so methods still cannot open upstream prediction files. Marked `# humanSLmapping PATCH`.

Diff against upstream: `diff -r <(git -C ../YEAST_SLs show 9b5dce9:src/ysl/...) third_party/ysl/...` — every changed line carries `humanSLmapping PATCH`.

## Files

| file | upstream SHA-256 (9b5dce9) | vendored SHA-256 | status |
|---|---|---|---|
| `__init__.py` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | unmodified |
| `_paths.py` | `—` | `1ebf24a818e5279ca89f523b410b1150fd11423f041abfd35c6e966144115734` | **new** (humanSLmapping) |
| `data/__init__.py` | `8ae03e258ddc933f275fed159a0cca996dfcbeedc27c5e47d88cf56427c92589` | `8ae03e258ddc933f275fed159a0cca996dfcbeedc27c5e47d88cf56427c92589` | unmodified |
| `data/costanzo.py` | `ff9a6c6811858cf53abc1e533383ca473ed3621deb0d3a574f7db161a155360f` | `b497f25c90376a3ea1c3b7ef9dff476122f3ab0573ab5a8dd9e882abf14c3dd4` | **patched** (path roots) |
| `data/dense.py` | `9dc58d90094f2cf3ce4cf932c92ddc45bfe77650a7b759fced10c10a8050ae0a` | `9dc58d90094f2cf3ce4cf932c92ddc45bfe77650a7b759fced10c10a8050ae0a` | unmodified |
| `data/expression_meta.py` | `97566351ea481d0d63cba433f937d0416025ac774bc7da6fadbb42033fa8c8f6` | `a60e0c799c60b84a2fc1edf19d94ca1cc115beae336450a7553cfb2f7b866dcd` | **patched** (path roots) |
| `data/features.py` | `8d904d93747be2105c2d91df6a8328c4096621d02c1027aa84f305053999976a` | `f392f3b33b992f42e9bfbdef6a8261e2bde8c6053f88b5871afcf37730b98bc9` | **patched** (path roots) |
| `data/sgd.py` | `4bcaa84da6d678758e71aef0be4417f574c6b00e7badd8c6627615289cf4518b` | `8229827142bca9eac999f4546d7d1cdc554a0717b81d96d7e900b33634b8c2e7` | **patched** (path roots) |
| `eval/__init__.py` | `03c7f2b8ac1f5abf2c3e6331caa7a2b0fe0ca71129990fbd9eaef01c1dc32f85` | `03c7f2b8ac1f5abf2c3e6331caa7a2b0fe0ca71129990fbd9eaef01c1dc32f85` | unmodified |
| `eval/evaluator.py` | `2d61915bb13194be46d1f33559d82c90d437070a825e4a1af81a6b71a8e9b273` | `2d61915bb13194be46d1f33559d82c90d437070a825e4a1af81a6b71a8e9b273` | unmodified |
| `eval/labels.py` | `b4ca5e5a86630c860b4e770b8a5402c67294b811e32678dcafb5f6583f842eab` | `b4ca5e5a86630c860b4e770b8a5402c67294b811e32678dcafb5f6583f842eab` | unmodified |
| `eval/metrics.py` | `0fee3c310009136c54b1b8752ee996c6cc7861390bd0515e9f29d989148a31aa` | `0fee3c310009136c54b1b8752ee996c6cc7861390bd0515e9f29d989148a31aa` | unmodified |
| `eval/splits.py` | `0011d94a748e454cf072dd8916f67b67711885a90899315b5e7ea4d49c48b336` | `3aea6ce53622177219df77b52aa43c6b69ec6013bc7b801b32a39bdee3997209` | **patched** (path roots) |
| `features/__init__.py` | `a81ff8bf5f32c06f3e07d2326f20ada87596734faeffc95f8f8490e780e8e95e` | `a81ff8bf5f32c06f3e07d2326f20ada87596734faeffc95f8f8490e780e8e95e` | unmodified |
| `features/cofit.py` | `abd778fe41bf3e194d06bb6111e12e565703b662f35671060258f69151107c77` | `4654d73c36af3c166b059e32042a9fee2be0f434c721545b0b8b068afd8483fa` | **patched** (path roots) |
| `features/genes.py` | `f374ed01f25b5ae4d7cd213cf5940007aba4163e548977af15fc7e88c44ccd96` | `f374ed01f25b5ae4d7cd213cf5940007aba4163e548977af15fc7e88c44ccd96` | unmodified |
| `features/pairs.py` | `f2e09e13f55545ce6f26699c143556d31ebed9a17993f7db0a9ca147830da2f2` | `6538ffc24458af97b342c13079fab92cd31c7a4f8c10423a28e3ea86331121d1` | **patched** (path roots) |
| `models/__init__.py` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855` | unmodified |
| `models/baselines.py` | `ee06051c96c5c7d6a328ed6227b87c350ad9982b69cf7011b7c62fe77c4b5b3d` | `ee06051c96c5c7d6a328ed6227b87c350ad9982b69cf7011b7c62fe77c4b5b3d` | unmodified |
| `models/combiner.py` | `33dbd8927d513520ae9d71ad11a22fa57d809a16b2d82b1fbc6f133daae0efd6` | `33dbd8927d513520ae9d71ad11a22fa57d809a16b2d82b1fbc6f133daae0efd6` | unmodified |
| `models/completion.py` | `61fce29235df19de9f6d88c72f2473647d2dd4dc12bfed32405cacd2aeff253f` | `61fce29235df19de9f6d88c72f2473647d2dd4dc12bfed32405cacd2aeff253f` | unmodified |
| `models/external.py` | `342ee266bdccd0c296653fcb4d7a39048125ae56e2fe3ffb473b4a5d8ccd3fd6` | `342ee266bdccd0c296653fcb4d7a39048125ae56e2fe3ffb473b4a5d8ccd3fd6` | unmodified |
| `models/foldin2.py` | `f4fd0f160ff03e65e27563698d2b31fa39dcef3d4938eb2876301a2422f80a43` | `f4fd0f160ff03e65e27563698d2b31fa39dcef3d4938eb2876301a2422f80a43` | unmodified |
| `models/m4d.py` | `ea0ca6c22e6c56595222fe55cc7187bc23ea5c5ea580c50432713badf2c09197` | `ea0ca6c22e6c56595222fe55cc7187bc23ea5c5ea580c50432713badf2c09197` | unmodified |
| `models/propensity.py` | `bb70b2db3966fbc2df712a5b007a8e266baf665d0ece645d82d4f6256603f3ae` | `bb70b2db3966fbc2df712a5b007a8e266baf665d0ece645d82d4f6256603f3ae` | unmodified |
| `models/registry.py` | `32c2fc95b1c431fc3b84f9b0dda08f2c2cb7ca9f3e4838543477a497f1ab0f56` | `32c2fc95b1c431fc3b84f9b0dda08f2c2cb7ca9f3e4838543477a497f1ab0f56` | unmodified |
| `models/tier2.py` | `d439536ec901d0953d206cc5bc546de2addd1adc6ab6b3ecbc988e4c7c037b89` | `4c52a5018212145a9e168006fb8f9112d44523cf1f0d482834f7772876a8c20a` | **patched** (path roots) |
| `models/view.py` | `bc670c91e299faa65061241fa1d5e36e58dfe4136a5b2aa4ef789db278d5fae1` | `bc670c91e299faa65061241fa1d5e36e58dfe4136a5b2aa4ef789db278d5fae1` | unmodified |
| `predict.py` | `0d2770d4bedc45d95facff126a76bec27ce4cce30610c9f90bab2a398b976b05` | `6dd97159a1f728884f2abd001dc2c3ea1d172ae63934a7d4c23151c93d9efe9e` | **patched** (path roots) |
| `runner.py` | `c4b637d386f907accd496b29858666125c61020b916e36322a1bddf6e9b2da40` | `589cc465dc425c9e066824065b5eeb8b44da441578b0c4f650510bce2e3d31bd` | **patched** (path roots) |
