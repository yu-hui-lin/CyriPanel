# Merge integration source review

Validated example source: `bd23768917ff889e68abc34fdb1bce8ebdaea3aa`.

Incorporated main source: `7a1b0ea9f502838df5a4796b9e7bbcac79f6f079`.

All files in the historical source manifest were compared with the validated snapshot. Only full-line comments changed in the two runtime files below. Their non-comment bytes are identical. The example runner changes only the expected R SHA256 and its diagnostic label.

`example/validation/source-sha256.json` identifies the earlier tested source, not the current merged source. Frozen outputs and validation records remain unchanged. This note records source review, not a new BAM analysis.

| File | Original SHA256 | Merged SHA256 |
| --- | --- | --- |
| `depth_calling/run_CNVPanelizer.R` | `8a864517cb066fcd45af6c21ae2c9a187fc57a8cb270ebbe1a097f2a27a308b4` | `d41cf4e1133fb4849020d8d83443bbc24b39ff6f24497be6c0d5b683b7412597` |
| `depth_calling/copy_number_call.py` | `81857a8877c7a5e0e984d84231df939fca40bc1fdeb9f26011e60c1d12490dc7` | `09ee2483044fc10d271d5356893fbee17683d98ef8ea8b3e1b88d7aad6da6b6c` |
| `example/run_example.py` | `4d32d28f6a0ad161b189beacb8f310f628fbce1726318f9c984f0bcc6b4f91dc` | `1ac4064061603bf50c64b8f648e3a59cbc5a89f652d99a3f0a12409cd5824498` |
