# Performance reproduction programs

These are byte-for-byte copies of the recorded benchmark and additive analysis
programs. The original experiment sources/identities and the separately hashed
postprocessor must remain distinct. See ../PERFORMANCE.md for supported verdicts.

For an exact historical rerun, use the companion upstream-delivery archive,
which supplies the frozen correctness prerequisites, gate, calibration/raw data,
monitors, original environment records and selection ledger. Follow its
performance-reproduction.md in a copy of that extracted evidence directory.
The scripts require those source-bound prerequisites and must not be described
as a self-contained clean benchmark without them. The archive preserves raw
records with original machine paths; portable repository review copies omit
installation locations without changing source hashes or physical findings.

Python3.11+ is required for the additive retry-aware analysis. Julia1.12.7 and
the exact recorded historical KS05d source are required for the original timing
experiment; the normal adapter installer now pins finalc905. All failed attempts
and every slow call remain visible. J2 wall-time performance is inconclusive.
