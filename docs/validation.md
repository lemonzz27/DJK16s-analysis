# Validation record

Checked on 16 October 2026.

- All 74 R source files parsed successfully; no explanatory R comments remain. Interpreter directives are retained.
- The 73 selected analysis R scripts have equivalent parsed calculations after documented project-path substitutions. The comparison normalizes redundant parentheses around callable expressions. All four code-fragment slices used by final figure helpers remain parseable and equivalent.
- All five Python files, including the new launcher, compile successfully. No Python comments or docstrings remain except the interpreter directive in MBPD.
- Both Bash scripts pass `bash -n`; this does not verify QIIME 2 execution or third-party program availability.
- The complete iCAMP call and five MBPD analysis functions were compared with the supplied sources and are unchanged. The MBPD FASTA reader accepts valid IDs and rejects blank/duplicate IDs; its nonempty-output guard was exercised.
- All 80 selected original source files retain their original SHA-256 checksums. Explicit local source dependencies are present in the code collection.
- The existing MBPD filtered-hit table contains 57,964 records and 2,755 ASVs, with minimum pident 97, minimum qcovs 95 and maximum E-value 1.2e-150. It does not recover the original command's HSP-level coverage parameter.

No full QIIME 2, BLAST, iCAMP, TITAN, random-forest or complete manuscript rerun was performed. Raw/reference data and cached result availability were not validated end to end. Numerical/figure equivalence of a new complete run therefore remains untested. Static code checks are not a substitute for that validation.

`provenance/release_sha256.tsv` records the delivered file hashes. `environment/RF_original_sessionInfo.txt` preserves the available original random-forest session record, while `R_packages.csv` is the separately observed local environment.
