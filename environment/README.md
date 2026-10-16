# Environment records

`R_packages.csv` contains package versions observed in the available R 4.5.1 library during code preparation. It is not an installation script or a claim that every historical analysis used those exact versions. Consult preserved original session records where available. `lmom` is required by the SRI sensitivity script, which also calculates Gamma L-moment comparisons; the manuscript selects Gamma MLE. Its locally observed version was found in the hydrological project library.

`requirements.txt` records locally available versions of the three Python packages used by the supporting extraction and PDF-layout scripts. The MBPD script uses only Python's standard library and calls the separately installed BLAST+ executables.

The QIIME 2 workflow targets the 2026.4 amplicon environment. Use the matching DADA2, feature-classifier, alignment and phylogeny plugins. Reference SILVA and MBPD database versions and original downloads must accompany a reproducible run.

R plotting requires a working graphics device and suitable fonts. sf additionally needs compatible geospatial libraries. Some source scripts request Windows locale names; on other operating systems these calls may warn. Use a UTF-8 locale and retain UTF-8 filenames. Do not copy the original project's installed R packages into this repository.

Computational demand varies substantially. The original iCAMP settings use 16 workers, a 50-GB allowance and 1,000 randomizations. No runtime estimate for a new machine has been measured.
