# Version-matched VitalDB source catalogue

`vitaldb_v1.0.0_track_index.tsv.gz` contains the filenames, recording bounds,
available track names, and file sizes read from the headers of all 6,388 public
`.vital` files in the PhysioNet VitalDB v1.0.0 release.

It contains no waveform samples, clinical values, or derived analysis results.
The downloader uses it with the official v1.0.0 `clinical_data.csv` to identify
the 1,577 records that must be downloaded before the volatile-anaesthetic
prescreen. This avoids downloading the complete 95.4 GB archive while retaining
the original result-independent cohort construction.

Source release: https://physionet.org/content/vitaldb/1.0.0/  
DOI: https://doi.org/10.13026/czw8-9p62

The catalogue is derived from the cited VitalDB release and remains subject to
its Creative Commons Attribution 4.0 licence. The repository software is
separately released under the MIT License.
