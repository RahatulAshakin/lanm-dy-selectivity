# Dataset Audit

## Configuration
- project_name: `lanm-dy-selectivity`
- target_metal: `Dy`
- competitors: `Nd, Y, Al, Fe`
- first_shell_cutoff_A: `3.2`
- second_sphere_cutoff_A: `6.0`

## Normalization
- normalized_file_count: `11`
- normalized_paths: `data/raw/local_bundle/8dq2.txt, data/raw/local_bundle/8dq2_atoms.csv, data/raw/local_bundle/8fns.txt, data/raw/local_bundle/8fns_atoms.csv, data/raw/local_bundle/project_description.docx, data/raw/local_bundle/readme.txt, data/raw/local_bundle/lanmodulin_lanthanide_structures.csv, data/raw/local_bundle/lanmodulin_sequences.csv, data/raw/local_bundle/lanthanide_dataset_manifest.json, data/raw/local_bundle/lanthanide_elements_basic_properties.csv, data/raw/local_bundle/lanthanide_repository_download_links.txt`

## Project Brief Extraction
- status: `extracted`
- source_path: `data/incoming/Project description .docx`
- output_path: `docs/project_brief_extracted.md`
- note: Extracted from incoming DOCX.

## Inventory Summary
- file_count: `11`
- total_size_bytes: `857800`

| relative_path | file_type | rows | size_bytes | chains | residue_ranges | hetatm_count |
| --- | --- | ---: | ---: | --- | --- | ---: |
| data/raw/local_bundle/8dq2.txt | text | 4613 | 373653 | - | - | 0 |
| data/raw/local_bundle/8dq2_atoms.csv | csv | 3641 | 218016 | A;B;C;D | A:24-374;B:24-371;C:24-362;D:24-356 | 305 |
| data/raw/local_bundle/8fns.txt | text | 2191 | 177471 | - | - | 0 |
| data/raw/local_bundle/8fns_atoms.csv | csv | 964 | 57165 | A | A:29-467 | 171 |
| data/raw/local_bundle/lanmodulin_lanthanide_structures.csv | csv | 4 | 910 | - | - | 0 |
| data/raw/local_bundle/lanmodulin_sequences.csv | csv | 3 | 1242 | - | - | 0 |
| data/raw/local_bundle/lanthanide_dataset_manifest.json | json | 4 | 1340 | - | - | 0 |
| data/raw/local_bundle/lanthanide_elements_basic_properties.csv | csv | 15 | 590 | - | - | 0 |
| data/raw/local_bundle/lanthanide_repository_download_links.txt | text | 27 | 1415 | - | - | 0 |
| data/raw/local_bundle/project_description.docx | docx | 0 | 25081 | - | - | 0 |
| data/raw/local_bundle/readme.txt | text | 18 | 917 | - | - | 0 |
