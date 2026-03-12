# Lanthanide datasets package

This package contains a small, ready-to-use local bundle of lanthanide-related datasets and manifests assembled from public sources.

Included files
- lanthanide_elements_basic_properties.csv
- lanmodulin_lanthanide_structures.csv
- lanthanide_dataset_manifest.json

Notes
- I could assemble and save structured CSV manifests locally in this environment.
- I could not reliably fetch every external binary file directly from the source websites from inside this session.
- For the external datasets that were not fetched as full raw downloads, the manifest includes direct public access locations.

Recommended next use
- Use lanmodulin_lanthanide_structures.csv to batch-download CIF files from RCSB.
- Use lanthanide_elements_basic_properties.csv as a base descriptor table for lanthanide-aware modeling.
- Use the PRIDE and NIST manifest entries if you want the larger raw datasets.
