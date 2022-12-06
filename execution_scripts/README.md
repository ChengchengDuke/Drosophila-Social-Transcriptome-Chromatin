# Execution workdir
- README: nextflow, or bash working directory (i.e. where a pipeline is executed), description of files
- Workflow control scripts. Can refer to ../data_raw/ or ../data_clean/ as input, and ../scripts/ for custom scripts. Output should be data_clean
- Workflow control configuration files.
- Out should be in data_clean
- Execution log is a directory where to store nextflow execution reports / or a log of history bash commands

| File name | Description | Location |
| --- |--- |--- |
| pipeline.sh | A pipeline | PBS cluster |


