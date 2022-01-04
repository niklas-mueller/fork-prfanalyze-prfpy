# PRFPY - Validation Framework

This repository is an extension of ![prfanalyze](https://github.com/vistalab/PRFmodel) implementing ![prfpy](https://github.com/VU-Cog-Sci/prfpy) into the validation framework. This way, we ensure our toolbox is consistent with the other toolboxes already validated in ![this paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007924).

---

To easify the following process, create a link to the directory that contains the data and should be used for the output:

``
export basedir=/path/to/folder
``

## Synthesize

For a detailed describing on how the synthesis of new data work please see ![prfsynth](https://github.com/vistalab/PRFmodel/wiki).

The main steps are:

1. Create a configuration file

``
docker run --rm -it \
        -v empty:/flywheel/v0/input/config.json \
        -v $basedir:/flywheel/v0/output \
        garikoitz/prfsynth
``

2. Rename the configuration file according to the subject number and session (date). Also fill in these numbers into the config file itself.

``
mv $basedir/prfsynth-config-defaults.json $basedir/prfsynth_sub-001_sess-20200320.json
``

3. Run the synthesis

``
docker run --rm -it \
        -v $basedir/prfsynth_sub-001_sess-20200320.json:/flywheel/v0/input/config.json \
        -v $basedir:/flywheel/v0/output \
           garikoitz/prfsynth
``

## Analysis

The integration into the validation framework is yet in progress. The analysis can be called when cloning the ![PRFModel](https://github.com/vistalab/PRFmodel) and executing inside the top level folder of prfanalyze-prfpy:

``
/PATH/TO/PRFmodel/gear/prfanalyze/run_prfanalyze.sh prfpy $basedir default_config.json
``

## Report

1. Create a report configuration file

``
docker run --rm -it \
        -v empty:/flywheel/v0/input/config.json \
        -v $basedir:/flywheel/v0/output \
        garikoitz/prfreport
``

2. Rename the file and edit the config file the above subject and session info:

``
mv $basedir/prfreport-configuration-defaults.json $basedir/prfreport-config_sub-001_sess-20200320.json
``

3. Run the creation of the report

``
docker run --rm -it \
        -v $basedir/prfreport-config_sub-001_sess-20200320.json:/flywheel/v0/input/config.json \
        -v $basedir:/flywheel/v0/output \
        garikoitz/prfreport
``

