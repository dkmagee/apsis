
#------------------------------------------------------------------------------
# User Settings
#------------------------------------------------------------------------------
# path to the Apsis root directory (where you installed Apsis)
export ACS_PIPE=${HOME}/apsis

# the root directory that contains the input datasets
export INGEST=${HOME}/apsis_runs/ingest

# the root directory of the Apsis output data products
export DATASETS=${HOME}/apsis_runs/datasets


#------------------------------------------------------------------------------
# Fixed Settings
#------------------------------------------------------------------------------
export PATH=${ACS_PIPE}/bin:${PATH}
export NUMERIX=numpy
export PIPELINE=${ACS_PIPE}/reffiles
export DUST_DIR=${PIPELINE}
export jref=${PIPELINE}/idctab/     # trailing slash required
export APSIS_PYTHONPATH=${ACS_PIPE}/python
export APSIS_PYTHONPATH=${APSIS_PYTHONPATH}:${ACS_PIPE}/python/apsis
export APSIS_PYTHONPATH=${APSIS_PYTHONPATH}:${ACS_PIPE}/python/utils
if [ ${PYTHONPATH} ]; then
   export PYTHONPATH=${APSIS_PYTHONPATH}:{PYTHONPATH}
else
   export PYTHONPATH=${APSIS_PYTHONPATH}
endif
unexport APSIS_PYTHONPATH

