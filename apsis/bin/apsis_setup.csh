
#------------------------------------------------------------------------------
# User Settings
#------------------------------------------------------------------------------
# path to the Apsis root directory (where you installed Apsis)
setenv ACS_PIPE ${HOME}/apsis

# the root directory that contains the input datasets
setenv INGEST ${HOME}/apsis_runs/ingest

# the root directory of the Apsis output data products
setenv DATASETS ${HOME}/apsis_runs/datasets


#------------------------------------------------------------------------------
# Fixed Settings
#------------------------------------------------------------------------------
setenv PATH ${ACS_PIPE}/bin:${PATH}
setenv NUMERIX numpy
setenv PIPELINE ${ACS_PIPE}/reffiles
setenv DUST_DIR ${PIPELINE}
setenv jref ${PIPELINE}/idctab/     # trailing slash required
setenv APSIS_PYTHONPATH ${ACS_PIPE}/python
setenv APSIS_PYTHONPATH ${APSIS_PYTHONPATH}:${ACS_PIPE}/python/apsis
setenv APSIS_PYTHONPATH ${APSIS_PYTHONPATH}:${ACS_PIPE}/python/utils
if (${?PYTHONPATH}) then
   setenv PYTHONPATH ${APSIS_PYTHONPATH}:${PYTHONPATH}
else
   setenv PYTHONPATH ${APSIS_PYTHONPATH}
endif
unsetenv APSIS_PYTHONPATH

