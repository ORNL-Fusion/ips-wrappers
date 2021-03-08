#!/bin/bash
rm -r simulation_setup
rm -rf work
rm -r simulation_results
rm -r simulation_log
rm -r ${TOKAMAK_ID}${SHOT_NUMBER}.${TIME_ID}
rm resource_usage
rm checklist.conf
rm *.zip
rm PORTAL_RUNID
rm *_log*
rm -r www
rm dask_preload.py
