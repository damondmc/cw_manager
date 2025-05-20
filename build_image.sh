#!/bin/bash

echo "Building image from"

apptainer build cw_manager_gc.sif cw_image.def

echo "Done"
