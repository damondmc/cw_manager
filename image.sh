#!/bin/bash

echo "Building image from"

apptainer build ../paws_v3.sif paws.def

echo "Done"
