#!/bin/bash

echo "Building image from"

apptainer build ../paws_gc.sif paws.def

echo "Done"
