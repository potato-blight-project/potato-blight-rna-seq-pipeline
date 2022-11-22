#!/bin/bash

DATE=`date +"%Y%m%dT%H%M%S"`
ARCHIVE="archive/${DATE}"

if [ ! -d "${ARCHIVE}" ]; then
	mkdir -p "${ARCHIVE}"
fi

mv slurm* "${ARCHIVE}"
