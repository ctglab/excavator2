#!/bin/bash
set -e

# Activate the micromamba environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate excavator2

# Execute the command passed to this script (the Docker CMD)
exec "$@"
