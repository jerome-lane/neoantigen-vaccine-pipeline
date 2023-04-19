!#/bin/bash

conda create -n benchmark python=3.10
eval "$(conda 'shell.bash' 'hook' 2> /dev/null)"
conda activate benchmark
pip install -r requirements.txt
