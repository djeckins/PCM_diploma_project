# setup / test / lint / run
# Override either variable if conda lives elsewhere or the environment has
# another name, e.g. make test CONDA_SH=$HOME/miniconda3/etc/profile.d/conda.sh
CONDA_SH  ?= $(HOME)/miniforge3/etc/profile.d/conda.sh
CONDA_ENV ?= chem
CONDA_RUN := source $(CONDA_SH) && conda activate $(CONDA_ENV) &&

setup:
	$(CONDA_RUN) pip install -e . --no-deps

test:
	$(CONDA_RUN) python -m pytest

lint:
	$(CONDA_RUN) python -m ruff check src tests tools

# make run CONFIG=configs/mthk.yaml STEP=all
run:
	$(CONDA_RUN) python -m pcm2.run $(STEP) --config $(CONFIG)

.PHONY: setup test lint run
