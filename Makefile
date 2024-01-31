.PHONY: all lint test test-cov install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	coverage run -m pytest
	coverage xml

install: all
	$(PYTHON) setup.py install

dev: all
	pip install -e .

prep-dev-container: all
	QIIME_VERSION=2023.9
	conda install mamba -qy -n base -c conda-forge
	mamba install -n "qiime2-amplicon-$QIIME_VERSION" -qy -c conda-forge -c bioconda -c defaults flake8 coverage wget pytest-xdist autopep8
	/opt/conda/envs/qiime2-amplicon-$QIIME_VERSION/bin/pip install -q https://github.com/qiime2/q2lint/archive/master.zip
	/opt/conda/envs/qiime2-amplicon-$QIIME_VERSION/bin/pip install -e .

clean: distclean

distclean: ;
