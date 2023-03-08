.ONESHELL:
SHELL := /bin/bash

VENV = venv
PYTHON = $(VENV)/bin/python
PIP = $(VENV)/bin/pip

ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

APP = $(PREFIX)/src/app


init: 
	git config core.hooksPath .githooks
	chmod +x .githooks/*

install: 
	@echo "Installing..."
	@sudo apt-get update
	@if [ "$(shell which build-essential)" = "" ]; then \
        $(MAKE) install-build-essential; \
    fi
	@if [ "$(shell which bedtools)" = "" ]; then \
        $(MAKE) install-bedtools; \
    fi
	@if [ "$(shell which python3.8-dev)" = "" ]; then \
        $(MAKE) install-python3.8-dev; \
    fi
	@if [ "$(shell which python3.8-venv)" = "" ]; then \
        $(MAKE) install-python3.8-venv; \
    fi
	@if [ "$(shell which libglib2.0-dev)" = "" ]; then \
        $(MAKE) install-libglib2.0-dev; \
    fi
	@if [ "$(shell which autoconf)" = "" ]; then \
        $(MAKE) install-autoconf; \
    fi
	@if [ "$(shell which ipcress)" = "" ]; then \
        $(MAKE) install-ipcress; \
    fi

install-bedtools:
	@echo "Installing bedtools..."
	@sudo apt-get -y install bedtools

install-build-essential:
	@echo "Installing build-essential..."
	@sudo apt-get -y install build-essential

install-python3.8-venv:
	@echo "Installing python3.8-venv..."
	@sudo apt-get -y install python3.8-venv

install-python3.8-dev: 
	@if [ "$(shell which python3)" = "" ]; then \
		# Python3 not installed.
		@echo "Installing python3.8-dev..."
		@sudo apt-get -y install python3.8-dev
	else
		PYTHONPATH = which python
		ver=$(python3 -V 2>&1 | sed 's/.* \([0-9]\).\([0-9]\).*/\1\2/')
		@if [ "${ver}" -ge "38" ]; then
			@echo "Using existing python3.8..."
			PYTHONPATH38 = which python3
		else
			@echo "Installing python3.8-dev..."
			@sudo apt-get -y install python3.8-dev
			PYTHONPATH38 = which python3.8
		fi
		@sudo update-alternatives --install ${PYTHONPATH} python ${PYTHONPATH38} 2 
		@sudo update-alternatives --config python 
	fi

install-libglib2.0-dev: 
	@echo "Installing libglib2.0-dev..."
	@sudo apt-get -y install libglib2.0-dev 

install-autoconf: 
	@echo "Installing autoconf..."
	@sudo apt-get -y install autoconf libtool

install-ipcress: 
	@echo "Installing ipcress..."
	@git clone https://github.com/nathanweeks/exonerate.git $(APP)/
	@cd $(APP)/exonerate \
		&& autoreconf -fi \
		&& ./configure \
		&& make -j 4 \
		&& make install
	@rm -rf $(APP)/exonerate

venv/bin/activate:
	@python -m venv venv

setup-venv: venv/requirements_run

venv/requirements_run: venv/bin/activate requirements.txt
	@./venv/bin/pip install -U pip wheel setuptools 
	@./venv/bin/pip install -r requirements.txt
	@./venv/bin/pip install -r sge-primer-scoring/requirements.txt
	@touch venv/requirements_run

clean-venv/requirements_run:
	@rm -f venv/requirements_run
	
test: setup-venv
	@. venv/bin/activate \
	&& python -m unittest

clean: 
	@rm -rf __pycache__
	@rm -rf venv