VENV = venv
PYTHON = $(VENV)/bin/python
PIP = $(VENV)/bin/pip

ifeq ($(PREFIX),)
    PREFIX := /usr/local
endif

APP = $(PREFIX)/src/app


init: 
	git config core.hooksPath .githooks

install: 
	@echo "Installing..."
	sudo apt-get update
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
	@echo "Installing python3.8-dev..."
	@sudo apt-get -y install python3.8-dev
	@sudo update-alternatives --install /usr/bin/python python /usr/bin/python3.8 2
	@sudo update-alternatives --config python

install-libglib2.0-dev: 
	@echo "Installing libglib2.0-dev..."
	sudo apt-get -y install libglib2.0-dev 

install-autoconf: 
	@echo "Installing autoconf..."
	sudo apt install autoconf

install-ipcress: 
	@echo "Installing ipcress..."
	git clone https://github.com/nathanweeks/exonerate.git $(APP)/
	cd $(APP)/exonerate \
		&& autoreconf -fi \
		&& ./configure \
		&& make -j 4 \
		&& make install
	rm -rf $(APP)/exonerate

setup-venv:
	echo "running venv"
	python -m venv venv
	./venv/bin/pip install -U pip wheel setuptools 
	./venv/bin/pip install -r requirements.txt
	./venv/bin/pip install -r sge-primer-scoring/requirements.txt

test:
	python -m unittest

clean: 
	rm -rf __pycache__
	rm -rf venv