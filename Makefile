.ONESHELL:

VENV = venv
PYTHON = $(VENV)/bin/python
PIP = $(VENV)/bin/pip
ENVIRONMENTAL_VARIABLE_FILE := .env

-include ${ENVIRONMENTAL_VARIABLE_FILE}

PREFIX ?= /usr/local

APP = $(PREFIX)/src/app
export DOCKER_ENV ?= test

MAKE_VERSION := $(shell make --version | grep '^GNU Make' | sed 's/^.* //g')
$(info "make version = ${MAKE_VERSION}, minimum version 3.82 required for multiline.")

# Docker
name = "primer-designer"
tag=$$DOCKER_ENV

BUILD_DOCKER ?= ${name}-${tag}
$(info $(BUILD_DOCKER))

init:
	git config core.hooksPath .githooks
	chmod +x .githooks/*

install:
	@echo "Installing..."
	@apt-get update
	@if [ "$(shell which build-essential)" = "" ]; then
		$(MAKE) install-build-essential;
	fi
	@if [ "$(shell which bedtools)" = "" ]; then
		$(MAKE) install-bedtools;
	fi
	@if [ "$(shell which python3.8-dev)" = "" ]; then
		$(MAKE) install-python3.8-dev;
	fi
	@if [ "$(shell which python3.8-venv)" = "" ]; then
		$(MAKE) install-python3.8-venv;
	fi
	@if [ "$(shell which libglib2.0-dev)" = "" ]; then
		$(MAKE) install-libglib2.0-dev;
	fi
	@if [ "$(shell which autoconf)" = "" ]; then
		$(MAKE) install-autoconf;		
	fi
	@if [ "$(shell which ipcress)" = "" ]; then
		$(MAKE) install-ipcress;
	fi

install-bedtools:
	@echo "Installing bedtools..."
	@apt-get -y install bedtools

install-build-essential:
	@echo "Installing build-essential..."
	@apt-get -y install build-essential git

install-python3.8-venv:
	@echo "Installing python3.8-venv..."
	@apt-get -y install python3.8-venv

install-python3.8-dev:
	@if [ "$(shell which python3)" = "" ]; then
		# Python3 not installed.
		@echo "Installing python3.8-dev..."
		@apt-get -y install python3.8-dev
	else
		@PYTHONPATH = which python
		@ver=$$(python3 -V 2>&1 | sed 's/.*([0-9]\).\([0-9]\).*/\1\2/')
		@if [ "$$ver" -ge 38 ]; then
			PYTHONPATH38 = which python3
		else
			@echo "Installing python3.8-dev..."
			@apt-get -y install python3.8-dev
			PYTHONPATH38 = which python3.8
		fi
		@update-alternatives --install $$PYTHONPATH python $$PYTHONPATH38 2 
		@update-alternatives --config python 
	fi

install-libglib2.0-dev:
	@echo "Installing libglib2.0-dev..."
	@apt-get -y install libglib2.0-dev 

install-autoconf:
	@echo "Installing autoconf..."
	@apt-get -y install autoconf libtool

install-ipcress:
	@echo "Installing ipcress..."
	@git clone https://github.com/nathanweeks/exonerate.git $(APP)/exonerate
	@cd $(APP)/exonerate \
		&& autoreconf -fi \
		&& ./configure \
		&& make -j 4 \
		&& make install
	@rm -rf $(APP)/exonerate

install-sudo:
	@echo "Installing sudo..."
	@apt-get update
	@apt-get -y install sudo

install-docker:
	@echo "Installing docker..."
	@curl -fsSL https://get.docker.com -o get-docker.sh
	@sh get-docker.sh
	@groupadd docker
	@usermod -aG docker $$USER
	@newgrp docker

install-basics:install-install-curl install-autoconf
	@apt-get -y install build-essential

venv/bin/activate:
	@python -m venv venv

setup-venv: venv/requirements_run

venv/requirements_run: venv/bin/activate requirements.txt
	@. venv/bin/activate
	@pip install -U pip wheel setuptools 
	@pip install -r requirements.txt
	@pip install -r sge-primer-scoring/requirements.txt
	@echo "Python requirements installed."
	@touch venv/requirements_run

clean-venv/requirements_run:
	@rm -f venv/requirements_run
	
test: setup-venv
	@. venv/bin/activate
	python -m unittest

$(BUILD_DOCKER): 
	@ver=$$(docker version --format '{{.Server.Version}}' 2>&1 | sed -E 's/([0-9]+).*/\1/')
	if [ "$$ver" -lt 23 ]; then
		echo Warning Docker engine version $$ver< 23, changing build to buildx.
		docker buildx install
		DOCKER_BUILDKIT=1
	fi
	@docker build --pull -t "${name}:${tag}" --target base .;
	@touch $(BUILD_DOCKER)

build-docker-test: $(BUILD_DOCKER)
	@docker build --pull -t "${name}:${tag}" --target unittest .;

run-docker: $(BUILD_DOCKER)
	@docker run --name "${name}" -p 8081:8081 -t "${name}:${tag}"

run-docker-test: build-docker-test run-docker

run-docker-interactive: $(BUILD_DOCKER)
	@docker run -i --name "${name}" -t "${name}:${tag}" bash

connect-docker-interactive: run-docker
	@docker exec -it ${name} bash

clean-docker:
	@docker builder prune -af
	@docker container prune -f
	@docker image prune -af
	@rm -f docker-touch

check-lint: 
	@echo "Running pycodestyle for src/"
	@. venv/bin/activate
	@pycodestyle --statistics -qq src || true
	@echo "Running pycodestyle for tests/"
	@pycodestyle --statistics -qq tests || true

auto-lint-tests:
	@. venv/bin/activate
	@python -m autopep8 -r -i tests/

auto-lint-src:
	@. venv/bin/activate
	@python -m autopep8 -r -i src/

clean:clean-docker
	@rm -rf __pycache__
	@rm -rf venv