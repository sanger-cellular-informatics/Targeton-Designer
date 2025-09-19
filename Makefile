.EXPORT_ALL_VARIABLES:
.ONESHELL:
SHELL := /bin/bash

VENV = venv
PYTHON = $(VENV)/bin/python

PYTHONPATH ?= /usr/bin/python
PYTHONPATH38 ?= /usr/bin/python3.8

PIP = $(VENV)/bin/pip
ENVIRONMENTAL_VARIABLE_FILE := .env

-include ${ENVIRONMENTAL_VARIABLE_FILE}

PREFIX ?= /usr/local

APP = $(PREFIX)/src/app
export DOCKER_ENV ?= test
	
# Docker
DOCKER_NAME ?= primer_designer
DOCKER_TAG ?=${DOCKER_ENV}
DOCKER_REPO ?=local
DOCKER_PORT ?=8081
DOCKER_IMAGE_NAME ?= ${DOCKER_REPO}/${DOCKER_NAME}:${DOCKER_TAG}

init: check-make
	git config core.hooksPath .githooks
	chmod +x .githooks/*


install:
	@echo "Installing..."
	@apt-get update 
	
	@echo "Installing build-essential..."
	@apt-get -y install build-essential
	
	@echo "Installing bedtools..."
	@apt-get -y install bedtools

	@echo "Installing curl..."
	@apt-get -y install curl

	$(MAKE) install-python3.8-dev
	$(MAKE) install-python-setuptools


install-python3.8-dev:
	@echo "Installing python3.8-dev..."
	@apt-get install -y python3.8-dev
	
	@echo "Install Python Virtual Env..."
	@apt-get install -y python3.8-venv
	
	@echo "Updating python alternatives list..."
	@update-alternatives --install $(PYTHONPATH) python $(PYTHONPATH38) 2
	@update-alternatives --config python
	
	@echo "Installing Python3 pip..."
	@apt-get install -y python3-pip

	@echo "Updated Python version:"
	@python3 --version
	@python --version



install-python-setuptools:
	@echo "Installing Python Setup Tools..."
	@apt-get install python3-setuptools



install-docker:
	@echo "Installing docker..."
	@curl -fsSL https://get.docker.com -o get-docker.sh
	@sh get-docker.sh
	@groupadd docker
	@usermod -aG docker $$USER
	@newgrp docker


check-make:
	@MAKE_VERSION=$$(make --version | grep '^GNU Make' | sed 's/^.* //g');
	echo "Detected make version: $$MAKE_VERSION";
	if (( $$(echo "$$MAKE_VERSION < 3.82" | bc -l) )); then
		echo "make version = $$MAKE_VERSION, minimum version 3.82 required for multiline.";
		sudo apt install make;
	fi

install-autopep8: venv/bin/activate
	@echo "Installing autopep8..."
	@./venv/bin/pip install autopep8


create-venv: 
	@echo "Creating Virtual Env..."
	@python -m venv venv


check-venv:
	@if [ ! -d "venv/bin/" ]; then \
		echo "Creating Virtual Env..."; \
		python -m venv venv; \
	fi

setup-venv: create-venv check-venv 
	@./venv/bin/pip install --upgrade pip
	@./venv/bin/pip install --upgrade pip setuptools wheel
	@./venv/bin/pip install -r requirements.txt
	@./venv/bin/pip install -r sge-primer-scoring/requirements.txt
	@echo "Python requirements installed."

	
test: setup-venv
	@. venv/bin/activate
	pip list
	python -m unittest discover --start-directory ./tests --top-level-directory .

download-kmers:
	bash download_kmer_lists.sh

build-docker:
	@ver=$$(docker version --format '{{.Server.Version}}' 2>&1 | sed -E 's/([0-9]+).*/\1/')
	@echo Docker version $$ver
	if [ "$$ver" -lt 23 ]; then
		echo "Warning Docker engine version <23, changing build to buildx."
		docker buildx install
		export DOCKER_BUILDKIT=1
	fi
	echo docker repo = ${DOCKER_REPO}
	echo docker image = ${DOCKER_IMAGE_NAME}
	if [ "$(docker images -q ${DOCKER_IMAGE_NAME} 2> /dev/null)" != "" ]; then
		@echo "docker image already exists. ${DOCKER_IMAGE_NAME}"
	else
		@echo "Building docker image..."
		@docker build --pull -t "${DOCKER_IMAGE_NAME}" .;
		if [[ ${DOCKER_REPO} != "local" ]]; then
			@docker push "${DOCKER_IMAGE_NAME}" 
		fi
	fi

build-docker-test: build-docker
	docker build --cache-from="${DOCKER_IMAGE_NAME}" -t "${DOCKER_IMAGE_NAME}" --target unittest .;

run-docker: build-docker
	@echo "Running Docker image =  $(DOCKER_IMAGE_NAME)"
	docker run --name "${DOCKER_NAME}" -p ${DOCKER_PORT}:${DOCKER_PORT} -t "${DOCKER_IMAGE_NAME}"

run-docker-test: build-docker
	@echo "Running Docker image =  $(DOCKER_IMAGE_NAME)"
	@docker run --name "${DOCKER_NAME}" -p ${DOCKER_PORT}:${DOCKER_PORT} -t "${DOCKER_IMAGE_NAME}" make test

run-docker-interactive: build-docker
	@echo "Running Docker image =  $(DOCKER_IMAGE_NAME)"
	@docker run -i --name "${DOCKER_NAME}" -t "${DOCKER_IMAGE_NAME}" bash

connect-docker-interactive:
	@docker exec -it ${DOCKER_NAME} bash

clean-docker-containers:
	@docker rm -f $$(docker ps -a -q)

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

auto-lint-tests: install-autopep8
	@echo "Linting tests..."
	@. venv/bin/activate
	@python -m autopep8 -r -i tests/

auto-lint-src: install-autopep8
	@echo "Linting src..."
	@. venv/bin/activate
	@python -m autopep8 -r -i src/

clean:clean-docker
	@rm -rf __pycache__
	@rm -rf venv
