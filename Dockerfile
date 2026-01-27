# Stage 1: Development Image
FROM python:3.10.12-slim-bullseye AS development

ARG GROUP_NAME=pd_group
ARG GROUP_ID=1001
ARG USER_NAME=pd_user
ARG USER_ID=1001

RUN apt-get update && apt-get install -y \
    build-essential \
    python3-dev \
    gfortran \
    bedtools \
    procps \
    curl \
 && rm -rf /var/lib/apt/lists/*

RUN groupadd -g $GROUP_ID $GROUP_NAME && \
    useradd -u $USER_ID -g $GROUP_NAME -m -d /home/$USER_NAME $USER_NAME

WORKDIR /targeton-designer

COPY --chown=$USER_NAME:$GROUP_NAME . /targeton-designer

RUN pip install --upgrade pip setuptools wheel
RUN pip install --no-cache-dir keeper-secrets-manager-cli
RUN pip install s3cmd
RUN pip install -r requirements.txt

USER $USER_NAME
ENV HOME=/home/$USER_NAME
WORKDIR /home/$USER_NAME

# Stage 2: Production Image
FROM python:3.10.12-slim-bullseye AS production

ARG GROUP_NAME=pd_group
ARG GROUP_ID=1001
ARG USER_NAME=pd_user
ARG USER_ID=1001

RUN apt-get update && apt-get install -y --no-install-recommends \
    bedtools \
    procps \
 && rm -rf /var/lib/apt/lists/*

RUN groupadd -g $GROUP_ID $GROUP_NAME && \
    useradd -u $USER_ID -g $GROUP_NAME -m -d /home/$USER_NAME -s /bin/bash $USER_NAME

WORKDIR /targeton-designer

COPY --from=development /usr/local/lib/python3.10/site-packages /usr/local/lib/python3.10/site-packages
COPY --from=development /usr/local/bin /usr/local/bin

COPY --chown=$USER_NAME:$GROUP_NAME . /targeton-designer

USER $USER_NAME

ENV HOME=/home/$USER_NAME \
    PATH=/usr/local/bin:$PATH \
    PYTHONUNBUFFERED=1

WORKDIR /targeton-designer
