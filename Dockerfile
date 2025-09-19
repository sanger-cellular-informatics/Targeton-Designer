FROM python:3.10-slim-bullseye AS build-stage

# Creat variables for a docker user with numerical id and named user and group
ARG GROUP_NAME=pd_group
ARG GROUP_ID=1001
ARG USER_NAME=pd_user
ARG USER_ID=1001

RUN apt-get update && apt-get install -y \
    build-essential \
    python3-dev \
    gfortran \
    bedtools \
    curl \
 && rm -rf /var/lib/apt/lists/*

# Add new user to a new group
RUN groupadd $GROUP_NAME && useradd -g $GROUP_NAME $USER_NAME

WORKDIR /targeton-designer

# Give access to the work directory to the new user
COPY --chown=$USER_NAME:$GROUP_NAME . /targeton-designer

RUN pip install --upgrade pip setuptools wheel
RUN pip install --no-cache-dir keeper-secrets-manager-cli
RUN pip install -r requirements.txt

# Switch to the new user inside docker container
USER $USER_NAME


FROM python:3.10-slim-bullseye AS dev-stage

# Re-create the user and group from the builder stage
ARG GROUP_NAME=pd_group
ARG GROUP_ID=1001
ARG USER_NAME=pd_user
ARG USER_ID=1001

RUN groupadd $GROUP_NAME && useradd -g $GROUP_NAME $USER_NAME

WORKDIR /targeton-designer

# Copy only the necessary files and installed dependencies from the builder stage
COPY --from=build-stage --chown=$USER_NAME:$GROUP_NAME /targeton-designer /targeton-designer
COPY --from=build-stage --chown=$USER_NAME:$GROUP_NAME /usr/local/lib/python3.10/site-packages /usr/local/lib/python3.10/site-packages

RUN pip install --no-cache-dir keeper-secrets-manager-cli
RUN pip install s3cmd


# Switch to the new user
USER $USER_NAME

# Set home and work directories
ENV HOME=/home/$USER_NAME
WORKDIR /home/$USER_NAME