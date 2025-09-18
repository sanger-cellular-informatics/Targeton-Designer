FROM python:3.10-slim-bullseye

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
RUN pip install numpy==1.26.4
RUN pip install pandas==2.0.3
RUN pip install --no-cache-dir keeper-secrets-manager-cli
RUN pip install s3cmd
RUN pip install -r requirements.txt

# Switch to the new user inside docker container
USER $USER_NAME

# Create Home Directory to utilize s3cmd config
ENV HOME=/home/$USER_NAME
WORKDIR /home/$USER_NAME