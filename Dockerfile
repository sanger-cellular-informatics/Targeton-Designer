FROM python:3.8.0

# Creat variables for a docker user with numerical id and named user and group
ARG GROUP_NAME=pd_group
ARG GROUP_ID=1001
ARG USER_NAME=pd_user
ARG USER_ID=1001

RUN apt-get update && \
    apt-get -y install build-essential && \
    apt-get -y install bedtools && \
    apt-get -y install curl && \
    apt-get -y install python3-setuptools

# Add new user to a new group
RUN groupadd $GROUP_NAME && useradd -g $GROUP_NAME $USER_NAME

WORKDIR /targeton-designer

# Give access to the work directory to the new user
COPY --chown=$USER_NAME:$GROUP_NAME . /targeton-designer

RUN pip install --no-cache-dir --upgrade pip setuptools wheel && \
    pip install --no-cache-dir -r requirements.txt && \
    pip install --no-cache-dir -r sge-primer-scoring/requirements.txt

# Switch to the new user inside docker container
USER $USER_NAME
