#!/bin/bash

usage() {
  echo "Usage: $0 -v <volume_name> -c <container_name>"
  exit 1
}

vol_name=""
container_name=""

while getopts ":v:c:" opt; do
  case ${opt} in
    v )
      vol_name=$OPTARG
      ;;
    c )
      container_name=$OPTARG
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      usage
      ;;
    : )
      echo "Option -$OPTARG requires an argument." 1>&2
      usage
      ;;
  esac
done

# Check if the required flags were provided
if [ -z "$vol_name" ] || [ -z "$container_name" ]; then
  echo "Error: Both -v and -c flags are required."
  usage
fi

# Run the docker command with the specified volume and container name
docker run -v $(pwd)/kmer/:/kmer -v ${vol_name}:/td_output --rm -it --entrypoint bash ${container_name}
