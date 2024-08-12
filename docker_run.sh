#!/bin/bash

# Function to display usage information
usage() {
  echo "Usage: $0 --vol <volume_name> --img <image_name> --cmd <command_to_exectute_inside_container>"
  exit 1
}

# Initialize variables
pd_vol=""
image_name=""
container_cmd=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --vol)
      pd_vol="$2"
      shift 2
      ;;
    --img)
      image_name="$2"
      shift 2
      ;;
    --cmd)
      shift
      container_cmd="$@"
      break
      ;;
    --help)
      usage
      ;;
    *)
      echo "Invalid option: $1" 1>&2
      usage
      ;;
  esac
done

# Check if the required flags were provided
if [ -z "$pd_vol" ] || [ -z "$image_name" ] || [ -z "$container_cmd" ]; then
  echo "Error: Arguments --vol, --img, and --cmd options are required."
  usage
fi

# Run the docker command with the specified volume and image name
docker run -v $(pwd)/kmer/:/kmer -v $(pwd)/${pd_vol}/:/td_output -it ${image_name} $container_cmd

