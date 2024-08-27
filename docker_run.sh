#!/bin/bash

# Function to display usage information
usage() {
  echo "Usage: $0 --img <image_name> --cmd <command_to_exectute_inside_container>"
  exit 1
}

# Initialize variables
pd_vol="docker_pd_output"
image_name=""
primer_cmd=""
args_flag_set=0

if [ ! -d "docker_pd_output/" ]; then \
    echo "$pd_vol local volume is created..."
		mkdir $pd_vol
    mkdir $pd_vol/logs/
fi



# Parse command line arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --img)
      if [ "$args_flag_set" -eq 1 ]; then
        echo "Error: --img must be provided before --cmd." 1>&2
        usage
      fi
      image_name="$2"
      shift 2
      ;;
    --cmd)
      if [ -z "$pd_vol" ] || [ -z "$image_name" ]; then
        echo "Error: --cmd must come after --img." 1>&2
        usage
      fi
      args_flag_set=1
      shift
      primer_cmd="$@"
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
if [ -z "$image_name" ] || [ -z "$primer_cmd" ]; then
  echo "Error: Arguments --img and --cmd options are required."
  usage
fi

# Run the docker command with the specified volume and image name
docker run --rm -v $(pwd)/kmer/:/kmer \
           -v $(pwd)/${pd_vol}/:/td_output \
           -v $(pwd)/${pd_vol}/logs/:/logs \
           --user $(id -u):$(id -g) \
           -it ${image_name} $primer_cmd

echo "Primer Designer output is generated in $pd_vol local volume."
echo "Contents in $pd_vol:"
ls -1 $pd_vol
