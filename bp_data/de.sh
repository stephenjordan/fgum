#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# 1. Build the target
echo "Building src:density..."
bazel build src:density

# 2. Loop over the k and D pairs and run the binary
echo "Running density evaluations..."
for k in {3..7}; do
  for (( D=k+1; D<=8; D++ )); do
    echo "Executing with k=$k, D=$D..."
    ./bazel-bin/src/density --k "$k" --D "$D" --bins 2000 --iterations 2000 --precision 0.000001
  done
done

echo "All runs completed successfully!"