#!/bin/sh
set -e


if [ -z "$GHCR_USERNAME" ]; then
  echo "GHCR_USERNAME not set" >&2
  exit 1
fi

if [ -z "$GHCR_TOKEN" ]; then
  echo "GHCR_TOKEN not set" >&2
  exit 1
fi

echo "$GHCR_TOKEN" | docker login ghcr.io -u "$GHCR_USERNAME" --password-stdin

docker build -t ghcr.io/${GHCR_USERNAME}/mde-r:latest .
docker push ghcr.io/${GHCR_USERNAME}/mde-r:latest
