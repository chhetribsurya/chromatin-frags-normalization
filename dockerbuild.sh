#!/usr/bin/bash

#This command builds a Docker image from the current directory’s Dockerfile, 
#targeting both Intel/AMD and ARM64 platforms, 
#tags it with both v1.3 and latest, and 
#then pushes those multi-architecture images to the Docker registry 
#under chhetribsurya/chromatin-frags-norm

docker buildx build --platform linux/amd64,linux/arm64 \
  -t chhetribsurya/chromatin-frags-norm:v1.3 \
  -t chhetribsurya/chromatin-frags-norm:latest \
  --push .
