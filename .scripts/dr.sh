#!/bin/sh
docker run  \
    --user $(id -u):$(id -g) \
    --rm -t \
    -v `pwd`:/builds/app \
    jiaqiknight/dealii:vscode-arm64 \
    /bin/sh -c "cd /builds/app; $@"
