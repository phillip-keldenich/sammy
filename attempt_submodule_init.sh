#!/bin/sh

if [ -d ".git" ]; then
  echo "git submodule update --init --recursive"
  git submodule update --init --recursive
fi
