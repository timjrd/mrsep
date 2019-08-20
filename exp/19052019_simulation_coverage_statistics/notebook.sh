#! /usr/bin/env nix-shell
#! nix-shell --pure -i bash
test -d results || tar -xf results.tar.xz
jupyter-notebook
