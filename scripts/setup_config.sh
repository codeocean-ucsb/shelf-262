#!/bin/bash

cp templates/*.template .

for f in *.template; do
    mv "$f" "${f%.template}"
done

echo "Configuration files created."
echo "Edit cppdefs.opt and *.opt before compiling."
