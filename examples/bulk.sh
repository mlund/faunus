#!/bin/bash
../yason.py bulk.yml | ../faunus
if [[ $? == 0 ]]; then
    python ../yason.py < out.json > out.yml
    rm -f out.json
fi
