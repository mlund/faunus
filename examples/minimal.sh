#!/bin/bash

echo -e "1\n\nA 0.0 0.0 0.0" > A.xyz
echo -e "1\n\nB 0.0 0.0 0.0" > B.xyz

../yason.py minimal.yml | ../faunus
if [[ $? == 0 ]]; then
    python ../yason.py < out.json > out.yml
    rm -f out.json
fi
