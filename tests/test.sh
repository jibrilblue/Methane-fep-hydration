#!/bin/bash
mkdir -p /logs/verifier

pytest /tests/test_fep.py -v -rA --json-ctrf=/logs/verifier/ctrf.json
status=$?

if [ $status -eq 0 ]; then
    echo "1" > /logs/verifier/reward.txt
else
    echo "0" > /logs/verifier/reward.txt
fi

exit $status
