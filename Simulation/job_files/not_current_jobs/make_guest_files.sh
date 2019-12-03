#!/bin/bash
for idx in *.sh; do
    cp $idx guest_$idx
done
