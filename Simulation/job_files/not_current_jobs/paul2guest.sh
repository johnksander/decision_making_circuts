#!/bin/bash

#this just replaces the stuff in guest*sh 
#job files that have paul-compute specifications
#with guest-compute info. Like... if you're making
#guest-compute job files from paul-compute ones etc.

sed -i 's/qos=medium/qos=low/;' guest*.sh
sed -i 's/paul-compute,neuro-compute,//;' guest*.sh
