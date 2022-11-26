#!/bin/bash

ls quantified/R* | sed -E -e 'p;s/R([[:digit:]])_/R0\1_/' | xargs -n2 mv
