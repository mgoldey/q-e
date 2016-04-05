#!/bin/bash
grep "E field correction" $1 | tail -n 1 | awk '{print $5}'
