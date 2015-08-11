#!/bin/bash
grep "!    total energy" $1 | tail -n 1 | awk '{print $5}'
