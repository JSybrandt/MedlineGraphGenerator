#!/bin/bash
lvg -f:q:q0:q1:q2:q3:q8:rs:T:l:t:P:g:C | awk 'BEGIN{FS="|"}{print $2}'
