#!/bin/bash
lvg -f:0:C:P:q0:q1:q2:rs:g:T:t:u | awk 'BEGIN{FS="|"}{print $2}'
