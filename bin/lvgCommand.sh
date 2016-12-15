#!/bin/bash
lvg -CR:o:oc:oe:oi -R:1 -f:q:q0:q1:q2:q3:q4:q5:q6:q7:q8:rs:T:l:t:g:P:g:C | awk 'BEGIN{FS="|"}{print $2}'
