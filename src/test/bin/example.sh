#! /bin/sh
#
# Bash script running a test case
# (c) 2016, Ralf Quast
iniseed=27182
inistep=0.5
accuracy=0.000001
stop=2000
trace=10

./rq-edfit ${iniseed} 10 20 ${inistep} ${accuracy} ${stop} ${trace} < ./src/test/resources/example.in > example.html
