#!/bin/bash

shopt -s extglob

#g++ -O2 A069675_tester.cpp -fopenmp -lgmp -lgmpxx --std=c++11
#nice -n 15 ./a.out &
#nice -n 10 python silly.py &
nice -n 15 ./ab_tester &
PID=$!

trap "echo \"killing PID($PID)\" && [ -n \"$PID\" ] && kill -SIGTERM $PID " EXIT

echo "PID: $PID"
echo -e "actions [irax]:\n"


#gdbus monitor -e -d com.canonical.Unity -o /com/canonical/Unity/Session


while read user_input
do
  case $user_input in
    [ir]?([1-9]*([0-9])s))
      letter=${user_input:0:1}
      case "$letter" in
        i)
          message="[I]nterupting"
          signal="-SIGSTOP"
          ;;
        r)
          message="[R]esuming"
          signal="-SIGCONT"
          ;;
      esac

      if [ ! -z "${user_input:1}" ]; then
        seconds=${user_input:1:-1}
        echo "$(date +\"%H:%M\"): $message in $seconds seconds"
        [ -n "$PID" ] && sleep $seconds && echo "$message at $(date +\"%H:%M\"))" && kill $signal $PID &
      else
        echo "$message"
        [ -n "$PID" ] && kill $signal $PID
      fi
      ;;

    [ax])
      echo "E[x]iting/[A]bort :("
      exit
      ;;
  esac

  echo
done
