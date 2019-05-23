#!/bin/bash

#curl http://www.worldofnumbers.com/topic1.htm#conu49 > worldofnumbers_49.html
cat worldofnumbers_49.html | sed 's#<br>##g' | grep -oP '[0-9]{1,3}(\.[0-9]{3})+(?![0-9])' | tr -d '.' | sort -n > factors
