# Colored Reset
COff='\033[0m'       
# Regular Colors
Green='\033[0;32m' #OK message
Red='\033[0;31m'
Blue='\033[0;34m'
BRed='\033[1;31m'
BGreen='\033[1;32m'
BBlue='\033[1;34m'

# $Color Message $COff
if [ "$(ls -A output/)" ];
then rm output/*.* ;
/bin/echo -e "$BGreen Files removed from folder output $COff";
else /bin/echo -e "$BRed Folder output is empty $COff" ;
fi

if [ "$(ls -A dump/)" ];
then rm dump/*.* ;
/bin/echo -e "$BGreen Files removed from folder dump $COff";
else /bin/echo "$BRed Folder dump is empty $COff" ;
fi

