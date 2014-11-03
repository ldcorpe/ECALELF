#!/bin/bash 

# temp folder to hold intermediary files
mkdir temp_and_json/
# fill default file in case no datasets are specified...
echo "/DoubleElectron/Run2012A-15Apr2014-v2/AOD" >  temp_and_json/temp.list
#defaults for options.
OUTPUTNAME="and.json"
DATASETLIST="temp_and_json/temp.list"
NOMINALJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt"

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then 
echo "[HELP] About: this tool is used to create a json file which is the logical AND (&&) of the lumisections available for a list of datasets in the DAS and some nominal good lumisection json. Used in the clustering validation to counter error expereinced if part of the dataset was not reprocessed."
echo "[HELP] Possible options"
echo "[HELP]   -d|--datasetlist "
echo "[HELP]     text file containing the list of datasets to use to create the and.json with the nominal json."
echo "[HELP]     dataset names should be in the format /DoubleElectron/Run2012A-15Apr2014-v2/AOD"
echo "[HELP]     DEFAULT: $DATASETLIST (created at run time, contains example dataset from prev line) "
echo "[HELP]   -n|--nominaljson "
echo "[HELP]     The nominal json file to be used as the list of a priori good lumis. "
echo "[HELP]     DEFAULT: $NOMINALJSON "
echo "[HELP]     (should be this one for any run I study, I think) "
echo "[HELP]   -o|--outputname "
echo "[HELP]     desired name of output file "
echo "[HELP]     DEFAULT: $OUTPUTNAME"
echo "[HELP] Created by Louie D. Corpe 3.11.14"
exit
fi

#Retrieve options
while [[ $# > 1 ]]
do
key="$1"
shift

case $key in
-d|--datasetlist)
DATASETLIST="$1"
shift
;;
-n|--nominaljson)
NOMINALJSON="$1"
shift
;;
-o|--outputname)
OUTPUTNAME="$1"
shift
;;
*)
# unknown option
;;
esac
done

#counter to loop over files and give different names depending on dataset
COUNTER=0

echo "[INFO] using nominal json: "
echo "$NOMINALJSON"

#reading each line of DATASETLIST
while read line; do
DATASET=$line

echo "[INFO] dataset $COUNTER, $DATASET " 

#retrieve list of available lumis
echo "[INFO] Running das client"
das_client.py --format=json --query="run,lumi dataset=$DATASET" --limit=0  > temp_and_json/temp.json

#make edits to turn into usable format
echo "[INFO] making edits"
  sed -i 's/, {/\n{/g' temp_and_json/temp.json #split into lines for each run
	sed -i '1s/\[{/{/g' temp_and_json/temp.json  #remove "{" at start ..
	sed -i '$s/}\]/}/g' temp_and_json/temp.json  #.. and end
	sed -i 's/ //g' temp_and_json/temp.json  # remove spaces

	grep ",\"run\"" temp_and_json/temp.json > temp_and_json/temp2.json  #pick only cases where run/lumi are the wrong order
	sed -i 's/,\"run\"/, \"run\"/g' temp_and_json/temp2.json #add one space before run, effectively creates two columns
	awk ' { t = $1; $1 = $2; $2 = t; print; } ' temp_and_json/temp2.json > temp_and_json/temp3.json # exchange columns
	sed -i 's/\"run\":\[{\"run_number\"://g' temp_and_json/temp3.json #remove useless stuff
	sed -i 's/\"run\":{\"run_number\"://g' temp_and_json/temp3.json 
	sed -i 's/}\]}//g' temp_and_json/temp3.json 
	sed -i 's/{\"lumi\":\[{\"number\"://g' temp_and_json/temp3.json 
	sed -i 's/{\"lumi\":{\"number\"://g' temp_and_json/temp3.json 
	sed -i 's/}\],/,/g' temp_and_json/temp3.json 
	sed -i 's/}}/,/g' temp_and_json/temp3.json 
	sed -i 's/},/,/g' temp_and_json/temp3.json 
	sed -i 's/\, \[/\[/g' temp_and_json/temp3.json 
	sed -i 's/^[0-9]*/\"&\":/g' temp_and_json/temp3.json #put run number in " " followed by :
	sed -i 's/,/, /g' temp_and_json/temp3.json # put spaces back in

	grep ",\"lumi\"" temp_and_json/temp.json > temp_and_json/temp2.json #now only look at cases where run, lumi in correct order
	sed -i 's/{\"run\":\[{\"run_number\"://g' temp_and_json/temp2.json # remove useless stuff
	sed -i 's/{\"run\":{\"run_number\"://g' temp_and_json/temp2.json 
	sed -i 's/}\],//g' temp_and_json/temp2.json 
	sed -i 's/},//g' temp_and_json/temp2.json 
	sed -i 's/\"lumi\":\[{\"number\":/ /g' temp_and_json/temp2.json 
	sed -i 's/\"lumi\":{\"number\":/ /g' temp_and_json/temp2.json 
#sed -i 's/{\"lumi\":{\"number\"://g' temp_and_json/temp2.json 
	sed -i 's/}\]}/,/g' temp_and_json/temp2.json 
	sed -i 's/}}/,/g' temp_and_json/temp2.json 
	sed -i 's/^[0-9]*/\"&\":/g' temp_and_json/temp2.json # put run in " " and followed by :
	sed -i 's/,/, /g' temp_and_json/temp2.json # put spaces back 

echo "[INFO] done editing"
# new file to hold merged output of the two cases above	
FILENAME="temp_and_json/temp_file$COUNTER.json"
#actually merge
echo "[INFO] merging into $FILENAME"
cat temp_and_json/temp3.json temp_and_json/temp2.json > $FILENAME
echo "[INFO] files merged"

echo "[INFO] finished with dataset$COUNTER"

# increment counter to get distinct filenames
((COUNTER++))

#end of loop
done < $DATASETLIST

#merge all of the different dataset temp files
echo "[INFO] merging edited jsons from each step above"
cat temp_and_json/temp_file*.json >  temp_and_json/temp.json

#add a { and } at start and end
	sed -i '1s/^/{/g' temp_and_json/temp.json 
	sed -i '$s/, $/}/g' temp_and_json/temp.json

echo "[INFO]creating AND json with nominal json: "
echo $NOMINALJSON
echo "and saving to:"
echo $OUTPUTNAME

#do the "AND" bit wioth nominal file
compareJSON.py --and $NOMINALJSON temp_and_json/temp.json > $OUTPUTNAME

echo "[INFO] clearing up"
rm temp_and_json/*
