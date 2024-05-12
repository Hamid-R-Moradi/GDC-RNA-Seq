#!/usr/bin/bash
# version 1.0 Written by Hamid Reza Moradi, With a little help from internet fellas ;)
# for any questions or bugs email hr-moradi@alumnus.tums.ac.ir
# this script might not be the best way to implement this process but it will do the job. 
# provide it with the path to the directory that you have stored all the samples
# downloaded form gdc. ( You have extract them, this action is not implemented
# right now)

# Next line is the place to provide it with the path
path='path/to/samples/location'

# The results will be in the parent folder. due to some presumptions about the
# location of this folder in the next lines, changing it will break the script.
if ! [ -d $path/../samples ] ; then
	mkdir $path/../samples
fi

if find $path/../samples -mindepth 1 -maxdepth 1 | read ; then
	echo "directory it's not empty"
else
	for folder in $path/* ; do
		file=$(echo $folder/* | cut -d '/' -f 11)
		echo moving $file to samples
		cp $folder/* $path/../samples
	done

	rm $path/../samples/annotations.txt # removing the annotation file since it doesn't state anything worth while
fi

touch test
awk -f $path/../../code/countM.awk $path/../samples/* > test

# In order to remove the unnecessey prefixes in column names we define another variable to
# correct the forward slashes and fix the column names by it

path2correct=$(echo $path | sed 's/\//\\\//g;s/$/samples\\\//' )
sed "s/$path2correct//g" $path/../test > countM.tsv
rm test
