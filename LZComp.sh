
## ======================== Welcome to LZComp ========================= ##
# LZComp is a wrapper bash script to run lzcomp_ann.m
# LZComp runs several programs. one after the other, to either train or to
# predict the answers using a pre-trained set of parameters. In either
# case it can be run from here.
# LZComp calls both python and octave scripts (You will need both installed
# before continuing)

# LZComp uses the output from fitting 1, 2, and 3 component fits with
# lzifu (Ho et al 2016)
# The idea is to make it fast and easy to work out which is the best
# fit for each spaxel.
# The initial testing data was created by checking each spaxel by eye
## =================================================================== ##

## ============================= Author ============================== ##
# Elise Hampton
# PhD Candidate Australian National University
# Research Schol of Astronomy and Astrophysics
# LZComp written for the S7 collaboration (PI: A/Prof Michael Dopita) and
# the SAMI Survey (PI: Prof Scott Croom).
# Last update: July 2016
## =================================================================== ##


## =============================== log =============================== ##
# 12th Jan 2015 - Finch created
# 12th Jan 2015 - user input tested - good!
# 12th Jan 2015 - running CreateInput and theMachine from inside script
# 12th Jan 2015 - need to update input and machine to take in arguments
# 13th Jan 2015 - update user interface to read out slowly (for fun)
# 13th Jan 2015 - Now reads out strings in telletyper fashion
# 19th Jan 2015 - Sends train and test as command line arguments to
# CreateInput, TheMachine, and reimage
# 21st Jul 2016 - Renamed to better identify after publication to LZComp
# 27th Sept 2016 - The ANN script renamed to lzcomp_ann.m
## =================================================================== ##
clear

START=$(date +%s.%N)

foo="Hello!"
echo ${foo}

foo="What are we doing today?"
echo ${foo}

foo="Training (train)? Testing (test)? Both (traintest)? or just running (run)?"
echo ${foo}

echo ' '

while :
do
    read tmp

    if [ "$tmp" == "run" ]; then
	train=0
	testt=0
	break
    elif [ "$tmp" == "train" ]; then
	train=1
	testt=0
	break
    elif [ "$tmp" == "test" ]; then
	foo="We can't test without training at the moment! Setting to Both."
	echo ' '
	train=1
	testt=1
	break
    elif [ "$tmp" == "traintest" ]; then
	train=1
	testt=1
	break
    else
	echo "invalid input.... Try again?"
    fi
done

# This is for taking out spaxels that have no information or don't have all 3 component fits from LZIFU
# This can be changed for whatever you expect for your input
echo "Calling createNanmasks.py on galaxies.txt."
printf $foo
python2.7 input-codes/createNanmasks.py

# lzcomp_ann takes input in a text file CreateInput will take a LZIFU cube and turn it into the
# input needed for it
#call create input
foo="Calling CreateInput.py with train($train) and test($testt)."
echo ${foo}
echo ' '
python2.7 input-codes/CreateInput.py ${train} ${testt}

#call the ANN - classifier with the train and test input
foo="Calling lzcomp_ann.m with train($train) and test($testt)."
echo ${foo}
echo ' '
START2=$(date +%s.%N) #take the time
/Applications/Octave.app/Contents/Resources/bin/octave --silent --eval "ANN/lzcomp_ann" ${train} ${testt}

# reimage turns the output, which is a text file like the input, back into the same shape and size of the IFU obs
#call reimage
foo="Calling reimage.py train($train) and test($testt)."
echo ${foo}
echo ' '
python2.7 output-codes/reimage.py ${train} ${testt}

END2=$(date +%s.%N)
DIFF2=$END2-$START2
echo ${DIFF2}

#call merge comp if not training or testing
foo="Calling merge_components.py"
echo ${foo}
echo ' '
python2.7 output-codes/merge_comp.py


foo="Task Complete. See you next time!"
echo ${foo}
echo ' '

END=$(date +%s.%N)
DIFF=$END-$START
echo ${DIFF}
