'''
execute_command.py

The purpose of this script is to make the matrix/map files determined by 
make_matrices.py. It takes in a single command parameter and executes the 
steps of that command in sequence. 

***This needs to be run on Python/3.9 as there's an incompatability with ***
***other versions of python.                                             ***
'''

#import modules that will be needed
import os
import argparse


###############################################################################
###########################  COLLECT AND CHECK ARGUMENTS  #####################
###############################################################################

## GLOBAL VARS

## MAIN ##

def main(argv): 
    #Create argument parser
    parser = argparse.ArgumentParser(description='Execute encoded commands.')
    parser.add_argument('-c', '--command_string', required=True, help='Encoded command string')
    args = parser.parse_args()

    #Call driver function
    driver(args.command_string)

###############################################################################
#############################  DRIVER  ########################################
###############################################################################

## drive the script ##
## ONE-TIME CALL -- called by main
def driver(command_string):
    
    #Decode triple underscore and carrots
    command_string = command_string.replace("___", " ")
    command_string = command_string.replace("^^^", ";")
    #Split command
    command_list = command_string.split(";")
    #Loop over the command_list and submit the commands
    for command in command_list:
        os.system(command)

###############################################################################
   
if __name__ == "__main__":
    main()

