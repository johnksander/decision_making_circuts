#!/usr/bin/env bash

#for getting the best stimulus values in stim optimization search

JNAME=$1 #sim name 
Rdir="../../Results/$JNAME/" #results directory relative to where this file should be...
Rfn=output_log_	#look for these	logs in	the results
Fout=best_stim_vals.txt
Fout=$Rdir$Fout
echo "$Fout"

for idx in $Rdir$Rfn*txt; do
     targ=$(grep target "$idx" | head -1) #take the first line mentioning target value
     targ=$(echo "$targ" | grep -o "for.*target") #text specifiying the target
     targ=$(echo "$targ" | grep -Eo "[0-9]+[\.][0-9]+[s]") #target time
     
     echo "log file: $idx" >> $Fout
     echo "target duration = $targ" >> $Fout
     
     #now get the results output
     header=$(grep -Eo "\|.*\|.*\|" "$idx" | head -1)
     results=$(grep -Eo "\|.*\|.*\|" "$idx" | grep Hz)
  
     #add to file
     echo "$header" >> $Fout
     echo "$results" >> $Fout
     
     #strip the units from the target value
     targ=$(echo "$targ" | sed 's/s//;')
     #from the result text, grab outcome values
     Rvals=$(echo "$results" | grep -Eo "[0-9]+\.[0-9]+[s]" | sed 's/s//')

     curr_best=$(echo "$Rvals" | head -1)
     best_ln=1
     curr_ln=0
     for v in $Rvals; do
         curr_ln=$((curr_ln+1))
         if (( $(echo "sqrt(($targ-$v)^2) < sqrt(($targ-$curr_best)^2)" |bc -l) ));
         then
             curr_best=$(echo "$v") 
             best_ln=$(echo "$curr_ln")
         fi
     done
     
     best_info=$(echo "$results" | sed -n "$best_ln p")
     echo "    duration closest to target: $curr_best" >> $Fout
     echo "    $best_info" >> $Fout
     
     
     
     #add two newlines for sperating results
     echo "" >> $Fout
     echo "" >> $Fout 
done



