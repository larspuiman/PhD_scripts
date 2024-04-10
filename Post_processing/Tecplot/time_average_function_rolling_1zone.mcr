#!MC 700
#
# script to calculate the time-average of a specified variable, denoted by its variable number (see lists of variables in Tecplot). 
#     important to set the number of zones in several lines (here it is 1, this was hardcoded...)
# 
$!GETUSERINPUT |AVEVAR|
  INSTRUCTIONS = "Enter the number of the variable to average."
$!GETUSERINPUT |USRNAME|
  INSTRUCTIONS = "Enter the name of the resulting variable"
$!Varset |NTIMES| = (|NUMZONES|/1 -1)
$!Varset |zindex| = 2
$!Varset |zindex_2| = 1
$!Varset |counter| = 2
$! AlterData
 Equation = "{RAVG} =  V|AVEVAR|"

$!LOOP |NTIMES|
$!LINEMAP [|loop|] ASSIGN{ZONE = |zindex|}
$!LINEMAP [|loop|] ASSIGN{ZONE = |zindex_2|}
$! AlterData 
 [|zindex|] Equation = "{RAVG} = {RAVG}[|zindex_2|] + (V|AVEVAR|[|zindex_2|]-{RAVG}[|zindex_2|])/|counter|"
$!Varset |zindex| += 1
$!Varset |zindex_2| += 1
$!Varset |counter| += 1
$!ENDLOOP

$!GETVARNUMBYNAME |RAVG|
NAME = "RAVG"
$!RENAMEDATASETVAR
VAR = |RAVG|
NAME = "|USRNAME|"
# $!ENDMACROFUNCTION