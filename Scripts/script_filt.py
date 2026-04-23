#!/usr/bin/python

f = open("Mouse_hitdata.txt_filt3","r")

prevLine = ""

for line in f:
  #if line.startswith("QUERY"):
    #if prevLine.startswith("QUERY"):
    
    #else:
     # print(prevLine)
  if not line.startswith("Q#"):
    
    if prevLine.startswith("Q#"):
    
      print(prevLine)
    #  print(line)
    print(line)
  prevLine = line
