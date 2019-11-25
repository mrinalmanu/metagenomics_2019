# I don't want my classifier to crash. I want to be safe.
# 
# I will use some script to prase the taxonomy file to remove bugging characters. I used this bit of code from here:
# 
# https://gist.github.com/walterst/0a4d36dbb20c54eeb952
#
# This code has been converted into Python 3

from sys import argv

taxa_mapping = open(argv[1], "U")

for line in taxa_mapping:
    curr_line = line.strip()

    try:
        curr_line.decode('ascii')
        if "*" in curr_line:
            raise(ValueError)
        print(curr_line)
    except:
        fixed_line = ""
        for n in curr_line:
            if ord(n) < 128 and n != "*":
                fixed_line += n
        print(fixed_line)
        
        

