try:
    # check to see if file is readable
    with open(filename) as tempFile:

except Exception as e:
    print e
    # here you can modify the error message to your liking


OR


try:
    with open(filename) as fp:

except IOError as err:
    print "Error reading the file {0}: {1}".format(filename, err)
    break


with open(filename) ... could be better as you don't have to explicitly close the file


