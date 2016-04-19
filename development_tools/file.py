

filename = '/'
filename = 'notes.md'


while(filename != None):
    
    filename=raw_input('Enter a file name : ')
    print 'filename = ',filename

    try:    

        with open(filename) as fp:
            print 'opened file'

    except IOError as err:
        print "Error reading the file {0}: {1}".format(filename, err)


        #filename = None

    print 'done'


