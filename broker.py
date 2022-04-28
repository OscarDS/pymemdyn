# This is a lame broker (or message dispatcher). When Gromacs enters a run, 
# it should choose a broker from here and dispatch messages through it.
#
# Depending on the broker, the messages may be just printed or something else

class Printing(object):
    def __init__(self):
        """
        This is a proxy to put a message directly to the stdout through
        *print* command
        """
        pass

    def dispatch(self, msg):
        """
        Simply print the msg passed
        """
        print (msg)